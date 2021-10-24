subroutine sld_cvgunk(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    Cauchy equations
   !-----------------------------------------------------------------------
   use MPI
   use typre
   use def_parame
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip), save       :: ipass=0, ipasscp=0
   integer(ip)             :: npoinLocal,ndime,idime,ipoin
   real(rp)                :: resld
   integer                 :: ierr
   logical                 :: isALE

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoinLocal(npoinLocal)

   select case(itask)

   !Check convergence of the inner iterations, always in the L2 norm:
   !|| u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
   case(1)
      
       call a%InnerResiduals(resld)

       if (a%MPIrank == a%MPIroot) then
         if(((resld<a%cotol).or.(a%itera>=a%maxit))) then

             a%kfl_goite = 0      

             !TODO:improve this, cleaner code
             !I'm using this as a force convergence loop
             if(a%sld_type== 'NONLI' ) then
                 if(a%force_factor < 0.9999999999_rp) then

                     a%kfl_goite = 1  
                     a%itera=0    
                     a%force_factor= a%force_factor + a%delta_force

                     call a%Mesh%DeallocExnorLpoty
                     call a%Mesh%DeallocVmass

                     !Recompute
                     call a%Mesh%ComputeVmass
                     call a%Mesh%ExtnorLpoty

                     call a%SpecificBegite

                 end if
             end if
         end if

         if(ipass==0) then
             ipass=1
             write(a%lun_conve,100)
             write(a%lun_outph,200)
         end if

         write(a%lun_conve,101) a%istep,a%itera,a%ctime,resld
         
         ! Write in the log file.
         if (a%kfl_flush == 1) then
            call flush(a%lun_outph)
            call flush(a%lun_conve)
         endif
      endif
      CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
   case(2)

      call vecresMPI(a%kfl_normc,ndime*npoinLocal,a%disp(:,:,1),a%disp(:,:,2),a%resid,zesld,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      
   !Check residual of the time evolution, always in the L2 norm:
   !|| u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
   case(3)

      call vecresMPI(two,ndime*npoinLocal,a%disp(:,:,1),a%disp(:,:,3),resld,zesld,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,resld
         if(resld.le.a%sstol) then
            a%kfl_stead = 1
            write(a%lun_conve,102) a%istep
         end if
         call a%Mesh%GetALE(isALE)
         if (isALE) a%kfl_stead=0
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

      !Check residual of case coupling, for example FSI
   case(4)

      call vecresMPI(two,ndime*npoinLocal,a%disp(:,:,1),a%disp_cp(:,:),a%cpres,zesld,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

      if (a%MPIrank == a%MPIroot) then

          if(ipasscp==0) then
              ipasscp=1
              write(a%lun_cpconve,100)
              write(a%lun_outph,200)
          end if

          write(a%lun_cpconve,101) a%istep,a%cpiter,a%ctime,a%cpres

          if (a%kfl_flush == 1) then
              call flush(a%lun_cpconve)
          endif     

          if(a%cpres.le.a%cptol) then
              a%kfl_coupconv= .true.
          end if
      endif

      CALL MPI_BCAST(a%kfl_coupconv, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
      
      
   end select       

   !Formats
100 format('$ ','       Time','      Inner',&
       &      '       Current','      Displacement',/,&
       & '$ ','       step','  iteration',&
       &      '          time','      residual')
103 format('$ ','       Time','     Global','      Inner',' Non-Linear',&
       &      '       Current','  Displacement',/,&
       & '$ ','       step','  iteration','  iteration','  iteration',&
       &      '          time','  residual','      residual')

101 format(4x,i9,2x,i9,20(2x,e12.6))
106 format(4x,i9,2x,i9,20(2x,e12.6),20(2x,e12.6))
104 format(2x,4(2x,i9),3(2x,e12.6),4x,a)
102 format('$ >>>  DISPLACEMENTS STATIONARY AT TIME STEP ',i5)

200 format(///,&
       5x, '>>>  EVOLUTION OF NUMERICAL PARAMETERS:',///,&
       10x,a)
300 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','      Max Disp','      Arg Max1',/,&
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual','      residual')
301 format(4x,i9,2x,i9,2x,i9,2x,e12.6,20(2x,e12.6,2x,i6))

end subroutine sld_cvgunk

subroutine sld_InnerResiduals(a,resld)
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   real(rp) :: resld
   
   integer(ip) :: ndime, npoinLocal
   
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   call vecresMPIHeterogeneous(2_ip,ndime,ndime,npoinLocal,1,1,ndime,a%unkno(:,:),a%disp(:,:,1),resld,zesld,a%MPIcomm,a%MPIroot)
   !call vecresMPIHeterogeneous(2_ip,ndime,ndime,npoinLocal,1,1,ndime,a%unkno(:,:)+a%disp(:,:,1),a%disp(:,:,1),resld,zesld,a%MPIcomm,a%MPIroot)

end subroutine
