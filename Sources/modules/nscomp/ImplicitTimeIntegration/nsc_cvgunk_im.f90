subroutine nsc_cvgunk_im(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    compressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_NSCompressibleImplicit
   use MPI
   use def_parame
   implicit none
   class(NSCompressibleImplicitProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip), save       :: ipass=0
   integer(ip)             :: npoinLocal,ndime
   real(rp)                :: drnsc,mrnsc,ernsc
   real(rp)                :: rdnsc,rmnsc,rensc
   integer                 :: ierr
   
   select case(itask)
      
   !Check convergence of the inner iterations, always in the L2 norm:
   !|| U(n,i,j) - U(n,i,j-1)|| / ||U(n,i,j)||
   case(1)

      call a%InnerResiduals(drnsc,mrnsc,ernsc)
     
      if (a%MPIrank == a%MPIroot) then
         if((((drnsc<a%cotol).and.(mrnsc<a%cotol).and.(ernsc<a%cotol)).or.(a%itera>=a%maxit))) a%kfl_goite = 0       

         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
         end if

         write(a%lun_conve,101) a%istep,a%itera,a%ctime,drnsc,mrnsc,ernsc

      endif

      if (a%kfl_flush == 1) call flush(a%lun_conve)
      
      CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| U(n,i,*) - U(n,i-1,*)|| / ||U(n,i,*)||
   case(2)

      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPI(a%kfl_normc,ndime*npoinLocal,a%momen(1,1,1),a%momen(1,1,2),a%resid,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
   
      !Outter Residuals is used by the coupling between modules (at master level)


   case(3)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPI(two,npoinLocal,      a%densf(1,1),  a%densf(1,3),  rdnsc,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,ndime*npoinLocal,a%momen(1,1,1),a%momen(1,1,3),rmnsc,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npoinLocal,      a%energ(1,1),  a%energ(1,3),  rensc,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,rdnsc,rmnsc,rensc
         if(rdnsc.le.a%sstol) then
            if(rmnsc.le.a%sstol) then
               if(rensc.le.a%sstol) then
                  a%kfl_stead = 1
                  write(a%lun_conve,102) a%istep
               end if
            end if
         end if
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   end select       

   !Formats
100 format('$ ','       Time','      Inner',&
       &      '       Current','      Density','      Momentum','      Energy',/,&
       & '$ ','       step','  iteration',&
       &      '          time','      residual','      residual','      residual')
101 format(4x,i9,2x,i9,20(2x,e12.6))
102 format('$ >>>  PROBLEM IS STATIONARY AT TIME STEP ',i5)

end subroutine nsc_cvgunk_im
