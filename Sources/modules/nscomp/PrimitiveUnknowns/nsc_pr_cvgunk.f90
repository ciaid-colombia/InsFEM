subroutine nsc_pr_cvgunk(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    compressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_NSCompressiblePrimitive
   use MPI
   use def_parame
   implicit none
   class(NSCompressiblePrimitiveProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip), save       :: ipass=0
   integer(ip)             :: npoinLocal,ndime
   real(rp)                :: prnsc,vrnsc,trnsc
   real(rp)                :: rpnsc,rvnsc,rtnsc
   integer                 :: ierr
   
   select case(itask)
      
   !Check convergence of the inner iterations, always in the L2 norm:
   !|| U(n,i,j) - U(n,i,j-1)|| / ||U(n,i,j)||
   case(1)

      call a%InnerResiduals(prnsc,vrnsc,trnsc)
     
      if (a%MPIrank == a%MPIroot) then
         if((((prnsc<a%cotol).and.(vrnsc<a%cotol).and.(trnsc<a%cotol)).or.(a%itera>=a%maxit))) a%kfl_goite = 0       

         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
         end if

         write(a%lun_conve,101) a%istep,a%itera,a%ctime,prnsc,vrnsc,trnsc

      endif

      if (a%kfl_flush == 1) call flush(a%lun_conve)
      
      CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| U(n,i,*) - U(n,i-1,*)|| / ||U(n,i,*)||
   case(2)

      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPI(a%kfl_normc,ndime*npoinLocal,a%veloc(1,1,1),a%veloc(1,1,2),a%resid,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

      !Outter Residuals is used by the coupling between modules (at master level)
       
   !Check residual of the time evolution, always in the L2 norm:
   !|| U(n,*,*) - U(n-1,*,*)|| / ||U(n,*,*)||
   case(3)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPI(two,npoinLocal,      a%press(1,1),  a%press(1,3),  rpnsc,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,ndime*npoinLocal,a%veloc(1,1,1),a%veloc(1,1,3),rvnsc,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npoinLocal,      a%tempe(1,1),  a%tempe(1,3),  rtnsc,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,rpnsc,rvnsc,rtnsc
         if(rpnsc.le.a%sstol) then
            if(rvnsc.le.a%sstol) then
               if(rtnsc.le.a%sstol) then
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
       &      '       Current','      Pressure','      Velocity','      Temperature',/,&
       & '$ ','       step','  iteration',&
       &      '          time','      residual','      residual','      residual')
101 format(4x,i9,2x,i9,20(2x,e12.6))
102 format('$ >>>  PROBLEM IS STATIONARY AT TIME STEP ',i5)

end subroutine nsc_pr_cvgunk
