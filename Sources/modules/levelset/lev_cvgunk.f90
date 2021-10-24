subroutine lev_cvgunk(a,itask)
   use typre
   use Mod_LevelSet
   use def_parame
   use MPI
   implicit none
   class(LevelSetProblem) :: a
   integer(ip) :: itask

   integer(ip), save       :: ipass=0
   real(rp)                :: ta,to,numer,denom,rilev
   integer(ip)             :: npoinLocal
   integer                 :: ierr
   
   select case(itask)
   
   !Check convergence of the inner iterations, always in the L2 norm:
   !|| L(n,i,j) - L(n,i,j-1)|| / ||L(n,i,j)||
   case(1)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call vecresMPI(a%kfl_normc,npoinLocal,a%unkno,a%level(1,1),rilev,zelev,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      
      if (a%MPIrank == a%MPIroot) then
         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
         end if
         write(a%lun_conve,101) a%istep,a%itera,a%ctime,rilev
         
         ! Write in the log file.
         if (a%kfl_flush == 1) call flush(a%lun_conve)
      endif
      
      if((rilev<a%cotol).or.(a%itera>=a%maxit)) a%kfl_goite = 0
      
      
      CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| L(n,i,*) - L(n,i-1,*)|| / ||L(n,i,*)||
   case(2)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call vecresMPI(a%kfl_normc,npoinLocal,a%level(1,1),a%level(1,2),rilev,zelev,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)      
   
   !Check residual of the time evolution, always in the L2 norm:
   !|| L(n,*,*) - L(n-1,*,*)|| / ||L(n,*,*)||
   case(3)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call vecresMPI(two,npoinLocal,a%level(1,1),a%level(1,3),rilev,zelev,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)      
      
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,rilev
         if(rilev<=a%sstol.and.a%kfl_conbc==1) then
            a%kfl_stead = 1
            write(a%lun_conve,102) a%istep
         end if
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   end select
   !
   ! Formats
   !
   100 format('$ ','       Time','      Inner','       Current','   LevelSet',/,&
            & '$ ','       step','  iteration','          time','      residual')
   101 format(4x,i9,2x,i9,10(2x,e12.6))
   102 format('$ >>>  LEVELSET IS STATIONARY AT TIME STEP ',i5)

end subroutine