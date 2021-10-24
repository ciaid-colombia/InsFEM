subroutine tem_cvgunk(a,itask)
   use typre
   use Mod_Temperature
   use def_parame
   use MPI
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: itask

   integer(ip), save       :: ipass=0
   real(rp)                :: ta,to,numer,denom,ritem,tenew 
   integer(ip)             :: npoinLocal
   integer                 :: ierr
   
   select case(itask)
   
   !Check convergence of the inner iterations, always in the L2 norm:
   !|| T(n,i,j) - T(n,i,j-1)|| / ||T(n,i,j)||
   case(1)
      call a%Mesh%GetNpoinLocal(npoinLocal)

      if (npoinLocal > 0) then
         call vecresMPI(2,npoinLocal,a%unkno,a%tempe(1,1),ritem,zetem,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      else   
         call vecresMPI(2,npoinLocal,a%unkno,0.0_rp,ritem,zetem,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      endif
      
      if (a%MPIrank == a%MPIroot) then
         if((ritem<a%cotol).or.(a%itera>=a%maxit)) a%kfl_goite = 0
         
         !Some postprocessing
         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
         end if
         write(a%lun_conve,101) a%istep,a%itera,a%ctime,ritem
         ! Write in the log file.
         if (a%kfl_flush ==1) call flush(a%lun_conve)
      endif
      CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| T(n,i,*) - T(n,i-1,*)|| / ||T(n,i,*)||
   case(2)
      call a%Mesh%GetNpoinLocal(npoinLocal)

      if (npoinLocal > 0) then
         call vecresMPI(a%kfl_normc,npoinLocal,a%tempe(1,1),a%tempe(1,2),a%resid,zetem,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      else
         call vecresMPI(a%kfl_normc,npoinLocal,0_rp,0_rp,a%resid,zetem,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      endif
   
   !Check residual of the time evolution, always in the L2 norm:
   !|| T(n,*,*) - T(n-1,*,*)|| / ||T(n,*,*)||
   case(3)
      call a%Mesh%GetNpoinLocal(npoinLocal)

      if (npoinLocal > 0) then
         call vecresMPI(two,npoinLocal,a%tempe(1,1),a%tempe(1,3),ritem,zetem,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      else
         call vecresMPI(two,npoinLocal,0_rp,0_rp,ritem,zetem,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      endif

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,ritem
         if(ritem<=a%sstol.and.a%kfl_conbc==1) then
            a%kfl_stead = 1
            write(a%lun_conve,102) a%istep
         end if
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   end select
   !
   ! Formats
   !
   100 format('$ ','       Time','      Inner','       Current','   Temperature',/,&
            & '$ ','       step','  iteration','          time','      residual')
   101 format(4x,i9,2x,i9,10(2x,e12.6))
   102 format('$ >>>  TEMPERATURE IS STATIONARY AT TIME STEP ',i5)

end subroutine
