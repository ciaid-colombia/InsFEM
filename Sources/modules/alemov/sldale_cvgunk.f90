subroutine sldale_cvgunk(a,itask)
   use MPI
   use typre
   use def_parame
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem) :: a
   integer(ip) :: itask

   real(rp)                :: riale,ridsp
   integer(ip)             :: npoinLocal,ndime,ipoin
   integer                 :: ierr

   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   
   select case(itask)
   
   !Inner iterations
   case(1)

       call vecresMPIHeterogeneous(2_ip,ndime,ndime,npoinLocal,1,1,ndime,&
               &a%unkno(:,:),a%Displacement(:,:,1),riale,0.0_rp,a%MPIcomm,a%MPIroot)

       call vecresMPIHeterogeneous(2_ip,ndime,ndime,npoinLocal,1,1,ndime,&
               &a%Velocity(:,:,1),a%Velocity(:,:,2),ridsp,0.0_rp,a%MPIcomm,a%MPIroot)


       !For ALE we dont converge the mesh, rather the attached module has to
       !converge
       a%kfl_goite = 0      

       if (a%MPIrank == a%MPIroot) then

           if(a%kfl_cvg) then
               a%kfl_cvg = .false.
               write(a%lun_conve,101)
           end if

           write(a%lun_conve,102) a%istep,a%itera,a%ctime,riale,ridsp

           ! Write in the log file.
           if (a%kfl_flush == 1) then
               call flush(a%lun_conve)
           endif
       endif

       CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

   !Check convergence of the outer iterations in the norm selected by the user:
   !|| u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
   !Outer iterations
   case(2)
      
   
   !End of step convergence
   case(3)

   !Coupling convergence
   case(4)

       call a%Mesh%GetNpoinLocal(npoinLocal)
       call a%Mesh%GetNdime(ndime)
       call vecnorGenericMPI(a%bres(:,:,1),npoinLocal,ndime,a%cpres,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)

       if (a%MPIrank == a%MPIroot) then

           if(a%kfl_cvgcp) then
               a%kfl_cvgcp = .false.
               write(a%lun_cpconve,100)
           endif

           write(a%lun_cpconve,102) a%istep,a%cpiter,a%ctime,a%cpres,a%sldale_dispRelax
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
          &      '       Current','  Displacement','    Relaxation',/,&
          & '$ ','       step','  iteration',&
          &      '          time','      residual','   Coefficient')
  101 format('$ ','       Time','      Inner',&
          &      '       Current','  Displacement','      Velocity',/,&
          & '$ ','       step','  iteration',&
          &      '          time','      residual','      Residual')
  102 format(4x,i9,2x,i9,20(2x,e12.6))

end subroutine
