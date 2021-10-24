subroutine plcd_cvgunk(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    Cauchy equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_PLCD
   use MPI
   use def_parame
   use Mod_plcd_Stages
   implicit none
   class(PLCDProblem), target :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip), save       :: ipass=0
   integer(ip)             :: ndime,idime,ipoin
   real(rp)                :: va,vo,numer,denom,vnume,vdeno,pnume,pdeno,rinsi
   real(rp)                :: rwa(4),rwa2(4)
   integer                 :: ierr
   logical                 :: isALE
   
   real(rp) :: ForcesResidualNorm
   integer(ip) :: npoin
   
   interface
      subroutine plcd_ComputeForcesResidualNorm(a,ForcesResidualNorm)
         import
         implicit none
         class(PLCDProblem), target :: a
         real(rp) :: ForcesResidualNorm
      end subroutine
   end interface
   
   select case(itask)
      
   case(1)
      !Convergence in forces 
      
      call plcd_ComputeForcesResidualNorm(a,ForcesResidualNorm)
      a%ForcesResidualNorm = ForcesResidualNorm

      if (a%MPIrank == a%MPIroot) then
         if((   ( ForcesResidualNorm<a%css%IterationsTolerance )    .or.(a%itera>=a%css%MaximumNonLinearIterations))) then 
            a%kfl_goite = 0      
         endif
         
         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
            write(a%lun_outph,200)
            !write(a%lun_nolin,300)
         end if
         
         write(a%lun_conve,101) a%istep,a%itera,a%ctime,ForcesResidualNorm
         
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
!       call a%Mesh%GetNpoinLocal(npoinLocal)
!       call a%Mesh%GetNdime(ndime)
!       call vecresMPI(a%kfl_normc,ndime*npoinLocal,a%Displacement(1,1,1),a%Displacement(1,1,2),a%resid,zeplcd,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      
   !Check residual of the time evolution, always in the L2 norm:
   !|| u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
   case(3)
!       call a%Mesh%GetNpoinLocal(npoinLocal)
!       call a%Mesh%GetNdime(ndime)
!       call vecresMPI(two,ndime*npoinLocal,a%Displacement(1,1,1),a%Displacement(1,1,3),rinsi,zeplcd,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
! 
!       if (a%MPIrank == a%MPIroot) then
!          write(a%lun_conve,101) a%istep,0,a%ctime,rinsi
!          if(rinsi.le.a%sstol) then
!             a%kfl_stead = 1
!             write(a%lun_conve,102) a%istep
!          end if
!          call a%Mesh%GetALE(isALE)
!          if (isALE) a%kfl_stead=0
!       endif
!       CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   end select       

   !Formats
100 format('$ ','       Time','      Inner',&
       &      '       Current','      Force',/,&
       & '$ ','       step','  iteration',&
       &      '          time','      residual','      residual')
103 format('$ ','       Time','     Global','      Inner',' Non-Linear',&
       &      '       Current','      Displacement',/,&
       & '$ ','       step','  iteration','  iteration','  iteration',&
       &      '          time','      residual','      residual')

101 format(4x,i9,2x,i9,20(2x,e12.6))
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

end subroutine plcd_cvgunk

subroutine plcd_InnerResiduals(a,rinsi)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem) :: a
   real(rp) :: rinsi
   
   integer(ip) :: ndime, npoinLocal
   
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   call vecresMPIHeterogeneous(2_ip,ndime,ndime,npoinLocal,1,1,ndime,a%unkno,a%Displacement(1,1,1),rinsi,zeplcd,a%MPIcomm,a%MPIroot)
   
end subroutine

subroutine plcd_ComputeForcesResidualNorm(a,ForcesResidualNorm)
   use typre
   use Mod_PLCD
   implicit none  
   class(PLCDProblem), target :: a
   real(rp) :: ForcesResidualNorm, ExternalForcesNorm
   integer(ip) :: npoinlocal
   
   call a%Mesh%GetNpoinLocal(npoinlocal)
   call vecnorMPI(a%ResidualForcesVector,a%ndofn*npoinlocal,ForcesResidualNorm,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)
   call vecnorMPI(a%ExternalForcesVector,a%ndofn*npoinlocal,ExternalForcesNorm,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)
   ForcesResidualNorm = ForcesResidualNorm/ExternalForcesNorm
   
end subroutine
   
   


