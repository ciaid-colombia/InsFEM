subroutine ale_solite(a)
   !-----------------------------------------------------------------------
   !> This is the main routine of the linear system of mesh displacements.
   !! The ALE flag is set to 0 so the displacement field is computed in the
   !! original coordinates. The currentbvess variable allows the decoupling 
   !! of the system for each dimension of the mesh. Finally, the ALE flag is
   !! reset to 1 so the computed displacements can be added to the coordinates
   !! when loading the mesh, as well as to update the mass matrix and the 
   !! exterior normal.
   !-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_Element
   use Mod_Alemov
   implicit none
   class(AlemovProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   type(TimeIntegratorDt1) :: Integrator
   
   integer(ip)  :: kfl_perio,idofn,ipoin,npoin,ielem,nelem,inode,ndime,currentbvess,nsteps,i,icomp
   real(rp)     :: LHSdtinv
   real(rp), allocatable :: valRHS(:,:),auxDisplacement(:,:,:)
   integer(ip)  :: kfl_HangingNodes
   
   logical :: AreFoldedALE 
   
   interface
      subroutine ale_BuildLinearSystem(a,currentbvess)
         use typre
         import AlemovProblem
         implicit none
         integer(ip) :: currentbvess
	 class(AlemovProblem) :: a
      end subroutine
   end interface

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   
   
   call a%Memor%alloc(ndime,npoin,1_ip,auxDisplacement,'auxDisplacement','ale_solite')
   
   !We compute the ALE movement with the undeformed mesh
   call a%Mesh%SetALE(0_ip)
   do currentbvess=1,a%ndofbc
      
      a%itera = a%itera + 1
      if (a%MPIrank == a%MPIroot) then
         if(a%itera==1) write(a%lun_solve,100) a%istep
         write(a%lun_solve,101) a%itera      
      endif
   
      call a%Timer%BuildLinearSystem%Tic
      !Construct the system matrix and right-hand-side.
      call ale_BuildLinearSystem(a,currentbvess)
      call a%Timer%BuildLinearSystem%Toc
   
      !Solve the algebraic system.
      call a%Timer%SolveLinearSystem%Tic
      call a%LinearSystem%Solve(a%unkno)
      
      !Ghost communicates
      call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%unkno)
   
      !HangingNodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(a%ndofn,a%unkno)
   
      !Periodic boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(a%ndofn,a%unkno)
   
      call a%Timer%SolveLinearSystem%Toc
      
      call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         auxDisplacement(currentbvess,ipoin,1) = a%unkno(1,ipoin)
      enddo
      
      
   end do

   
      !Formats. 
      100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
      101 format('------------------------------------------------------------', &
            /,'   INNER ITERATION NUMBER: ',i5)
   
   call a%Mesh%SetALE(1_ip)
   
   !At this point, even if the elements are folded, we proceed, although we will project the results if necessary
   !Update displacement
   do i = 1,nsteps-1 
      a%Displacement(:,:,nsteps+1-i) = a%Displacement(:,:,nsteps-i)   
   enddo
   a%Displacement(:,:,1) = auxDisplacement(:,:,1)

   !Compute Velocity
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Memor%alloc(ndime,npoin,valRHS,'valRHS','ale_solite')
   call Integrator%GetRHS(ndime*npoin,a%Displacement(:,:,2:nsteps),valRHS)

   if (a%istep>nsteps-1 .or. a%kfl_inter==0) then
      a%Velocity(:,:,1) = LHSDtinv*a%Displacement(:,:,1) - valRHS*a%dtinv
   endif

   call a%Memor%dealloc(ndime,npoin,valRHS,'valRHS','ale_solite')
   call a%Mesh%SetVelocities(a%Velocity)
   
!    call a%FilePostpr%postpr(a%Velocity(:,:,1)-a%NSVeloc,'LevelSETAdvection',a%istep,a%ctime,a%Mesh)
!    call a%FilePostpr%postpr(a%Displacement(:,:,1),'ALEDisplacement',a%istep,a%ctime,a%Mesh)
   
   
   !Remesh/Project only when folded elements
   if (a%kfl_RemeshingCriteria == 0) then
      !Check if the computed displacement leads to element folding
      call a%Mesh%SetDisplacements(auxDisplacement)
      call a%Mesh%ComputeCheckFoldedALE
      call a%Mesh%GetCheckFoldedALE(AreFoldedALE)
      call a%Mesh%SetDisplacements(a%Displacement)
      
      if (AreFoldedALE .eqv. .true.) then
         a%DoRemesh = 1
      else
         a%DoRemesh = 0
      endif
      
   !Remesh/project at each time step
   elseif (a%kfl_RemeshingCriteria == 2) then
      a%DoRemesh = 1
   
   
   endif
   
   
  
   
   if (a%DoRemesh == 0) then
   
      !Dealloc and recompute ExtnorLpoty and Vmass
      call a%Mesh%DeallocExnorLpoty
      call a%Mesh%DeallocVmass
      
      !Recompute
      call a%Mesh%ComputeVmass
      call a%Mesh%ExtnorLpoty
   else
      
      !We need to revert the mesh displacement to the previous step
      !By now we do nothing, it will be done in the mesh projection step
   
   
   endif
   
   
   call a%Memor%dealloc(ndime,npoin,1_ip,auxDisplacement,'auxDisplacement','ale_solite')
   
   
   
   

end subroutine ale_solite
