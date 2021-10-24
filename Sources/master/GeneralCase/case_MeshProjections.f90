module Mod_caseMeshProjections
   use typre
   use Mod_GeneralCase
   use Mod_DistributedContainer
   use Mod_PhysicalProblem
   use Mod_DriverInterface
   use Mod_caseVariables
   use Mod_PhysicalProblemDriver
   use Mod_DC_Driver
   implicit none
   
   type(caseVariables), pointer :: c => NULL()
   
contains
   subroutine LoopMeshProjections1(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%MeshProjections(c,1_ip)
   end subroutine
   
   subroutine LoopMeshProjections2(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%MeshProjections(c,2_ip)
   end subroutine   
   
   subroutine LoopMeshAdvection1(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%MeshAdvections(c,1_ip)
   end subroutine
   
   subroutine LoopMeshAdvection2(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%MeshAdvections(c,2_ip)
   end subroutine 
   
end module
   
   

subroutine case_MeshProjections(a)
   use MPI
   use Mod_caseMeshProjections
   implicit none
   class(GeneralCase), target :: a
   
   logical :: AreFoldedALE
   real(rp), pointer :: coord(:,:) => NULL()
   
   integer(ip) :: DoRemesh
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
   class(PhysicalProblem), pointer :: myPhysical => NULL()
   
   
   type(domainVariables), pointer :: d => NULL()
   type(masterVariables), pointer :: m => NULL()
   
   character(DCHashCharSize) :: DoFirstKeys(5)
   character(5) :: auxKey
   integer :: ierr
   real(rp), pointer :: meshve(:,:,:) => NULL()
   
   logical :: isALE
   
   c => a%caseVars
   m => a%caseVars%masterVars
   d => a%caseVars%domainVars
   
   call d%cpu_MeshProjections%Tic
   
   !It makes no sense if no ALE
   call a%DriverList%GetFromKey('ALEPR',myDC)
   if (.not. associated(myDC)) return
   
   !Check if we need to ReMesh, project etc
   !The ALE physical problem knows it
   call ExtractDriver(myDC,myDriver)
   !Extract the physical problem from the physicalDriver
   call myDriver%GetRemeshingCriteria(DoRemesh)
   
   if (DoRemesh == 1) then
      
      !Fixedd%MeshALE
      if (m%kfl_ReMeshingStrategy == 2) then
         !If Fixed d%Mesh ALE, we need to project the results onto the background d%Mesh
         !Setup the FMALE interpolator
         call d%Mesh%GetCoord(coord)
         call d%FMALEInterpolator%SetOutputFiles(m%lun_memor,m%lun_outpu)
         call d%FMALEInterpolator%Initialize_Interp(d%Mesh,coord)
      
         !Initializations.
         call a%DriverList%LoopList(LoopMeshProjections1)
      
         !ALE must be the first one of the second pass, so that the displacements are set to zero
         !Also ALE takes care of telling the d%Mesh to recompute the normal and vmass
         DoFirstKeys(1) = 'ALEPR'
         call a%DriverList%LoopList(LoopMeshProjections2,DoFirstKeys(1:1))
         
         !3rd loop (only required for ALEmov, resetting updbcs
         !call a%DriverList%LoopList(LoopMeshProjections3)
      
         call d%FMALEInterpolator%Finalize
         
      elseif (m%kfl_ReMeshingStrategy == 3) then
         !In this case, instead of doing a projection we are going to advect all the fields
         !This is only valid if it is done at each time step
         !It would equivalent to separating advection from the rest of the equations if ALEveloc = veloc
         !First we need to build an "Advector" object
         !Secondly, we need to advect all the fields in each problem
         call d%Mesh%GetAle(isALE)
         call d%Mesh%SetALE(0)
         
         call d%FMALEAdvector%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         call d%FMALEAdvector%SetMesh(d%Mesh)
         call d%FMALEAdvector%SetParallelLibrary(m%ParallelLibrary)
         call d%FMALEAdvector%SetMemor(d%Mesh%Memor)
         call d%FMALEAdvector%Initialize
         
         call d%Mesh%GetMeshVeloc(meshve)
         call d%FMALEAdvector%BuildMatrix(m%dtime,meshve(:,:,1))
         
         
         !ALE must be the last one here so that the advection velocity is the last one to be advected
         DoFirstKeys(1) = 'NSTIN'
         DoFirstKeys(2) = 'LEVSE'
         
         !Initializations.
         call a%DriverList%LoopList(LoopMeshAdvection1,DoFirstKeys(1:2))
         
         !ALE must be the first one of the second pass, so that the displacements are set to zero
         !Also ALE takes care of telling the d%Mesh to recompute the normal and vmass
         DoFirstKeys(1) = 'ALEPR'
         call a%DriverList%LoopList(LoopMeshAdvection2,DoFirstKeys(1:1))
         
         call d%FMALEAdvector%Finalize
      
         if (isALE .eqv. .true.) then
            call d%Mesh%SetALE(1_ip)
         else
            call d%Mesh%SetALE(0_ip)
         endif
      else
         if (m%MPIrank == m%MPIroot) then
             call myDC%getKey(auxKey)
             write(*,*) '****: ',adjustl(trim(m%namda)),'/',adjustl(trim(auxKey)),' Case/Driver detected and error'
             write(*,*) '****:Folded Elements, attempting to shutdown THIS case smoothly'
         endif
         !call runend('Folded Elements')
         !End time step
         call flush
         !call MPI_Barrier(m%MPIcomm, ierr)
         call a%Endste

         !Case output
         call flush
         call a%Output
         !call MPI_Barrier(m%MPIcomm, ierr)

         call a%Turnof

     endif
   
   endif
   
   call d%cpu_MeshProjections%Toc
   
end subroutine
