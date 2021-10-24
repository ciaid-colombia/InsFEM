subroutine rom_Refine(a,itask)
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   character(6) :: itask
   character(150) :: optfile
   character(3)   :: exmod
   integer(ip) :: oldnpoin,newnpoin,idofr,idofn,npoinLocal,gnpoin
   real(rp), allocatable :: auxBasis(:,:,:),auxMean(:,:),auxFOM(:,:)

   oldnpoin = size(a%Basis,2)
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetGnpoin(gnpoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   exmod =adjustl(trim(a%exmod))
   optfile = adjustl(trim(a%InputFolder))//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exmod))//'.rom.sol'
   
   call a%Memor%alloc(a%ndofn,newnpoin,a%ndofr,auxBasis,'Basis','rom_Refine')
   do idofr = 1,a%ndofr
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(a%ndofn,a%Basis(:,:,idofr),auxBasis(:,:,idofr))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(a%ndofn,a%Basis(:,:,idofr),auxBasis(:,:,idofr))
      endif   
   enddo
   call move_alloc(auxBasis,a%Basis)
   call a%Memor%deallocObj(0,'Basis','rom_Refine',rp*a%ndofn*a%ndofr*oldnpoin)

   call a%Memor%alloc(a%ndofn,newnpoin,auxMean,'SnapMean','rom_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(a%ndofn,a%SnapMean,auxMean)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(a%ndofn,a%SnapMean,auxMean)
   endif   
   call move_alloc(auxMean,a%SnapMean)
   call a%Memor%deallocObj(0,'SnapMean','rom_Refine',rp*a%ndofn*oldnpoin)

   call a%Memor%alloc(a%ndofn,newnpoin,auxFOM,'FOMSolution','rom_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(a%ndofn,a%FOMSolution,auxFOM)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(a%ndofn,a%FOMSolution,auxFOM)
   endif   
   call move_alloc(auxFOM,a%FOMSolution)
   call a%Memor%deallocObj(0,'FOMSolution','rom_Refine',rp*a%ndofn*oldnpoin)

   call a%EigenSystem%DeallocateRun

   if (a%kfl_updateBasis .eqv. .true. .and. a%kfl_inter == 1) then
      if (mod(a%Problem%istep,a%UpdateFreq) == 0) call a%Interpolate
   end if
   
   call a%EigenSystem%Init(gnpoin,newnpoin,npoinLocal,a%ndofn,optfile,a%exmod,a%Problem%lun_solve,a%Memor)
   call a%EigenSystem%SetLinearSystem(a%Problem%LinearSystem)
   call a%EigenSystem%SetOrderingPointers
   call a%EigenSystem%SetMatrixPointers
   call a%EigenSystem%SetRHSPointers
   call a%EigenSystem%InitSystem
   call a%EigenSystem%InitOperators
   
   call a%EigenSystem%SetNConvergence(a%nconv)
   call a%EigenSystem%SetBasisSize(a%ndofr)
   call a%EigenSystem%SetMean(a%SnapMean(:,1:npoinLocal))
   call a%EigenSystem%SetSigma(a%sigma)
   do idofr=1,a%ndofr
      call a%EigenSystem%SetBasis(idofr,a%Basis(:,1:npoinLocal,idofr))
   end do

   if (a%Mesh%kfl_ProjectorReady == 0) call a%ComputeMassMatrix
end subroutine 

subroutine rom_Interpolate(a)
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   integer(ip)  :: npoin,oldnpoin,idofr
   character(3) :: exmod

   exmod =adjustl(trim(a%exmod))
   if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Basis update, Begin'

   call a%Mesh%GetNpoin(npoin)
   call a%Problem%OldMesh%GetNpoin(oldnpoin)

   do idofr = 1, a%ndofr
      call a%OldMesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%OldMeshBasis(:,:,idofr))
      call a%Int_Restart%Interpolate(a%ndofn,a%OldMeshBasis(:,:,idofr),a%Basis(:,:,idofr))
   end do
   call a%OldMesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%OldMeshSnapMean)
   call a%Int_Restart%Interpolate(a%ndofn,a%OldMeshSnapMean,a%SnapMean) 

   if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Basis update, End'

end subroutine
