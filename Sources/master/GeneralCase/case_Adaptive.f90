module Mod_CaseAdaptive
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_DriverInterface
   implicit none


   type(caseVariables), pointer :: c => NULL()
   type(masterVariables), pointer :: m => NULL()
   type(domainVariables), pointer :: d => NULL()
   type(adaptiveVariables), pointer :: ad => NULL()

contains

   subroutine LoopPreRefine(myDC)
      implicit none

      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()

      call ExtractDriver(myDC,myDriver)
      call myDriver%PreRefine(c)
   end subroutine


   subroutine LoopRefine(myDC)
      implicit none

      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()

      call ExtractDriver(myDC,myDriver)
      call myDriver%Refine(c)
   end subroutine

   subroutine LoopRebalance(myDC)
      implicit none

      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()

      call ExtractDriver(myDC,myDriver)
      call myDriver%Rebalance(c)
   end subroutine
end module

subroutine case_Adaptive(a,itask)
   !Deals with the adaptive process at the beginning of each time step
   use typre
   use Mod_GeneralCase
   use Mod_caseVariables
   use def_parame
   use Mod_int2str
   use Mod_Mesh_AR_Interface
   use Mod_Refiners
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_DriverInterface
   use Mod_CaseAdaptive
   use MPI
   implicit none
   class(GeneralCase), target :: a
   integer(ip) :: itask

   integer(ip) :: nelem
   integer(ip), allocatable :: level(:)
   character(150) :: fil_outpu_dom, fil_postp,fil_postp_graf,extra
   logical :: kfl_DoRebalancing
   type(Mesh_IniData) :: MyIniData
   integer(ip), allocatable :: processor(:)
   integer(ip) :: npoin,npoinLocal
   integer(ip), allocatable :: PointGlobNumber(:),PointProcNumber(:)
   real(rp), pointer :: coord(:,:) => NULL()
   integer(ip) :: iaux,maxmark
   integer(ip) :: ierr

   character(DCHashCharSize) :: DoFirstKeys(1)


   character(150) :: namda

   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer      :: myDriver => NULL()

   interface
      subroutine RebalancingCriteria(a,kfl_DoRebalancing)
         use Mod_GeneralCase
         class(GeneralCase), target :: a
         logical :: kfl_DoRebalancing
      end subroutine

      subroutine case_EnsureMaxLevelLimit(a)
         use typre
         use Mod_GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine

   end interface

   c   => a%caseVars
   m   => a%caseVars%masterVars
   d   => a%caseVars%domainVars
   ad  => a%caseVars%adaptiveVars


   call ad%cpu_adaptive%Tic

   if (itask == 1) then
      !Refiner initializations

      if (ad%kfl_AdaptiveRefinement == 1) then
         call MyIniData%Initialize(d%Mesh)

         !Only one type of adaptive refiner available
         ad%Refiner => IsotropicAdaptiveRefiner_Const()

         call ad%Refiner%SetMPI(m%MPIcomm,m%MPIsize,m%MPIroot,m%MPIrank)
         call ad%Refiner%SetOutputFiles(m%lun_memor,m%lun_outpu)
         call ad%Refiner%Initialize(MyIniData)
         call ad%Refiner%Set_2_1_Balancing(ad%kfl_21Balancing)

         call d%Mesh%InitializeAdaptiveRefinement(ad%Refiner)

      endif
   elseif (itask == 2) then

      if (ad%kfl_AdaptiveRefinement == 1 .and. c%masterVars%istep > ad%AdaptiveDelay) then
         !Refine
         call d%Mesh%GetNelem(nelem)
         call m%MasterMemo%alloc(nelem,ad%RefinerMarkel,'RefinerMarkel','Begste')

         if (ad%RefinementLeader == 'ALLRE') then
            ad%RefinerMarkel = 1
         else
            !Get Refinement Criteria from the Refinement leader
            call a%DriverList%GetFromKey(ad%RefinementLeader,myDC)
            if (.not. associated(myDC)) call runend('Adaptive: RefinementLeader not found')
            call ExtractDriver(myDC,myDriver)
            !This sets the RefinerMarkel array, we should change the interface (not the full c)
            call myDriver%GetRefinementCriteria(c)
            !Limit the maximum refinement level if required
            call case_EnsureMaxLevelLimit(a)
         endif

         !Do not refine if there is nothing to do
         iaux = maxval(abs(ad%RefinerMarkel))
         call MPI_ALLREDUCE(iaux,maxmark,1,MPI_INTEGER4,MPI_MAX,m%MPIcomm,ierr)
         if (maxmark == 0) then
            call m%MasterMemo%dealloc(nelem,ad%RefinerMarkel,'RefinerMarkel','Begste')
            return
         endif

         !Continue with the preparation of the adaptive refinement process
         call ad%Refiner%PrepareOldExternalBoundariesInfoToKeep
         call d%Mesh%MatchOldBoundariesToRefiner

         !Prerefine
         call a%DriverList%LoopList(LoopPreRefine)

         !Coord for setting tetra variations
         call d%Mesh%GetCoord(coord)
         call ad%Refiner%SetCoordArray(coord)

         !Refine the refiner and the mesh
         call ad%Refiner%Refine(ad%RefinerMarkel)
         call d%Mesh%Refine('Refine')

         call m%MasterMemo%dealloc(nelem,ad%RefinerMarkel,'RefinerMarkel','Begste')

         if (ad%UpdateIO) then
            if (mod(m%istep,ad%UpdateFreq) == 0) then
               call a%SetupIO(1)
               call a%SetupIO(0)
            end if
         end if

         !For each driver, do their adaptive refinement process
         !ALE needs to be done first
         DoFirstKeys(1) = 'ALEPR'
         call a%DriverList%LoopList(LoopRefine,DoFirstKeys(1:1))

         !Load Rebalancing Criteria
         call RebalancingCriteria(a,kfl_DoRebalancing)

         if (kfl_DoRebalancing .eqv. .true.) then

            call ad%Refiner%PrepareOldExternalBoundariesInfoToKeep
            call d%Mesh%MatchOldBoundariesToRefiner
            call d%Mesh%GetNpoin(npoin)

            call m%MasterMemo%alloc(npoin,PointGlobNumber,'PointGlobNumber','Adaptive')
            call m%MasterMemo%alloc(npoin,PointProcNumber,'PointProcNumber','Adaptive')

            call d%Mesh%ComputeRebalancingNumbering(PointGlobNumber,PointProcNumber)
            call ad%Refiner%SetRebalancingNumbering(PointGlobNumber,PointProcNumber)

            call m%MasterMemo%dealloc(npoin,PointGlobNumber,'PointGlobNumber','Adaptive')
            call m%MasterMemo%dealloc(npoin,PointProcNumber,'PointProcNumber','Adaptive')

            call ad%Refiner%LoadRebalance
            call d%Mesh%Refine('Rebalance')

            !For each driver, do the load rebalancing process
            !ALE needs to be done first
            DoFirstKeys(1) = 'ALEPR'
            call a%DriverList%LoopList(LoopRebalance,DoFirstKeys(1:1))
         endif

         !For each driver, reset the pointers for the outChannel, so that anybody can find them
         call a%DriverList%GetFirst(myDC)
         do while (associated(myDC))
            call ExtractDriver(myDC,myDriver)
            call myDriver%SetOutChannels('UPDATE')
            call a%DriverList%GetNext(myDC)
         enddo

         !PostProcessing
         !close postprocess file and open a new file
         call m%FilePostpr%ClosePostFile
         call m%FilePostpr%CloseGrafFile
         m%meshCounter = m%meshCounter+1
         fil_postp_graf = trim(m%PostProcessFolder)//'/'//adjustl(trim(m%namda))//adjustl(trim(int2str(m%MPIrank)))
         fil_postp_graf = trim(fil_postp_graf)//'mesh'//adjustl(trim(int2str(m%meshCounter)))//'.post.grf'
         call m%FilePostpr%OpenPostFile(m%PostProcessFolder,m%namda,m%meshCounter)
         call m%FilePostpr%OpenGrafFile(fil_postp_graf)
         call m%FilePostpr%Postpr(d%Mesh)

      endif

   elseif (itask == 100) then
      !Deallocate Refiner
      if (ad%kfl_AdaptiveRefinement == 1) then
         call ad%Refiner%Dealloc
      endif
   elseif (itask == 101) then
      if (ad%kfl_AdaptiveRefinement == 1) then
         !Deallocate it
         deallocate(ad%Refiner)
      endif
   endif

   call ad%cpu_adaptive%Toc

contains

end subroutine


subroutine RebalancingCriteria(a,kfl_DoRebalancing)
   use typre
   use Mod_CaseVariables
   use Mod_GeneralCase
   use MPI
   implicit none
   class(GeneralCase), target :: a
   logical :: kfl_DoRebalancing

   integer(ip) :: npoinLocal,maxnpoinLocal,minnpoinLocal,gnpoin
   integer :: ierr
   real(rp) :: ratio
   real(rp) :: optimumnpoin

   type(domainVariables), pointer :: d => NULL()
   type(masterVariables), pointer :: m => NULL()
   type(adaptiveVariables), pointer :: ad => NULL()

   d => a%caseVars%domainVars
   m => a%caseVars%masterVars
   ad => a%caseVars%adaptiveVars

   kfl_DoRebalancing = .false.

   call d%Mesh%GetNpoinLocal(npoinLocal)
   call d%Mesh%GetGnpoin(gnpoin)

   optimumnpoin = real(gnpoin)/m%MPIsize
   call MPI_REDUCE( npoinLocal, maxnpoinLocal, 1, MPI_INTEGER4, MPI_MAX, m%MPIroot,m%MPIcomm, ierr )

   if (m%MPIrank == m%MPIroot) then
      ratio = real(maxnpoinLocal)/optimumnpoin
      if (ratio > ad%RefinerRebalancingRatio) kfl_DoRebalancing = .true.
   endif

   CALL MPI_BCAST(kfl_DoRebalancing, 1, MPI_LOGICAL, m%MPIroot, m%MPIcomm, ierr)
end subroutine

subroutine case_EnsureMaxLevelLimit(a)
   use typre
   use Mod_CaseVariables
   use Mod_AdaptiveVariables
   use Mod_GeneralCase
   use Mod_Element
   implicit none
   class(GeneralCase), target :: a


   integer(ip) :: npoin
   integer(ip), allocatable :: level(:)
   class(FiniteElement), pointer :: e => NULL()

   type(caseVariables), pointer :: c => NULL()
   type(domainVariables), pointer :: d => NULL()
   type(masterVariables), pointer :: m => NULL()
   type(adaptiveVariables), pointer :: ad => NULL()

   integer(ip) :: ielem,nelem,maxlevel
   integer(ip) :: CurrentAdaptiveMaxLevel

   c => a%caseVars
   d => a%caseVars%domainVars
   m => a%caseVars%masterVars
   ad => a%caseVars%adaptiveVars
   
   currentAdaptiveMaxLevel = ad%RefinerMaxLevel
   if (currentAdaptiveMaxlevel == 0) currentAdaptiveMaxLevel = -1

   if (ad%ProgressiveMaxLevel > 0) then
      currentAdaptiveMaxLevel = min(currentAdaptiveMaxLevel,m%istep/ad%ProgressiveMaxLevel)
   endif


   !if (currentAdaptiveMaxLevel > 0) then
      call d%Mesh%GetNpoin(npoin)
      call m%MasterMemo%alloc(npoin,level,'level','ApplyErrorCriteria')
      call ad%Refiner%GetLevel(level)

      call d%Mesh%ElementAlloc(e,m%MasterMemo,'DefaultRule','nsm_EnditeElmope')


      call d%Mesh%GetNelem(nelem)
      do ielem = 1,nelem
         !Load Element
         call d%Mesh%ElementLoad(ielem,e)
         maxlevel = maxval(level(e%lnods(1:e%pnode)))

         if ((currentAdaptiveMaxLevel > -1) .and. (maxlevel >= currentAdaptiveMaxLevel) .and. (ad%RefinerMarkel(ielem) == 1)) ad%RefinerMarkel(ielem) = 0
         if (ad%RefinerMarkel(ielem) == -1 .and. (maxlevel == 0)) ad%RefinerMarkel(ielem) = 0

      enddo


      call m%MasterMemo%dealloc(npoin,level,'level','ApplyErrorCriteria')
      call d%Mesh%ElementDeAlloc(e,m%MasterMemo,'DefaultRule','nsm_EnditeElmope')


   !endif


end subroutine
