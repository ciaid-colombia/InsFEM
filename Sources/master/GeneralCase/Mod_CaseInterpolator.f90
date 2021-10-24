module Mod_CaseInterpolator
   use typre
   use Mod_MeshInterpolator
   use Mod_Listen
   use Mod_BroadCastBuffer
   use Mod_GeneralCase
   use Mod_Mesh
   use Mod_Memor
   use Mod_Element
   implicit none
   private
   public CaseInterpolator, CaseInterpolatorDataFromListener, AddCaseInterpolatorToBroadCastBuffer,CaseInterpolator_Const
   
   type :: CaseInterpolator
      character(5) :: InterpolationKind = 'VOLUM'    !Can also be 'BOUND'
      character(5) :: InterpInitSearch = 'DEFOR'     !Can also be 'ORIGI'
      integer(ip)  :: ibody = 0                      !If boundary, default is the whole boundary
      
      logical                    :: IsInitializedInterpolator = .false.
      type(GeneralCase), pointer  :: NeighborCase => NULL()
      class(GeneralCase), pointer :: MyCase => NULL()
      
      integer(ip) :: npoin
      real(rp)    :: InterpolatorTolerance
      type(Interpolator)         :: Interpolator
      
   
contains
      procedure :: SetNeighborCase
      procedure :: SetMyCase
      procedure :: SetInterpolatorTolerance
      procedure :: SetupInterpolator
      procedure :: ReSetupInterpolator
      procedure :: FinalizeInterpolator
      procedure :: GetInterpolatorAndNpoin
   
   end type
   
   interface CaseInterpolator_Const
      procedure constructor
   end interface
   
contains

   function constructor()
      type(CaseInterpolator), pointer :: constructor
      
      allocate(constructor)
   end function
   
   subroutine SetMyCase(a,MyCase)
      class(CaseInterpolator) :: a
      class(GeneralCase), target :: MyCase
      
      a%MyCase => MyCase
   end subroutine
   
   subroutine SetInterpolatorTolerance(a,Tolerance)
      class(CaseInterpolator) :: a
      real(rp) :: Tolerance   
      
      call a%Interpolator%SetTol(Tolerance)
   end subroutine
   
   subroutine SetNeighborCase(a,NeighborCase)
      class(CaseInterpolator) :: a
      type(GeneralCase), target :: NeighborCase
      
      a%NeighborCase => NeighborCase
   end subroutine
   
   subroutine CaseInterpolatorDataFromListener(a,Listener)
      type(CaseInterpolator) :: a
      type(ListenFile) :: Listener
      
      if (Listener%nnwor > 2) then
         a%InterpolationKind = Listener%words(3)
         if (a%InterpolationKind /= 'VOLUM' .and. a%InterpolationKind /= 'BOUND') then 
            call runend('Wrong type of Caseinterpolator type')
         endif

         a%InterpInitSearch= Listener%words(4)
         if (a%InterpInitSearch /= 'DEFOR' .and. a%InterpInitSearch /= 'ORIGI') then 
             a%InterpInitSearch= 'DEFOR' !Use default
         endif
         a%ibody = Listener%param(3)
      endif
   end subroutine
   
   subroutine AddCaseInterpolatorToBroadCastBuffer(a,bb)
      type(CaseInterpolator) :: a
      type(BroadCastBuffer) :: bb
      
      call bb%Add(5,a%InterpolationKind)
      call bb%Add(5,a%InterpInitSearch)
      call bb%Add(a%ibody)
   end subroutine
   
   subroutine ReSetupInterpolator(a,myMesh,Memor)
      class(CaseInterpolator) :: a
      type(FemMesh)     :: MyMesh
      type(MemoryMan)   :: Memor
      
      type(FemMesh), pointer :: NeighborMesh => NULL()
      logical :: isALE = .false. , isAdaptive = .false.,needstobereset = .false.
      
      !Now I ask my Neighbor to Initialize an interpolator for me
      NeighborMesh => a%NeighborCase%caseVars%domainVars%Mesh
      call NeighborMesh%GetALE(isALE)
      
      !Also I check if the neighbor is adaptive
      !or if I am adaptive
      if (a%NeighborCase%caseVars%adaptiveVars%kfl_AdaptiveRefinement == 1) then
         isAdaptive = .true.
      endif
      if (a%MyCase%caseVars%adaptiveVars%kfl_AdaptiveRefinement == 1) then
         isAdaptive = .true.
      endif
      
      needstobereset = .false.
      if (isALE .and. a%InterpolationKind == 'VOLUM') then
         needstobereset = .true.
      endif
      if (isAdaptive) needstobereset = .true.
      
      
      !If it is a volume interpolation, and the neighbor is ALE: resetup interpolators
      if (needstobereset) then
         !if (a%InterpolationKind == 'VOLUM') then
            call a%setupInterpolator(myMesh,Memor)
         !endif
      endif
   end subroutine
      
      
   
   subroutine SetupInterpolator(a,myMesh,Memor)
      class(CaseInterpolator) :: a
      type(FemMesh)     :: MyMesh
      type(MemoryMan)   :: Memor
      
      type(FemMesh), pointer :: NeighborMesh => NULL()
      real(rp), pointer :: coord(:,:) => NULL()
      integer(ip) :: npoin,nboun,ipoin,inodb,iboun
      
      integer(ip), pointer :: ibody => NULL()
      
      class(FiniteElement), pointer :: e => NULL()
      logical, allocatable :: DoFlag(:)
      integer(ip) :: lumem
      logical::isNeighborALE=.false.
      logical::turnNeighborALEOn=.false.
      
      if (a%IsInitializedInterpolator) then
         call a%Interpolator%Finalize
         call Memor%deallocObj(0,'Interpolator','SetupInterpolator',1)
      endif
      a%IsInitializedInterpolator = .true.
      call Memor%allocObj(0,'Interpolator','SetupInterpolator',1)
      
      call myMesh%GetNpoin(npoin)
      call myMesh%GetCoord(coord)

      a%npoin=npoin
      
      !Now I ask my Neighbor to Initialize an interpolator for me
      NeighborMesh => a%NeighborCase%caseVars%domainVars%Mesh

      !ask if neighbor mesh is ALE, if it is we ask for undeformed coords
      call NeighborMesh%GetALE(isNeighborALE)
      if(isNeighborALE .and. a%InterpInitSearch=='ORIGI') then
          call NeighborMesh%SetALE(0)
          turnNeighborALEOn=.true.
      endif


      call Memor%GetLumem(lumem)
      call a%Interpolator%SetLOutputFile(lumem,0)
      
      if (a%InterpolationKind == 'VOLUM') then
         call a%Interpolator%Initialize_Interp(NeighborMesh,coord)
      
      elseif (a%InterpolationKind == 'BOUND') then
         call Memor%alloc(npoin,doflag,'doflag','SetupInterpolator')
         call myMesh%ElementAlloc(e,Memor,'DefaultRule','SetupInterpolator')      
         
         ! Loop over boundaries, mark the points corresponding to the selected body
         call myMesh%GetNboun(nboun)
         boundaries: do iboun=1,nboun
            !Load Element
            call myMesh%BoundaryLoad(iboun,e)
            
            call myMesh%GetLbody(ibody,iboun)   
            if (ibody == a%ibody .or. a%ibody == 0) then
               !We mark the nodes as to be interpolated
               do inodb = 1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  
                  DoFlag(ipoin) = .true.
                enddo
             endif
         enddo boundaries
         
         call a%Interpolator%Initialize_Interp(NeighborMesh,coord,DoFlag)

         if(turnNeighborALEOn) then
             call NeighborMesh%SetALE(1)
             turnNeighborALEOn=.false.
         endif
         
         call myMesh%ElementDeAlloc(e,Memor,'DefaultRule','SetupInterpolator')      
         call Memor%dealloc(npoin,doflag,'doflag','SetupInterpolator')
      endif
      
   end subroutine
   
   subroutine FinalizeInterpolator(a,Memor)
      class(CaseInterpolator) :: a
      type(MemoryMan) :: Memor
      
      if (a%IsInitializedInterpolator) then
         call a%Interpolator%Finalize
      endif
      a%IsInitializedInterpolator = .false.
      call Memor%deallocObj(0,'Interpolator','SetupInterpolator',1)
      
      
   end subroutine
   
   subroutine GetInterpolatorAndNpoin(a,InterpolatorPointer,npoinPointer)
      class(CaseInterpolator), target :: a
      type(Interpolator), pointer :: InterpolatorPointer
      integer(ip), pointer :: npoinPointer
      
      InterpolatorPointer => a%Interpolator
      npoinPointer => a%npoin
   end subroutine
      
   
end module
