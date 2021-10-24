module Mod_ElementalSubroutine
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   use Mod_SupExacso    
   use Mod_SupOperations
   implicit none
    
   type :: PHook
      procedure(), NOPASS, pointer :: PreAllocate => NULL()
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: PreGaussOpt => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: PostGauss => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: Assembly => NULL()
      procedure(), NOPASS, pointer :: AllocateBaseElmopeArrays => NULL()
      procedure(), NOPASS, pointer :: AllocateBaseElmopeMatrices => NULL()
      procedure(), NOPASS, pointer :: DeallocateBaseElmopeArrays => NULL()
      procedure(), NOPASS, pointer :: DeallocateBaseElmopeMatrices => NULL()
      procedure(), NOPASS, pointer :: SetPointers => NULL()
      procedure(), NOPASS, pointer :: PreLoop => NULL()
      procedure(), NOPASS, pointer :: PostLoop => NULL()
   end type
   type(PHook) :: ProcHookE
   

   class(FiniteElement), pointer :: e => NULL()
   type(FemMesh) :: ES_Mesh
   type(MemoryMan) :: ES_Memor
   integer(ip) :: ielem, nelem, ES_ndofn
   integer(ip) :: igaus
   real(rp)    :: dvol
      
   real(rp), allocatable   :: elmat(:,:,:,:), elrhs(:,:)
      
   logical :: flag_ComputeCGDerivatives, flag_ComputeElmlen, flag_ComputeDerivatives,&
              flag_GaussLoop, flag_Assembly, flag_Elmats, flag_element
      
#include "COMPOSEPROCS_POINTERS_50.i90"      
contains
#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   subroutine SetGeneralPointersAndHooksToNULLSUB
      use typre 
      procedure() :: NULLSUB   
      ProcHookE%PreAllocate => NULLSUB
      ProcHookE%Initializations => NULLSUB
      ProcHookE%PreGauss => NULLSUB
      ProcHookE%InGauss => NULLSUB
      ProcHookE%PostGauss => NULLSUB
      ProcHookE%Finalizations => NULLSUB
      ProcHookE%Assembly => NULLSUB
      ProcHookE%AllocateBaseElmopeArrays => NULLSUB
      ProcHookE%AllocateBaseElmopeMatrices => NULLSUB
      ProcHookE%DeallocateBaseElmopeArrays => NULLSUB
      ProcHookE%DeallocateBaseElmopeMatrices => NULLSUB
      ProcHookE%PreGaussOpt => NULLSUB
      ProcHookE%SetPointers => NULLSUB
      ProcHookE%PreLoop => NULLSUB
      ProcHookE%PostLoop => NULLSUB
   end subroutine
   
   subroutine ResetOptions
      flag_ComputeCGDerivatives =.false.
      flag_ComputeElmlen=.false.
      flag_ComputeDerivatives=.false.
      flag_Elmats=.false.
      flag_GaussLoop=.false.
      flag_Assembly=.false. 
      flag_element =.true.
   end subroutine   

   subroutine SetMesh(ExternalMesh)
      type(FemMesh) :: ExternalMesh
      ES_Mesh = ExternalMesh
   end subroutine
   
   subroutine SetMemor(ExternalMemor)
      type(MemoryMan) :: ExternalMemor
      ES_Memor = ExternalMemor
   end subroutine
   
   subroutine Setndofn(Externalndofn)
      integer(ip) :: Externalndofn
      ES_ndofn = Externalndofn
   end subroutine  
   
   subroutine SetComputeElmlen(flag)
      logical :: flag
      flag_ComputeElmlen = flag
      if (flag_ComputeElmLen) then
         call SetComputeCGDerivatives(.true.)
      endif
   end subroutine
   
   subroutine SetComputeCGDerivatives(flag)
      logical :: flag
      flag_ComputeCGDerivatives = flag
   end subroutine
   
   subroutine SetComputeDerivatives(flag)
      logical :: flag
      flag_ComputeDerivatives = flag
   end subroutine
   
   subroutine ComputeElmlen
      call e%elmlen
   end subroutine
   
   subroutine ComputeCGDerivatives
      call e%elmdcg
   end subroutine
   
   subroutine ComputeDerivatives
      call e%elmdel
   end subroutine  
   
   subroutine SetComputeGaussLoop(flag)
      logical :: flag
      flag_GaussLoop = flag
   end subroutine
   
   subroutine SetComputeAssembly(flag)
      logical :: flag
      flag_Assembly = flag
   end subroutine
   
   subroutine SetElmatsToZero(flag)
      logical :: flag
      flag_Elmats = flag
   end subroutine
   
   subroutine ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
   end subroutine   
   
   subroutine AllocateElmats
      call ES_Memor%alloc(ES_ndofn,e%mnode,ES_ndofn,e%mnode,elmat,'elmat','sub_elem')
      call ES_Memor%alloc(ES_ndofn,e%mnode,elrhs,'elrhs','sub_elem')
   end subroutine
   
   subroutine DeallocateElmats
      call ES_Memor%dealloc(ES_ndofn,e%mnode,ES_ndofn,e%mnode,elmat,'elmat','sub_elem')
      call ES_Memor%dealloc(ES_ndofn,e%mnode,elrhs,'elrhs','sub_elem')
   end subroutine
   
   subroutine AllocateElement
      call ES_Mesh%ElementAlloc(e,ES_Memor,'DefaultRule','sub_elem')
   end subroutine
   
   subroutine DeallocateElement
      call ES_Mesh%ElementDealloc(e,ES_Memor,'DefaultRule','sub_elem')
   end subroutine
   
   subroutine SetElementAlloc(flag)
      logical :: flag
      flag_element = flag
   end subroutine   
   

   subroutine ResetPointersAndOptions
       call SetGeneralPointersAndHooksToNULLSUB
       call ResetOptions
   end subroutine
   

   subroutine Initialize
      ielem = 1
      
      if (flag_element .eqv. .true.) then
         call AllocateElement
      end if   
      
      !Set the Pointers for execution
      call ProcHookE%SetPointers
      
      call ProcHookE%PreAllocate

      if (flag_ComputeCGDerivatives .eqv. .true.) then
         call ConcatenateProcedures(ProcHookE%PreGaussOpt, ComputeCGDerivatives)
      end if
      
      if (flag_ComputeElmlen .eqv. .true.) then
         call ConcatenateProcedures(ProcHookE%PreGaussOpt,ComputeElmlen)
      end if
      
      if (flag_ComputeDerivatives .eqv. .true.) then
         call ConcatenateProcedures(ProcHookE%PreGaussOpt, ComputeDerivatives)
      end if
      
     if (flag_Elmats .eqv. .true.) then
        call ConcatenateProcedures(ProcHookE%PreGaussOpt, ElmatsToZero)
        call ConcatenateProcedures(ProcHookE%AllocateBaseElmopeMatrices, AllocateElmats)
        call ConcatenateProcedures(ProcHookE%DeallocateBaseElmopeMatrices, DeallocateElmats)
     end if   
       
      call ProcHookE%AllocateBaseElmopeMatrices
      call ProcHookE%AllocateBaseElmopeArrays
      call ProcHookE%PreLoop
      call ProcHookE%Initializations 
   end subroutine
   
   subroutine GaussLoop
      do igaus = 1,e%pgaus
         e%igaus = igaus
         dvol = e%weigp(e%igaus)*e%detjm
         call ProcHookE%InGauss
      enddo
   end subroutine

   subroutine DoLoop
      call ES_Mesh%GetNelem(nelem)

      do ielem = 1,nelem
         !Load Element
         call ES_Mesh%ElementLoad(ielem,e)
            
         call ProcHookE%PreGaussOpt
         call ProcHookE%PreGauss
         
         if (flag_GaussLoop .eqv. .true.) then
            call GaussLoop
         end if
          
         call ProcHookE%PostGauss
         
         if (flag_Assembly .eqv. .true.) then
            call ProcHookE%Assembly
         end if
      enddo   
   end subroutine
   
   subroutine Finalize
      call ProcHookE%Finalizations
      call ProcHookE%DeallocateBaseElmopeMatrices
      call ProcHookE%DeallocateBaseElmopeArrays
      call ProcHookE%PostLoop
      if (flag_element .eqv. .true.) then 
         call DeallocateElement
      end if  
   end subroutine

end module
