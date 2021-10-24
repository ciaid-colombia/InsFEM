module Mod_ElementalBound
   use typre
   use Mod_Memor
   use Mod_Element   
   use Mod_Mesh    
   implicit none
   
   type :: PHook
      procedure(), NOPASS, pointer :: AllocationArraysAndMatrices => NULL()
      procedure(), NOPASS, pointer :: DeallocationArraysAndMatrices => NULL()
      procedure(), NOPASS, pointer :: Assembly => NULL()
      procedure(), NOPASS, pointer :: DoLoopInitialize => NULL()
      procedure(), NOPASS, pointer :: BoundaryTerms => NULL()
      procedure(), NOPASS, pointer :: BoundaryConditions => NULL()
   end type
   type(PHook) :: ProcHookE
   
   
   class(FiniteElement), pointer :: e => NULL()
   type(FemMesh) :: EB_Mesh
   type(MemoryMan) :: EB_Memor
   integer(ip) :: EB_ndofn
   real(rp), allocatable   :: elmat(:,:,:,:), elrhs(:,:)
   real(rp), allocatable :: welmat(:,:,:,:)
   real(rp), allocatable :: welrhs(:,:)   
   integer(ip) :: iboun,nboun,nelem,ndime,bcstart

   logical :: flag_Elmats
   
#include "COMPOSEPROCS_POINTERS_50.i90"      
contains
#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   subroutine ResetPointers
      procedure() :: NULLSUB
      ProcHookE%AllocationArraysAndMatrices => NULLSUB
      ProcHookE%DeallocationArraysAndMatrices => NULLSUB
      ProcHookE%Assembly => NULLSUB
      ProcHookE%DoLoopInitialize => NULLSUB
      ProcHookE%BoundaryTerms => NULLSUB
      ProcHookE%BoundaryConditions => NULLSUB
   end subroutine
   
   subroutine ResetOptions
      flag_Elmats = .false.
   end subroutine
   
   subroutine SetMesh(ExternalMesh)
      type(FemMesh)  :: ExternalMesh
      EB_Mesh = ExternalMesh
   end subroutine
   
   subroutine SetMemor(ExternalMemor)
      type(MemoryMan) :: ExternalMemor
      EB_Memor = ExternalMemor
   end subroutine
   
   subroutine Setndofn(Externalndofn)
      integer(ip) :: Externalndofn
      EB_ndofn = Externalndofn
   end subroutine   
   
   subroutine ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
   end subroutine   
  
   subroutine AllocateElmats
      call EB_Memor%alloc(EB_ndofn,e%mnode,EB_ndofn,e%mnode,elmat,'elmat','sub_bouelem')
      call EB_Memor%alloc(EB_ndofn,e%mnode,elrhs,'elrhs','sub_bouelem')
   end subroutine
   
   subroutine DeallocateElmats
      call EB_Memor%dealloc(EB_ndofn,e%mnode,EB_ndofn,e%mnode,elmat,'elmat','sub_bouelem')
      call EB_Memor%dealloc(EB_ndofn,e%mnode,elrhs,'elrhs','sub_bouelem')
   end subroutine
   
   subroutine SetElmats(flag)
      logical :: flag
      flag_Elmats = flag
   end subroutine   
      
   
   subroutine Initialize
  
      !Initializations
      call EB_Mesh%GetNelem(nelem)
      call EB_Mesh%GetNdime(ndime)
      call EB_Mesh%GetNboun(nboun)
      
      !Memory allocation
      allocate(e)
      call EB_Mesh%ElementSetSizes(e)
      call e%alloc(EB_Memor,'sub_bouelem')
      call EB_Mesh%ElementSetPointers(e)
      call EB_Mesh%BoundarySetPointers(e)
      !call EB_Mesh%ElementAlloc(e,EB_Memor,'DefaultRule','EB')
      
   
      if (flag_Elmats .eqv. .true.) then
         call AllocateElmats 
         ProcHookE%DeallocationArraysAndMatrices => DeallocateElmats
         ProcHookE%DoLoopInitialize => ElmatsToZero
         
         
      end if
      
      call ProcHookE%AllocationArraysAndMatrices
   end subroutine
   
   subroutine DoLoop
      
      !Loop over boundaries
      boundaries: do iboun=1,nboun
         !Load Element
         call EB_Mesh%BoundaryLoad(iboun,e)

         call ProcHookE%DoLoopInitialize
      
         call ProcHookE%BoundaryTerms
         call ProcHookE%BoundaryConditions
         call ProcHookE%Assembly
      
      end do boundaries
   
   end subroutine
   
   subroutine Finalize
      call ProcHookE%DeallocationArraysAndMatrices
      
      call e%dealloc(EB_Memor,'sub_bouelem')
      !call EB_Mesh%ElementDeAlloc(e,EB_Memor,'DefaultRule','EB')
   end subroutine

end module
