module Mod_ElementWithDataStructures
   use typre
   use Mod_Element
   use Mod_ElementDataStructures
   implicit none
   private
   public FiniteElementWDS
   
   type, extends(FiniteElement) ::  FiniteElementWDS
  
      type(ElementDataStructure), pointer :: EDS => NULL()
      type(BoundaryDataStructure), pointer :: BDS => NULL()

contains

      procedure :: alloc  => ElementAlloc
      procedure :: dealloc=> ElementDealloc
      procedure :: elmdel => ElementElmdel
      procedure :: elmder => ElementElmder
      procedure :: elmdcg => ElementElmdcg
      procedure :: elmhes => ElementElmhes
      procedure :: elmderb=> ElementElmderb
      procedure :: elmlen => ElementElmlen
      procedure :: bounor => ElementBounor
      procedure :: elenor => ElementElenor
      
   end type

contains   
   
   subroutine ElementAlloc(e,Memor,runam)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElementWDS) :: e
      character(*) :: runam
      type(MemoryMan) :: Memor

      call Memor%palloc(e%mnode,e%lnods,'e%lnods',runam)
      call Memor%alloc(e%ndime,e%mnode,e%elcod,'e%elcod',runam)
      
      call Memor%alloc(e%mnodb,e%lnodb,'e%lnodb',runam)
      call Memor%alloc(e%ndime,e%mnodb,e%bocod,'e%bocod',runam)
      
      
   end subroutine
   
   
   
   subroutine ElementDealloc(e,Memor,runam)
      use typre
      use Mod_Memor
      implicit none
      class(FiniteElementWDS) :: e
      type(MemoryMan) :: Memor
      character(*) :: runam
      
      call Memor%pdealloc(size(e%lnods,1),e%lnods,'e%lnods',runam)
      call Memor%dealloc(e%ndime,size(e%elcod,2),e%elcod,'e%elcod',runam)
      
      call Memor%dealloc(e%mnodb,e%lnodb,'e%lnodb',runam)
      call Memor%dealloc(e%ndime,e%mnodb,e%bocod,'e%bocod',runam)
      
   end subroutine
   
   subroutine ElementElmdel(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      if (e%linea == 1) then
         e%detjm = e%EDS%detjmg
         e%cartd => e%EDS%cartg
      endif

   end subroutine
   
   subroutine ElementElmdcg(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      e%detjm = e%EDS%detjmg
      e%cartd => e%EDS%cartg
      
   end subroutine
  
  
   subroutine ElementElmder(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      e%detjm = e%EDS%detjm(e%igaus)
      e%cartd => e%EDS%cartd(:,:,e%igaus)
      
      
   end subroutine
   
   subroutine ElementElmhes(e)
      implicit none
      class(FiniteElementWDS) :: e
   
      if (e%linea == 1) then
      else   
         e%hessi => e%EDS%hessi(:,:,e%igaus)
      endif   
   end subroutine
   
   !Calculate cartesian derivates at boundary gauss points
   subroutine ElementElmderb(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      e%cartb => e%BDS%cartb(:,:,e%igaub)
      
   end subroutine
   
   subroutine ElementBounor(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      e%baloc => e%BDS%baloc(:,:,e%igaub)
      e%eucta = e%BDS%eucta(e%igaub)
      
   end subroutine   
   
   subroutine ElementElenor(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      e%baloc => e%EDS%elvec(:,:,e%iface,e%igaub)
      e%eucta = e%EDS%elnor(e%iface,e%igaub)
      
   end subroutine   
   
   subroutine ElementElmlen(e)
      implicit none
      class(FiniteElementWDS) :: e
      
      e%hleng => e%EDS%hleng
      e%tragl => e%EDS%tragl
   
   end subroutine
   
end module
