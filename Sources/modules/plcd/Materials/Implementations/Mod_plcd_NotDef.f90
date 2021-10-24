module Mod_plcd_Notdef
   use typre
   use Mod_plcd_Material
   use Mod_Listen
   implicit none
   private
   public Notdef
   
   type, extends(PLCDMaterial) :: Notdef
   
contains
      procedure :: SpecificCreateElementData => SpecificCreateElementData
      procedure :: ComputeStressTensor => CSLinElas
      procedure :: ComputeTangentTensor => CTLinElas
      procedure :: SpecificReadData => LE_SpecificReadData
      procedure :: SpecificScatterData => LE_SpecificScatterData
   end type

contains
   
   subroutine CSLinElas(a,S)
      class(Notdef):: a
      real(rp) :: S(:)
      call runend('Notdef material cannot compute stresses')
   end subroutine
   
   subroutine CTLinElas(a,T)
      class(Notdef):: a
      real(rp) :: T(:,:)
      call runend('Notdef material cannot compute tangent matrix')
   end subroutine
   
   subroutine LE_SpecificReadData(a,Listener)
      class(Notdef) :: a
      type(ListenFile) :: Listener
   end subroutine
   
   subroutine LE_SpecificScatterData(a)
      class(Notdef) :: a
   end subroutine
   
   subroutine SpecificCreateElementData(a,pgaus,ElementData)
      implicit none
      class(NotDef), target :: a
      integer(ip) :: pgaus
      class(ElementMaterialData), pointer :: ElementData
      
      call runend('Notdef material cannot create gausspointdata')
   end subroutine
   
end module