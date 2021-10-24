module Mod_plcd_SIMP_ChiFilter
   use typre
   use Mod_PLCD
   use Mod_Element
   use Mod_Postpr
   use Mod_int2str
   implicit none
   private
   public :: ElementCenterChiFilter,ChiFilterInterface,UpdateChiValue


   type, abstract :: ChiFilterInterface
      class(PLCDProblem), pointer :: plcd => NULL()


contains
      procedure(Initialize0), deferred :: Initialize
      procedure(Finalize0), deferred :: Finalize
      procedure(AddGaussPointContributionToLambda0), deferred :: AddGaussPointContributionToStrainEnergy
      procedure(PostGaussLambda0), deferred :: PostGaussStrainEnergy
      procedure(FinalizeLambdaComputations0), deferred :: FinalizeStrainEnergyComputations
      procedure(ResetChi0), deferred :: ResetChi
      procedure(AddGaussPointContributionToChi0), deferred :: AddGaussPointContributionToChi
      procedure(PostGaussChi0), deferred :: PostGaussChi
      procedure(FinalizeChiComputations0), deferred :: FinalizeChiComputations
      procedure(UpdateChi0), deferred :: UpdateChi
      procedure(GetGaussPointChi0), deferred :: GetGaussPointChi
      procedure(GetNodalChi0), deferred :: GetNodalChi
      procedure(PostprocessChi0), deferred :: PostprocessChi


   end type

   type, extends(ChiFilterInterface) :: ElementCenterChiFilter
      real(rp), allocatable :: EWLambda(:), EWChi(:)

contains
      procedure :: Initialize
      procedure :: Finalize
      procedure :: AddGaussPointContributionToStrainEnergy
      procedure :: PostGaussStrainEnergy
      procedure :: FinalizeStrainEnergyComputations
      procedure :: ResetChi
      procedure :: AddGaussPointContributionToChi
      procedure :: PostGaussChi
      procedure :: FinalizeChiComputations
      procedure :: UpdateChi
      procedure :: GetGaussPointChi
      procedure :: PostprocessChi
      procedure :: GetNodalChi
   end type


   abstract interface
      subroutine Initialize0(a,plcd)
         import
         class(ChiFilterInterface) :: a
         class(PLCDProblem), target :: plcd
         integer(ip) :: nelem
      end subroutine

      subroutine Finalize0(a)
         import
         class(ChiFilterInterface) :: a
      end subroutine

      subroutine AddGaussPointContributionToLambda0(a,ielem,e,StrainEnergy)
         import
         class(ChiFilterInterface) :: a
         integer(ip) :: ielem
         class(FiniteElement) :: e
         real(rp) :: StrainEnergy
      end subroutine

      subroutine PostGaussLambda0(a,ielem,e)
         import
         class(ChiFilterInterface) :: a
         integer(ip) :: ielem
         class(FiniteElement) :: e
      end subroutine

      subroutine FinalizeLambdaComputations0(a)
         import
         class(ChiFilterInterface) :: a
      end subroutine

      subroutine ResetChi0(a)
         import
         class(ChiFilterInterface) :: a
      end subroutine

      subroutine AddGaussPointContributionToChi0(a,ielem,e,chi)
         import
         class(ChiFilterInterface) :: a
         integer(ip) :: ielem
         class(FiniteElement) :: e
         real(rp) :: chi
      end subroutine

      subroutine PostGaussChi0(a,ielem,e)
         import
         class(ChiFilterInterface) :: a
         integer(ip) :: ielem
         class(FiniteElement) :: e
      end subroutine

      subroutine FinalizeChiComputations0(a)
         import
         class(ChiFilterInterface) :: a
      end subroutine

      subroutine UpdateChi0(a,lambda)
         import
         class(ChiFilterInterface) :: a
         real(rp) :: lambda
      end subroutine

      subroutine GetGaussPointChi0(a,ielem,e,Chi)
         import
         class(ChiFilterInterface) :: a
         integer(ip) :: ielem
         class(FiniteElement) :: e
         real(rp) :: Chi
      end subroutine

      subroutine PostprocessChi0(a,FilePostpr)
         import
         class(ChiFilterInterface) :: a
         class(PostprFile) :: FilePostpr
      end subroutine

      subroutine GetNodalChi0(a,NodalChi)
         import
         class(ChiFilterInterface), target :: a
         real(rp), pointer :: NodalChi(:)
      end subroutine

   end interface




contains

   subroutine Initialize(a,plcd)
      class(ElementCenterChiFilter) :: a
      class(PLCDProblem), target :: plcd
      integer(ip) :: nelem

      a%plcd => plcd

      call a%plcd%Mesh%GetNelem(nelem)
      call a%plcd%Memor%alloc(nelem,a%EWLambda,'EWLambda','EW')
      call a%plcd%Memor%alloc(nelem,a%EWChi,'EWChi','EW')



   end subroutine

   subroutine Finalize(a)
      class(ElementCenterChiFilter) :: a

      integer(ip) :: nelem


      call a%plcd%Mesh%GetNelem(nelem)
      call a%plcd%Memor%dealloc(nelem,a%EWLambda,'EWLambda','EW')
      call a%plcd%Memor%dealloc(nelem,a%EWChi,'EWChi','EW')

   end subroutine

   subroutine AddGaussPointContributionToStrainEnergy(a,ielem,e,StrainEnergy)
      class(ElementCenterChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: StrainEnergy

      a%EWLambda(ielem) = a%EWLambda(ielem) + StrainEnergy*e%weigp(e%igaus)/sum(e%weigp(1:e%pgaus))
   end subroutine

   subroutine PostGaussStrainEnergy(a,ielem,e)
      class(ElementCenterChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
   end subroutine

   subroutine FinalizeStrainEnergyComputations(a)
      class(ElementCenterChiFilter) :: a

   end subroutine

   subroutine ResetChi(a)
      class(ElementCenterChiFilter) :: a

      a%EWChi = 0.0_rp
   end subroutine

   subroutine AddGaussPointContributionToChi(a,ielem,e,chi)
      class(ElementCenterChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: chi

      a%EWChi(ielem) = chi
   end subroutine

   subroutine PostGaussChi(a,ielem,e)
      class(ElementCenterChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
   end subroutine

   subroutine FinalizeChiComputations(a)
      class(ElementCenterChiFilter) :: a

   end subroutine

   subroutine UpdateChi(a,lambda)
      class(ElementCenterChiFilter) :: a
      real(rp) :: lambda

      integer(ip) :: nelem,ielem
      real(rp) :: chi,BK
      real(rp) :: newchi,compareValue

      call a%plcd%Mesh%GetNelem(nelem)
      do ielem = 1,nelem
         BK = a%EWLambda(ielem)/lambda
         chi = a%EWChi(ielem)

         call UpdateChiValue(a,BK,chi,newchi)

         a%EWChi(ielem) = newchi
      enddo
   end subroutine

   subroutine UpdateChiValue(a,BK,chi,newchi)
      implicit none
      class(ChiFilterInterface) :: a
      real(rp) :: chi,BK
      real(rp) :: newchi

      real(rp) :: compareValue

      compareValue = chi*(BK**a%plcd%SIMPData%eta)
      if (compareValue <= max((1-a%plcd%SIMPData%xi)*chi,a%plcd%SIMPData%chimin)) then
         newchi = max((1-a%plcd%SIMPData%xi)*chi,a%plcd%SIMPData%chimin)

      !elseif (compareValue >= min((1+a%plcd%SIMPData%xi)*chi,a%plcd%SIMPData%chimax)) then
      !   newchi = min((1+a%plcd%SIMPData%xi)*chi,a%plcd%SIMPData%chimax)

      elseif (compareValue >= (1+a%plcd%SIMPData%xi)*chi) then
         newchi = (1+a%plcd%SIMPData%xi)*chi
      else
         newchi = compareValue;
      endif
   end subroutine

   subroutine GetGaussPointChi(a,ielem,e,Chi)
      class(ElementCenterChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: Chi

      Chi = a%EWChi(ielem)
   end subroutine

   subroutine PostprocessChi(a,FilePostpr)
      class(ElementCenterChiFilter) :: a


      class(PostprFile) :: FilePostpr

      call a%plcd%FilePostpr%postgp(a%EWChi,'ChiElementWise'//adjustl(trim(int2str(a%plcd%itera))),a%plcd%istep,a%plcd%ctime,a%plcd%Mesh)
   end subroutine

   subroutine GetNodalChi(a,NodalChi)
      class(ElementCenterChiFilter), target :: a
      real(rp), pointer :: NodalChi(:)

      NodalChi => NULL()
   end subroutine








end module
