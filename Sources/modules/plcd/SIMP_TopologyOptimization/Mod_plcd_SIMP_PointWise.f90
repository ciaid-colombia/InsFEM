module Mod_plcd_SIMP_PointWise
   use typre
   use Mod_PLCD
   use Mod_Element
   use Mod_Postpr
   use Mod_int2str
   use Mod_plcd_SIMP_ChiFilter
   implicit none
   private
   public :: LumpedMassMatrixChiFilter


   type, extends(ChiFilterInterface) :: LumpedMassMatrixChiFilter

      real(rp), allocatable :: PWLambda(:), PWChi(:)
      class(FiniteElement), pointer :: eClosedRule

      real(rp) :: GaussLambda(1,50), GaussChi(1,50)

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

contains

   subroutine Initialize(a,plcd)
      class(LumpedMassMatrixChiFilter) :: a
      class(PLCDProblem), target :: plcd
      integer(ip) :: npoin

      a%plcd => plcd

      call a%plcd%Mesh%Getnpoin(npoin)
      call a%plcd%Memor%alloc(npoin,a%PWLambda,'PWLambda','PW')
      call a%plcd%Memor%alloc(npoin,a%PWChi,'PWChi','PW')

      call a%plcd%Mesh%ElementAlloc(a%eClosedRule,a%plcd%Memor,'ForceClosedRule','SIMP_TopOpt')

   end subroutine

   subroutine Finalize(a)
      class(LumpedMassMatrixChiFilter) :: a

      integer(ip) :: npoin


      call a%plcd%Mesh%GetNpoin(npoin)
      call a%plcd%Memor%dealloc(npoin,a%PWLambda,'PWLambda','PW')
      call a%plcd%Memor%dealloc(npoin,a%PWChi,'PWChi','PW')

      call a%plcd%Mesh%ElementDeAlloc(a%eClosedRule,a%plcd%Memor,'ForceClosedRule','SIMP_TopOpt')

   end subroutine

   subroutine AddGaussPointContributionToStrainEnergy(a,ielem,e,StrainEnergy)
      class(LumpedMassMatrixChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: StrainEnergy

      a%GaussLambda(1,e%igaus) = StrainEnergy


   end subroutine

   subroutine PostGaussStrainEnergy(a,ielem,e)
      class(LumpedMassMatrixChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e

      real(rp) :: NodeLambda(1,e%pnode)
      real(rp) :: elStrainEnergy(e%pnode),interp_lambda(1),StrainEnergy
      real(rp) :: dvol
      integer(ip) :: igaus

      call e%GaussToNodes(a%GaussLambda,NodeLambda)

      !Ensure that everything makes sense
      where (NodeLambda<0.0_rp) NodeLambda = 0.0_rp
      !Closed Rule of Integration, for assembling smoothed GaussChi and Lambda
      call a%plcd%Mesh%ElementLoad(ielem,a%eClosedRule)
      do igaus = 1,a%eClosedRule%pgaus
         a%eClosedRule%igaus = igaus

         call a%eClosedRule%elmder
         dvol = a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%detjm

         call a%eClosedRule%interpg(1_ip,NodeLambda,interp_lambda)
         StrainEnergy = interp_lambda(1)

         elStrainEnergy(1:a%eClosedRule%pnode) = a%eClosedRule%detjm*a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%shape(1:a%eClosedRule%pnode,a%eClosedRule%igaus)*StrainEnergy

         call a%plcd%Mesh%AssemblyToArray(a%eClosedRule,1_ip,elStrainEnergy,a%PWLambda)
      enddo
   end subroutine



   subroutine FinalizeStrainEnergyComputations(a)
      class(LumpedMassMatrixChiFilter) :: a

      call a%plcd%Mesh%Smooth(1,a%PWLambda)

      call a%plcd%FilePostpr%postpr(a%PWLambda,'LambdaPointWise'//adjustl(trim(int2str(a%plcd%itera))),a%plcd%istep,a%plcd%ctime,a%plcd%Mesh)
   end subroutine


   subroutine ResetChi(a)
      class(LumpedMassMatrixChiFilter) :: a

      a%PWChi = 0.0_rp
   end subroutine

   subroutine AddGaussPointContributionToChi(a,ielem,e,Chi)
      class(LumpedMassMatrixChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: Chi

      a%GaussChi(1,e%igaus) = Chi
   end subroutine

   subroutine PostGaussChi(a,ielem,e)
      class(LumpedMassMatrixChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e

      real(rp) :: NodeChi(1,e%pnode)
      real(rp) :: elChi(e%pnode), interp_chi(1),chi,dvol
      integer(ip) :: igaus

      call e%GaussToNodes(a%GaussChi,NodeChi)

      !Ensure that everything makes sense
      where (NodeChi<a%plcd%SIMPData%chimin) NodeChi = a%plcd%SIMPData%chimin
      where (NodeChi>a%plcd%SIMPData%chimax) NodeChi= a%plcd%SIMPData%chimax

      !Closed Rule of Integration, for assembling smoothed GaussChi and Lambda
      call a%plcd%Mesh%ElementLoad(ielem,a%eClosedRule)
      do igaus = 1,a%eClosedRule%pgaus
         a%eClosedRule%igaus = igaus

         call a%eClosedRule%elmder
         dvol = a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%detjm

         call a%eClosedRule%interpg(1_ip,NodeChi,interp_chi)
         chi = interp_chi(1)

         elChi(1:e%pnode) = a%eClosedRule%detjm*a%eClosedRule%weigp(a%eClosedRule%igaus)*a%eClosedRule%shape(1:a%eClosedRule%pnode,a%eClosedRule%igaus)*chi
         call a%plcd%Mesh%AssemblyToArray(a%eClosedRule,1_ip,elChi,a%PWChi)
      enddo
   end subroutine




   subroutine FinalizeChiComputations(a)
      class(LumpedMassMatrixChiFilter) :: a

      call a%plcd%Mesh%Smooth(1,a%PWChi)
   end subroutine

   subroutine UpdateChi(a,lambda)
      class(LumpedMassMatrixChiFilter) :: a
      real(rp) :: lambda

      integer(ip) :: npoin,ipoin
      real(rp) :: chi,BK
      real(rp) :: newchi,compareValue

      call a%plcd%Mesh%Getnpoin(npoin)
      do ipoin = 1,npoin
         BK = a%PWLambda(ipoin)/lambda
         chi = a%PWChi(ipoin)

         call UpdateChiValue(a,BK,chi,newchi)
         a%PWChi(ipoin) = newchi
      enddo


   end subroutine

   subroutine GetGaussPointChi(a,ielem,e,Chi)
      class(LumpedMassMatrixChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: Chi

      real(rp) :: auxChi(1)

      real(rp) :: elchi(e%pnode)

      call e%gather(1_ip,elchi,a%PWChi)
      !call e%interpg(1_ip,elchi,auxchi)
      call e%interpc(1_ip,elchi,auxchi)
      chi = auxchi(1)
   end subroutine

   subroutine GetCenterPointChi(a,ielem,e,Chi)
      class(LumpedMassMatrixChiFilter) :: a
      integer(ip) :: ielem
      class(FiniteElement) :: e
      real(rp) :: Chi

      real(rp) :: auxChi(1)

      real(rp) :: elchi(e%pnode)

      call e%gather(1_ip,elchi,a%PWChi)
      !call e%interpg(1_ip,elchi,auxchi)
      call e%interpc(1_ip,elchi,auxchi)
      chi = auxchi(1)
   end subroutine


   subroutine PostprocessChi(a,FilePostpr)
      class(LumpedMassMatrixChiFilter) :: a
      class(PostprFile) :: FilePostpr

      call a%plcd%FilePostpr%postpr(a%PWChi,'ChiPointWise'//adjustl(trim(int2str(a%plcd%itera))),a%plcd%istep,a%plcd%ctime,a%plcd%Mesh)
   end subroutine

   subroutine GetNodalChi(a,NodalChi)
      class(LumpedMassMatrixChiFilter),target :: a
      real(rp), pointer :: NodalChi(:)

      NodalChi => a%PWChi
   end subroutine



end module
