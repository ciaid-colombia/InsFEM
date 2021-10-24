module Mod_lmn_BaseElmope
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_Element
   use Mod_LowMachElement
   use Mod_ConvectiveElement
   use Mod_LowMach
   use Mod_lmn_elmdir

   implicit none
   
   class(LowMachProblem), pointer :: a => NULL()
   
   character(6) :: itask

   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: ExternalForces => NULL()
      procedure(), NOPASS, pointer :: ComputeAdvectionVelocity => NULL()
      procedure(), NOPASS, pointer :: ComputeTemperature => NULL()
      procedure(), NOPASS, pointer :: ComputeDensity => NULL()
      procedure(), NOPASS, pointer :: TimeIntegrationToElext => NULL()
     
      !Additional pointers (for Elmope only)
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(lmn_elmbuv), NOPASS, pointer :: lmn_elmbuv => NULL()
      procedure(lmn_elmrhu), NOPASS, pointer :: lmn_elmrhu => NULL()
      procedure(lmn_elmrhp), NOPASS, pointer :: lmn_elmrhp => NULL()
      procedure(lmn_elmbuq), NOPASS, pointer :: lmn_elmbuq => NULL()
      procedure(lmn_elmbtw), NOPASS, pointer :: lmn_elmbtw => NULL()
      procedure(lmn_elmrht), NOPASS, pointer :: lmn_elmrht => NULL()
   end type
   type(PPointer) :: ProcPointer
      
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: AllocateArrays => NULL()
      procedure(), NOPASS, pointer :: PreLoop => NULL()
      procedure(), NOPASS, pointer :: PreDirichlet => NULL()
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: ElmatsToZero => NULL()
      procedure(), NOPASS, pointer :: Gathers => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: Interpolates => NULL()
      procedure(), NOPASS, pointer :: ComputeTaus => NULL()
      procedure(), NOPASS, pointer :: ComputeTaus_seg => NULL()
      procedure(), NOPASS, pointer :: ComputeTestf => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats_seg => NULL()
      procedure(), NOPASS, pointer :: InGaussElmatsAssembly => NULL()
      procedure(), NOPASS, pointer :: ComputeResidual => NULL()
      procedure(), NOPASS, pointer :: AssemblyEndite => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PostLoop => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()   
      procedure(), NOPASS, pointer :: Assembly => NULL()
      procedure(), NOPASS, pointer :: DeallocateArrays => NULL()
   end type
   type(PHook) :: ProcHook


   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ielem,nelem
  
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv

   integer(ip)           :: ielty0 = 0 !Previous element type
      
   real(rp), allocatable :: elvel(:,:,:),elpre(:,:),eltem(:,:),elsou(:)
   real(rp), allocatable :: elext_mom(:),elext_pre(:)
   real(rp), allocatable :: AGradN(:)
   real(rp), allocatable :: gpvel(:,:),gpden(:),acpth(:),gpadv(:),gppre(:),gptem(:),gpsou(:)
      
   real(rp) :: chale(2),timom,ticon,tiene
   real(rp) :: acden,acvis,actco,acsph,actex,acrea,acsou!,acrcp
   real(rp) :: gpvno,gtemp,divvel,elext_ene
   real(rp), allocatable :: grpre(:,:),grvel(:,:),grtem(:,:)
   
   integer(ip) :: igaus
   integer(ip) :: kfl_GoIteInGauss
   
   !Test function
   real(rp), allocatable :: testf_mom(:,:),testf_ene(:),testf_mom_DSS(:,:),testf_ene_DSS(:)
   
   !Matrices
   real(rp), allocatable :: elmuv(:,:,:,:),elmpv(:,:,:,:),elmtv(:,:,:,:)
   real(rp), allocatable :: elmpq(:,:,:,:),elmuq(:,:,:,:)
   real(rp), allocatable :: elmtw(:,:,:,:), elmuw(:,:,:,:)
   real(rp), allocatable :: elrhu(:,:),elrhp(:,:),elrht(:,:)
   real(rp), allocatable :: wrmat1(:,:),wrmat2(:,:)
   
   real(rp), allocatable :: elmuvJ(:,:,:,:),elmuwJ(:,:,:,:),elmtvJ(:,:,:,:),elmuqJ(:,:,:,:)
   
   real(rp) :: dvolt0,dvolt1
   
   !Reference Dtinv
   !This one is for things which have their own integration scheme (taus, dynamic subscales)
   !It is necessary only for CrankNicolson
   real(rp) :: ReferenceDtinv
   
   !Mass Matrix
   real(rp), allocatable :: elmat_mass(:,:,:,:),elrhs_mass(:,:)
   real(rp), allocatable :: elm_mass(:,:)
   
   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   subroutine SetPointersAndHooksToNULLSUB
      implicit none
      !External Procedures
      procedure() :: NULLSUB
   
      !Pointers
      ProcPointer%ExternalForces => NULLSUB
      ProcPointer%ComputeAdvectionVelocity => NULLSUB
      ProcPointer%ComputeTemperature => NULLSUB
      ProcPointer%ComputeDensity => NULLSUB
      ProcPointer%TimeIntegrationToElext => NULLSUB
   
      !Hooks
      ProcHook%AllocateArrays => NULLSUB
      ProcHook%PreLoop => NULLSUB
      ProcHook%PreDirichlet => NULLSUB
      ProcHook%Initializations => NULLSUB
      ProcHook%OnIeltyChange => NULLSUB
      ProcHook%PreGauss => NULLSUB
      ProcHook%ElmatsToZero => NULLSUB
      ProcHook%Gathers => NULLSUB
      ProcHook%InGauss => NULLSUB
      ProcHook%Interpolates => NULLSUB
      ProcHook%ComputeTaus => NULLSUB
      ProcHook%ComputeTaus_seg => NULLSUB
      ProcHook%ComputeTestf => NULLSUB
      ProcHook%InGaussElmats => NULLSUB
      ProcHook%InGaussElmats_seg => NULLSUB
      ProcHook%InGaussElmatsAssembly => NULLSUB
      ProcHook%ComputeResidual => NULLSUB
      ProcHook%AssemblyEndite => NULLSUB
      ProcHook%Finalizations => NULLSUB
      ProcHook%PostLoop => NULLSUB
      ProcHook%PhysicalProp => NULLSUB
      ProcHook%Assembly => NULLSUB
      ProcHook%DeallocateArrays => NULLSUB
   end subroutine
   
   subroutine AllocateBaseElmopeArrays
      !Set Time Integrator
      call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
      ReferenceDtinv = a%dtinv
      !If CrankNicolson the reference is 1/2 dt, used in dynamic subscales, taus etc
      if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = 2*a%dtinv
      
      !Other arrays alloc
      call a%Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','lmn_elmope')
      call a%Memor%alloc(        e%mnode,a%ncomp-1,eltem,'eltem','lmn_elmope')
      call a%Memor%alloc(        e%mnode,a%ncomp-1,elpre,'elpre','lmn_elmope')
      call a%Memor%alloc(        e%mnode,          elsou,'elsou','lmn_elmope')
      call a%Memor%alloc(e%ndime,elext_mom,'elext_mom','lmn_elmope')
      call a%Memor%alloc(e%ndime,elext_pre,'elext_pre','lmn_elmope')
      call a%Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','lmn_elmope')
      call a%Memor%alloc(        a%ncomp-1,gpden,'gpden','lmn_elmope')
      call a%Memor%alloc(        1,gpsou,'gpsou','lmn_elmope')
      call a%Memor%alloc(        a%ncomp-1,gptem,'gptem','lmn_elmope')
      call a%Memor%alloc(          a%ncomp,acpth,'acpth','lmn_elmope')
      call a%Memor%alloc(        a%ncomp-1,gppre,'gppre','lmn_elmope')
      call a%Memor%alloc(e%ndime,gpadv,'gpadv','lmn_elmope')
      call a%Memor%alloc(e%mnode,AGradN,'AGradN','lmn_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','lmn_EndElmope')
      call a%Memor%alloc(1_ip   ,e%ndime,grpre,'grpre','lmn_EndElmope')
      call a%Memor%alloc(1_ip   ,e%ndime,grtem,'grtem','lmn_EndElmope')
      
      !Physical Parameters
      call a%GetPhysicalParameters(acvis=acvis,acsph=acsph,actco=actco,actex=actex,acpth=acpth,acrea=acrea,acsou=acsou)
      actco = actco/acsph
      actex = actex/acsph
!      acrcp = acrea/acsph 
   end subroutine
   
   subroutine AllocateElmopeArrays
      !Matrices Alloc
      call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','lmn_elmope')
      call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','lmn_elmope')
      call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','lmn_elmope')
      call a%Memor%alloc(e%mnode,e%mnode,wrmat2,'wrmat2','lmn_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,1,e%mnode,elmtw,'elmtw','lmn_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','lmn_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmtv,'elmtv','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','lmn_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrhu,'elrhu','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,elrhp,'elrhp','lmn_elmope')
      call a%Memor%alloc(1,e%mnode,elrht,'elrht','lmn_elmope')
   end subroutine
   
   subroutine ElementGathers
      implicit none
      integer(ip) :: itime
      
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      call e%gather(1,eltem(:,1),a%tempe(:,1))
      call e%gather(1,eltem(:,2),a%tempe(:,3))
      call e%gather(1,elpre(:,1),a%press(:,1))
      call e%gather(1,elpre(:,2),a%press(:,3))
      if(a%kfl_sourc==2) call e%gather(1,elsou,a%PointwiseSource)
      do itime = 3,nsteps !Time bdf2 and others
         call e%gather(e%ndime,elvel(:,:,itime),a%veloc(:,:,itime+1)) 
         call e%gather(1,eltem(:,itime),a%tempe(:,itime+1))
         call e%gather(1,elpre(:,itime),a%press(:,itime+1))
      enddo 
   end subroutine  
   
   subroutine ElementInterpolates
      implicit none
      integer(ip) :: itime
      
      call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1))
      call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2))
      call e%interpg(1,eltem(:,1),gptem(1))
      call e%interpg(1,eltem(:,2),gptem(2))
      call e%interpg(1,elpre(:,1),gppre(1))
      call e%interpg(1,elpre(:,2),gppre(2))
      if(a%kfl_sourc==2) call e%interpg(1,elsou,gpsou)
      do itime = 3,nsteps !Time bdf2 and others
         call e%interpg(e%ndime,elvel(:,:,itime),gpvel(:,itime))
         call e%interpg(1,eltem(:,itime),gptem(itime))
         call e%interpg(1,elpre(:,itime),gppre(itime))
      enddo
   end subroutine
  
   subroutine lmn_TimeIntegrationToElext 
         call lmn_TimeIntegrationToElext_mom(e,Integrator,acden,a%dtinv,gpvel,elext_mom)
         call lmn_TimeIntegrationToElext_ene(e,Integrator,acden,a%dtinv,gptem,elext_ene)
   end subroutine

   subroutine DeallocateElmopeArrays
      !Matrix Deallocations
      call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','lmn_elmope')
      call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','lmn_elmope')
      call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','lmn_elmope')
      call a%Memor%dealloc(e%mnode,e%mnode,wrmat2,'wrmat2','lmn_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmtw,'elmtw','lmn_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','lmn_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmtv,'elmtv','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','lmn_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrhu,'elrhu','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,elrhp,'elrhp','lmn_elmope')
      call a%Memor%dealloc(1,e%mnode,elrht,'elrht','lmn_elmope')
   end subroutine

   subroutine DeallocateBaseElmopeArrays
      !Other arrays alloc
      call a%Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','lmn_elmope')
      call a%Memor%dealloc(        e%mnode,a%ncomp-1,eltem,'eltem','lmn_elmope')
      call a%Memor%dealloc(        e%mnode,a%ncomp-1,elpre,'elpre','lmn_elmope')
      call a%Memor%dealloc(        e%mnode,          elsou,'elsou','lmn_elmope')
      call a%Memor%dealloc(e%ndime,elext_mom,'elext_mom','lmn_elmope')
      call a%Memor%dealloc(e%ndime,elext_pre,'elext_pre','lmn_elmope')
      call a%Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','lmn_elmope')
      call a%Memor%dealloc(        a%ncomp-1,gpden,'gpden','lmn_elmope')
      call a%Memor%dealloc(        1,gpsou,'gpsou','lmn_elmope')
      call a%Memor%dealloc(        a%ncomp-1,gptem,'gptem','lmn_elmope')
      call a%Memor%dealloc(          a%ncomp,acpth,'acpth','lmn_elmope')
      call a%Memor%dealloc(        a%ncomp-1,gppre,'gppre','lmn_elmope')
      call a%Memor%dealloc(e%ndime,gpadv,'gpadv','lmn_elmope')
      call a%Memor%dealloc(e%mnode,AGradN,'AGradN','lmn_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','lmn_EndElmope')
      call a%Memor%dealloc(1_ip   ,e%ndime,grpre,'grpre','lmn_EndElmope')
      call a%Memor%dealloc(1_ip   ,e%ndime,grtem,'grtem','lmn_EndElmope')
   end subroutine
end module









