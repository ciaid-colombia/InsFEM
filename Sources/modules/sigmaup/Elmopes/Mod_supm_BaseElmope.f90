module Mod_supm_BaseElmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_ThreeField
   use Mod_Element
   use Mod_NavierStokesElement
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_nsm_elmdir
   use Mod_ThreeFieldElement
   use Mod_LogarithmicElement
   use Mod_SupExacso    
   use Mod_SupOperations
   implicit none
    
   class(ThreeFieldNSProblem), pointer :: a  
   character(6) :: itask
   type(SupExacso) :: exacso
   procedure() :: NULLSUB
   
   !the elemental matrix take the form
   !************************************************
   !elmbst  elmbut  elmbpt  S    elmrhc
   !elmbsv  elmbuv  elmbpv  u =  elmrhu
   !elmbsq  elmbuq  elmbpq  p    elmrhp
   !************************************************
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer ::  PostGaussElmats_sup
      procedure(nsm_elmbuv), NOPASS, pointer ::  supm_elmbuv
      procedure(nsm_elmrhu), NOPASS, pointer ::  supm_elmrhu
      procedure(nsm_elmrhp), NOPASS, pointer ::  supm_elmrhp
      procedure(nsm_elmbuq), NOPASS, pointer ::  supm_elmbuq
      procedure(nsm_elmbpv), NOPASS, pointer ::  supm_elmbpv
      procedure(), NOPASS, pointer ::  ConstitutiveComponents   
      procedure(), NOPASS, pointer ::  ExternalForces_sup
      procedure(), NOPASS, pointer ::  ComputeAdvectionVelocity_sup  
      procedure(), NOPASS, pointer ::  ResProRHS 
      procedure(), NOPASS, pointer ::  TauConstitutive
      procedure(), NOPASS, pointer ::  TemporalDerivatives
      procedure(), NOPASS, pointer ::  TwoFieldTerms
      procedure(), NOPASS, pointer ::  ViscoelasticGalerkin
      procedure(), NOPASS, pointer ::  ViscoelasticEstab1
      procedure(), NOPASS, pointer ::  ViscoelasticEstab1Lapla      
      procedure(), NOPASS, pointer ::  ViscoelasticEstab2
      procedure(), NOPASS, pointer ::  PenaltyTerm
      procedure(), NOPASS, pointer ::  VelocityProjection
      procedure(), NOPASS, pointer ::  DiscontinuityCapturing
      procedure(), NOPASS, pointer ::  PreAssembly_sup
      procedure(), NOPASS, pointer ::  StokesOutside
      procedure(), NOPASS, pointer ::  FEResPro
      procedure(), NOPASS, pointer ::  ComputeAGradV_sup
      procedure(), NOPASS, pointer ::  ComputeAdvectionVelocityNorm_sup
      procedure(nsm_elmrhp_dss), NOPASS, pointer ::  supm_elmrhp_dss
      procedure(nsm_elmrhu_dss), NOPASS, pointer ::  supm_elmrhu_dss
   end type
   
   type(PPointer) :: ProcPointer
   
     !Hooks
      procedure(), pointer :: ProcHook_PreLoop           => NULLSUB
      procedure(), pointer :: ProcHook_PreAllocate       => NULLSUB
      procedure(), pointer :: ProcHook_Initializations   => NULLSUB
      procedure(), pointer :: ProcHook_OnIeltyChange     => NULLSUB
      procedure(), pointer :: ProcHook_PreGauss          => NULLSUB
      procedure(), pointer :: ProcHook_ElmatsToZero      => NULLSUB
      procedure(), pointer :: ProcHook_Gathers           => NULLSUB
      procedure(), pointer :: ProcHook_InGauss           => NULLSUB
      procedure(), pointer :: ProcHook_Interpolates      => NULLSUB
      procedure(), pointer :: ProcHook_PostInterpolate   => NULLSUB
      procedure(), pointer :: ProcHook_ComputeTaus       => NULLSUB
      procedure(), pointer :: ProcHook_ComputeTestf      => NULLSUB
      procedure(), pointer :: ProcHook_InGaussElmats     => NULLSUB
      procedure(), pointer :: ProcHook_InGaussElmatsAssembly => NULLSUB
      procedure(), pointer :: ProcHook_AssemblyEndite    => NULLSUB
      procedure(), pointer :: ProcHook_Finalizations     => NULLSUB
      procedure(), pointer :: ProcHook_PostLoop          => NULLSUB
      procedure(), pointer :: ProcHook_PhysicalProp      => NULLSUB
      procedure(), pointer :: ProcHook_PreDirichlet      => NULLSUB
      procedure(), pointer :: ProcHook_Elext             => NULLSUB
      procedure(), pointer :: ProcHook_GradientCapturing => NULLSUB
      procedure(), pointer :: ProcHook_HessianComponents => NULLSUB
      
   class(FiniteElement) , pointer     :: e => NULL()
   integer(ip)             :: ielem, nelem

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: nsteps
   real(rp)                :: LHSdtinv
   integer(ip)             :: ielty0 = 0 !Previous element type    
   real(rp), allocatable   :: elmat(:,:,:,:)
   real(rp), allocatable   :: elrhs(:,:)   
   integer(ip)             :: igaus
   real(rp)                :: dvol
   real(rp), allocatable   :: elvel(:,:,:), elsig(:,:,:), bovel(:,:)
   real(rp), allocatable   :: elpre(:,:)
   real(rp), allocatable   :: elext(:),elextC(:),elextS(:), elextEstab(:), elextEstab2(:),  elextEstab3(:), elextEstab4(:), elextEstab5(:)
   real(rp), allocatable   :: elextSEstab(:), elextSEstab2(:), elextSEstab3(:), elextSEstab4(:), elextSMat(:,:), elextSEstabMat(:,:)
   real(rp), allocatable   :: AGradV(:), testf(:), AGradV_VE(:) 
   real(rp), allocatable   :: gpvel(:,:),gpadvec(:),gpsig(:,:), gpadvec_VE(:)
   real(rp)                :: chale(2),timom,tidiv, dvolt0,dvolt1,dvolt2,dvolt3,gprhs(3)
   integer(ip)             :: nmean,auxiter !Stabilization
   integer(ip)             :: iboun
   real(rp)                :: dsurf, acden,acvis,beta,auxVE,auxG,lambda, auxL,lambda0
   real(rp)                :: reyno
   real(rp)                :: gpvno,gppre(1,1),divvel, gpvno_VE
   real(rp), parameter     :: zensi = 0.0_rp   
   integer(ip)             :: ipnode0 = 0 !no change pnode
   real(rp), allocatable   :: elmuv(:,:,:,:), elmpq(:,:,:,:), elmpv(:,:,:,:),elmuq(:,:,:,:)
   real(rp), allocatable   :: elmst(:,:,:,:), elmut(:,:,:,:), elmpt(:,:,:,:)
   real(rp), allocatable   :: elmsv(:,:,:,:), elmsq(:,:,:,:), elmsv2(:,:,:,:)
   real(rp), allocatable   :: elrhu(:,:), elrhp(:,:), elrhc(:,:),elrhu2(:,:)
   real(rp), allocatable   :: wrmat1(:,:), wrmat2(:,:,:,:)
   real(rp)                :: tisig, tisig_static
   !Residual Projections   
   real(rp), allocatable   :: elrep(:,:),gprep(:)
   real(rp), allocatable   :: gpres(:),elres(:,:)
   real(rp), allocatable   :: grpre(:,:),grvel(:,:),grsig(:,:)
   integer(ip)             :: ibopo,npoin
   real(rp), pointer       :: exnor(:,:)
   !Oss Split
   real(rp), allocatable   :: gpdivs(:),elrepdivs(:,:)
   real(rp), allocatable   :: gpconv(:),elrepconv(:,:)
   real(rp), allocatable   :: gpgrap(:),elrepgrap(:,:) 
   real(rp), allocatable   :: elresdivs(:,:), elresconv(:,:), elresgrap(:,:)  , elreslapl(:,:)
   !Dynamic S-OSS
   real(rp), allocatable   :: gpdivu(:),elrepdivu(:,:)
   real(rp), allocatable   :: gpgrau(:),elrepgrau(:,:)
   !Oss Split LCR
   real(rp), allocatable   :: gplapl(:),elreplapl(:,:)  
   real(rp), allocatable   :: gpgrav(:), elrepgrav(:,:)
   real(rp), allocatable   :: gpconvsigma(:), elrepconvsigma(:,:)
   real(rp), allocatable   :: gpdeform(:), elrepdeform(:,:)
   real(rp), allocatable   :: gpexp(:), elrepexp(:,:)
   !Discontinuity Capturing
   real(rp), allocatable   :: gpgrad(:),elrepGrad(:,:)
   real(rp), allocatable   :: grsigRP(:,:),grsigRPO(:,:)
   !Gradient Projection
   real(rp), allocatable   :: elrepSGrad(:,:),gpsgrad(:)
   integer(ip)             :: idime,itime,auxtens,auxGrad,auxSGrad,auxdim,auxrepro
   
   real(rp), allocatable   :: grvelRP(:,:)
   integer(ip)             :: auxSkip   
   integer(ip)             :: jdime,auxntens,auxndim,kdime,ldime
   !Level Set
   integer(ip)             :: inode,ipoin
   integer(ip)             :: elemStatus,ngauss_minus,ngauss_plus,ngaus_total,poinStatus
   real(rp), allocatable   :: weigp(:),xloc(:,:)
   !Gravity forces
   real(rp), allocatable   :: elext3(:)
   real(rp)                :: elext2(3)
   !Physical properties
   real(rp)                :: auxvis
   !Number of materials
   integer(ip)             :: imat=1 
   !Stabilization
   integer(ip)             :: auxpena  
   !shock-capturing parameters
   real(rp)                :: cshock(2) 
   real(rp)                :: auxp1,auxp2,auxp3,auxp4,kdisc,facdisc
   !Free Surface Stokes Outside
   real(rp), allocatable   :: elmuqFS(:,:,:,:),elmpvFS(:,:,:,:),auxtestf(:),elmstFS(:,:,:,:),elmutFS(:,:,:,:),elmuvFS(:,:,:,:)
   !Logarithm conformation Reformulation
   real(rp), allocatable   :: ExpGiesekus_Matrix(:,:)
   real(rp), allocatable   :: ExpGpPsi_Matrix(:,:), ExpGpPsi(:)
   real(rp), allocatable   :: DivExpPsi(:),  GrPsiMatrix(:,:,:), GrExpPsi(:,:)
   real(rp), allocatable   :: ExpGradV(:,:), GrExpMatrix(:,:,:), RHSconst(:,:),RHSconstEstab(:,:), RHSmom(:)
   real(rp), allocatable   :: RHS1momGal(:), RHS2momGal(:,:), RHSGiesekus(:,:)
   real(rp), allocatable   :: RHSconst_conv(:,:), RHSconst_deform(:,:), ConvExpMatrix(:,:), GrvelExp(:,:)
   !Vorticity
   real(rp), allocatable   :: vorti(:,:)
   
   !Reference Dtinv
   !This one is for things which have their own integration scheme (taus, dynamic subscales)
   !It is necessary only for CrankNicolson
   real(rp) :: ReferenceDtinv

#include "COMPOSEPROCS_POINTERS_50.i90"      
contains
#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   subroutine SetPointersAndHooksToNULL
      use typre
      implicit none
      procedure() :: NULLSUB   
      !Pointers
      ProcPointer%PostGaussElmats_sup => NULLSUB 
      ProcPointer%supm_elmbuv => nsm_elmbuv
      ProcPointer%supm_elmrhu => nsm_elmrhu
      ProcPointer%supm_elmrhp => nsm_elmrhp
      ProcPointer%supm_elmbuq => nsm_elmbuq
      ProcPointer%supm_elmbpv => nsm_elmbpv
      
      ProcPointer%ConstitutiveComponents => NULLSUB 

      ProcPointer%ResProRHS => NULLSUB
      ProcPointer%TauConstitutive => NULLSUB
      ProcPointer%TwoFieldTerms => NULLSUB
      ProcPointer%ViscoelasticGalerkin => NULLSUB
      ProcPointer%ViscoelasticEstab1  => NULLSUB
      ProcPointer%ViscoelasticEstab2 => NULLSUB

      ProcPointer%VelocityProjection => NULLSUB
      ProcPointer%DiscontinuityCapturing => NULLSUB
      ProcPointer%PreAssembly_sup => NULLSUB
      ProcPointer%FEResPro => NULLSUB
      ProcPointer%TemporalDerivatives => NULLSUB
      
      ProcPointer%ViscoelasticEstab1Lapla => NULLSUB      
      ProcPointer%PenaltyTerm  => NULLSUB
      ProcPointer%StokesOutside => NULLSUB
      ProcPointer%ExternalForces_sup => NULLSUB
      ProcPointer%ComputeAdvectionVelocity_sup => NULLSUB
      ProcPointer%ComputeAGradV_sup => NULLSUB
      ProcPointer%ComputeAdvectionVelocityNorm_sup => NULLSUB
      
      ProcPointer%supm_elmrhp_dss => nsm_elmrhp_dss
      ProcPointer%supm_elmrhu_dss => nsm_elmrhu_dss
      
   call SetHooksToNULL       
   end subroutine
   
   
   subroutine  SetHooksToNULL
      implicit none
      procedure() :: NULLSUB
      !Hooks
      ProcHook_Initializations   => NULLSUB
      ProcHook_Gathers           => NULLSUB
      ProcHook_OnIeltyChange     => NULLSUB
      ProcHook_PreGauss          => NULLSUB
      ProcHook_Interpolates      => NULLSUB
      ProcHook_InGauss           => NULLSUB
      ProcHook_InGaussElmats     => NULLSUB
      ProcHook_Finalizations     => NULLSUB
      ProcHook_PhysicalProp      => NULLSUB
      ProcHook_HessianComponents => NULLSUB  
      ProcHook_ComputeTestf      => NULLSUB
      ProcHook_PreAllocate       => NULLSUB
      ProcHook_PostInterpolate   => NULLSUB
      ProcHook_PreLoop           => NULLSUB
      ProcHook_ElmatsToZero      => NULLSUB
      ProcHook_Elext             => NULLSUB
      ProcHook_AssemblyEndite    => NULLSUB
      ProcHook_PostLoop          => NULLSUB
      ProcHook_GradientCapturing => NULLSUB
   end subroutine 
   
   
   subroutine AllocateBaseElmopeArrays
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
      ReferenceDtinv = a%dtinv
      
      !Constitutive case
      call a%Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','supm_EnditeElmope')
      call a%Memor%alloc(      e%mnode,a%ncomp-1,elpre,'elpre','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supm_EnditeElmope')   
      call a%Memor%alloc(e%ndime,elext,'elext','supm_EnditeElmope')
      call a%Memor%alloc(1,elextC,'elextC','supm_EnditeElmope')  
      call a%Memor%alloc(auxtens,elextS,'elextS','supm_EnditeElmope')  

      call a%Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,a%ncomp-1,gpsig,'gpsig','supm_EnditeElmope')   
      call a%Memor%alloc(e%ndime,gpadvec,'gpadvec','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpadvec_VE,'gpadvec_VE','supm_EnditeElmope')
      call a%Memor%alloc(e%mnode,AGradV,'AGradV','supm_EnditeElmope')
      call a%Memor%alloc(e%mnode,AGradV_VE,'AGradV_VE','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,e%ndime,grsig,'grsig','supm_EnditeElmope')
      
      !Allocate elextEstab arrays - only analytical solution
      call a%Memor%alloc(e%ndime,e%ndime,elextSMat,'elextSMat','supm_EnditeElmope') 
      call a%Memor%alloc(e%ndime,e%ndime,elextSEstabMat,'elextSEstabMat','supm_EnditeElmope') 
      call a%Memor%alloc(auxtens,elextSEstab,'elextSEstab','supm_EnditeElmope') 
      call a%Memor%alloc(auxtens,elextSEstab2,'elextSEstab2','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,elextSEstab3,'elextSEstab3','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,elextSEstab4,'elextSEstab4','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,elextEstab,'elextEstab','supm_EnditeElmope') 
      call a%Memor%alloc(e%ndime,elextEstab2,'elextEstab2','supm_EnditeElmope') 
      call a%Memor%alloc(e%ndime,elextEstab3,'elextEstab3','supm_EnditeElmope') 
      call a%Memor%alloc(e%ndime,elextEstab4,'elextEstab4','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,elextEstab5,'elextEstab5','supm_EnditeElmope')
   end subroutine 

   subroutine AllocateBaseElmopeMatrices
      call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supm_elmope')
      call a%Memor%alloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','supm_elmope')
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrhu,'elrhu','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrhu2,'elrhu2','supm_elmope')
      call a%Memor%alloc(1,e%mnode,elrhp,'elrhp','supm_elmope')
   
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      call a%Memor%alloc(auxtens,e%mnode,auxtens,e%mnode,elmst,'elmst','supm_elmope')
      call a%Memor%alloc(auxtens,e%mnode,e%ndime,e%mnode,elmut,'elmut','supm_elmope')
      call a%Memor%alloc(auxtens,e%mnode,1,e%mnode,elmpt,'elmpt','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,auxtens,e%mnode,elmsv,'elmsv','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,auxtens,e%mnode,elmsv2,'elmsv2','supm_elmope')
      call a%Memor%alloc(1,e%mnode,auxtens,e%mnode,elmsq,'elsq','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,wrmat2,'wrmat2','supm_elmope')
      call a%Memor%alloc(auxtens,e%mnode,elrhc,'elrhc','supm_elmope')

   end subroutine
   
   
   subroutine VelocityAndSigmaGathers
      implicit none
      integer(ip) :: itime
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      call e%gather(auxtens,elsig(:,:,1),a%sigma(:,:,1))
      call e%gather(auxtens,elsig(:,:,2),a%sigma(:,:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(e%ndime,elvel(:,:,itime),a%veloc(:,:,itime+1))
         call e%gather(auxtens,elsig(:,:,itime),a%sigma(:,:,itime+1))  
      enddo   
   end subroutine  

   subroutine InterpolateGpVelocityAndSigma
      implicit none
      integer(ip) :: itime
       call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1)) !i-1 iteration
       call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2)) !j-1 time
       call e%interpg(auxtens,elsig(:,:,1),gpsig(:,1)) !i-1 iteration
       call e%interpg(auxtens,elsig(:,:,2),gpsig(:,2)) !i-1 time         
       do itime = 3,nsteps ! Time bdf2 and others
         call e%interpg(e%ndime,elvel(:,:,itime),gpvel(:,itime))
         call e%interpg(auxtens,elsig(:,:,itime),gpsig(:,itime))
       enddo
   end subroutine
  
  
   subroutine DeallocateBaseElmopeMatrices
    !Matrix Deallocations
      call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supm_elmope')
      call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','supm_elmope')
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrhu,'elrhu','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrhu2,'elrhu2','supm_elmope')
      call a%Memor%dealloc(1,e%mnode,elrhp,'elrhp','supm_elmope')
   
      !Constitutive case
      call a%Memor%dealloc(auxtens,e%mnode,auxtens,e%mnode,elmst,'elmst','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,e%ndime,e%mnode,elmut,'elmut','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,1,e%mnode,elmpt,'elmpt','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,auxtens,e%mnode,elmsv,'elmsv','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,auxtens,e%mnode,elmsv2,'elmsv2','supm_elmope')
      call a%Memor%dealloc(1,e%mnode,auxtens,e%mnode,elmsq,'elsq','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,wrmat2,'wrmat2','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elrhc,'elrhc','supm_elmope') 
      
   end subroutine
   
   subroutine DeallocateBaseElmopeArrays
      call a%Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supm_EnditeElmope')   
      call a%Memor%dealloc(      e%mnode,a%ncomp-1,elpre,'elpre','supm_EnditeElmope')   
      call a%Memor%dealloc(e%ndime,elext,'elext','supm_EnditeElmope')
      call a%Memor%dealloc(1,elextC,'elextC','supm_EnditeElmope')  
      call a%Memor%dealloc(auxtens,elextS,'elextS','supm_EnditeElmope')  

      call a%Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,a%ncomp-1,gpsig,'gpsig','supm_EnditeElmope')   
      call a%Memor%dealloc(e%ndime,gpadvec,'gpadvec','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpadvec_VE,'gpadvec_VE','supm_EnditeElmope')
      call a%Memor%dealloc(e%mnode,AGradV,'AGradV','supm_EnditeElmope')
      call a%Memor%dealloc(e%mnode,AGradV_VE,'AGradV_VE','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,e%ndime,grsig,'grsig','supm_EnditeElmope') 
      
      !Deallocate elextEstab arrays - only analytical solution
      call a%Memor%dealloc(e%ndime,e%ndime,elextSMat,'elextSMat','supm_EnditeElmope') 
      call a%Memor%dealloc(e%ndime,e%ndime,elextSEstabMat,'elextSEstabMat','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,elextSEstab,'elextSEstab','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,elextSEstab2,'elextSEstab2','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,elextSEstab3,'elextSEstab3','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,elextSEstab4,'elextSEstab4','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,elextEstab,'elextEstab','supm_EnditeElmope') 
      call a%Memor%dealloc(e%ndime,elextEstab2,'elextEstab2','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,elextEstab3,'elextEstab3','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,elextEstab4,'elextEstab4','supm_EnditeElmope')
      call a%Memor%dealloc(auxtens,elextEstab5,'elextEstab5','supm_EnditeElmope')
   end subroutine
   
  
end module
   
