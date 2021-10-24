module Mod_nsc_pr_BaseElmope
   use Mod_Mesh
   use Mod_Memor
   use Mod_Element
   use Mod_php_SetTimeIntegrator
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_NSCompressibleElement
   use Mod_NSCompressiblePrimitive
   use Mod_NSCompressiblePrimitiveElement
   use Mod_NSCompressibleImplicitElement
   use Mod_NSCompressibleSubroutines
   use Mod_NscExacso   
   use Mod_NscDamping   
   implicit none
   
   class(NSCompressiblePrimitiveProblem), pointer :: a => NULL()
   
   character(6) :: itask

   !Pointers
      procedure(), pointer :: ProcPointer_nsc_pr_ExternalForces => NULL()
      procedure(), pointer :: ProcPointer_nsc_pr_ComputeVariables => NULL()
      procedure(), pointer :: ProcPointer_nsc_pr_ComputeTransientCoefficients => NULL()
      procedure(), pointer :: ProcPointer_nsc_pr_ComputeConvectionCoefficients => NULL()
      procedure(), pointer :: ProcPointer_nsc_pr_ComputeDiffusionCoefficients => NULL()
      procedure(), pointer :: ProcPointer_nsc_pr_ComputeTransportCoefficients => NULL()
      !Additional pointers (for Elmope only)
      !the elemental matrix takes the form
      !************************************************
      !q·  elmbdq  elmbmq  elmbeq  rho =  elmrhd
      !n·  elmbdn  elmbmn  elmben  mom =  elmrhm
      !g·  elmbdg  elmbmg  elmbeg  ene =  elmrhe
      !************************************************
         procedure(), pointer :: ProcPointer_nsc_pr_PreAssembly => NULL()
         procedure(), pointer :: ProcPointer_nsc_pr_PostGaussElmats => NULL()
         procedure(nsc_pr_elmbdq), pointer :: ProcPointer_nsc_pr_elmbdq => NULL()
         procedure(nsc_pr_elmbdn), pointer :: ProcPointer_nsc_pr_elmbdn => NULL()
         procedure(nsc_pr_elmbdg), pointer :: ProcPointer_nsc_pr_elmbdg => NULL()
         procedure(nsc_pr_elmbmq), pointer :: ProcPointer_nsc_pr_elmbmq => NULL()
         procedure(nsc_pr_elmbmn), pointer :: ProcPointer_nsc_pr_elmbmn => NULL()
         procedure(nsc_pr_elmbmg), pointer :: ProcPointer_nsc_pr_elmbmg => NULL()
         procedure(nsc_pr_elmbeq), pointer :: ProcPointer_nsc_pr_elmbeq => NULL()
         procedure(nsc_pr_elmben), pointer :: ProcPointer_nsc_pr_elmben => NULL()
         procedure(nsc_pr_elmbeg), pointer :: ProcPointer_nsc_pr_elmbeg => NULL()
         procedure(nsc_pr_elmrhd), pointer :: ProcPointer_nsc_pr_elmrhd => NULL()
         procedure(nsc_pr_elmrhm), pointer :: ProcPointer_nsc_pr_elmrhm => NULL()
         procedure(nsc_pr_elmrhe), pointer :: ProcPointer_nsc_pr_elmrhe => NULL()
      
   
   !Hooks
      procedure(), pointer :: ProcHook_nsc_pr_PreLoop => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_PreAllocate => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_Initializations => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_OnIeltyChange => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_PreGauss => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_ElmatsToZero => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_Gathers => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_InGauss => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_Interpolates => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_PostInterpolate => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_ComputeTaus => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_ComputeTestf => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_InGaussElmats => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_InGaussElmatsAssembly => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_AssemblyEndite => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_Finalizations => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_PostLoop => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_PhysicalProp     => NULL()
      procedure(), pointer :: ProcHook_nsc_pr_PreDirichlet => NULL()
                                       
                                       
   type(NscExacso) :: exacso           
                                       
   type(NscDamping) :: damping           

   class(FiniteElement) , pointer     :: e        => NULL()
   integer(ip) :: ielem,nelem
  
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv

   integer(ip)           :: ielty0 = 0  !Previous element type
   integer(ip)           :: ipnode0 = 0 !Previous pnode number
      
   real(rp), allocatable :: elpre(:,:)
   real(rp), allocatable :: elvel(:,:,:)
   real(rp), allocatable :: eltem(:,:)
   real(rp), allocatable :: eltemv(:)
   real(rp), allocatable :: elexv(:)
   real(rp), allocatable :: AGradV(:)
   real(rp), allocatable :: gppre(:)  
   real(rp), allocatable :: gpvel(:,:)
   real(rp), allocatable :: gptem(:)  
   real(rp), allocatable :: gpadv(:)
   real(rp), allocatable :: vgvel(:)
   real(rp), allocatable :: grpre(:), grvel(:,:), grtem(:)
   real(rp), allocatable :: hess(:,:,:), lapl(:) 
      
   real(rp)    :: chale(2),timom(3), tidiv
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: gpadp,gpadt
   real(rp)    :: gpvno,divvel
   real(rp)    :: vgpre,vgtem
   real(rp)    :: aux,acalpha,acbeta
   real(rp)    :: gpden,gpspd,gpvst
   real(rp)    :: eltemp(1),eltemt(1)
   real(rp)    :: elexp(1),elext(1)
   real(rp)    :: Erhd,Erhe
   
   integer(ip) :: igaus
   integer(ip) :: kfl_GoIteInGauss
   
   !Euler contribution Matrix
   real(rp), allocatable :: Edd(:),   Edm(:,:),   Ede(:)
   real(rp), allocatable :: Emd(:,:), Emm(:,:,:), Eme(:,:)
   real(rp), allocatable :: Eed(:),   Eem(:,:),   Eee(:)
   real(rp), allocatable :: Erhm(:)

   !Transient Matrix coefficients
   real(rp), allocatable :: Atdd(:),   Atde(:)
   real(rp), allocatable :: Atmd(:,:), Atmm(:,:,:), Atme(:,:)
   real(rp), allocatable :: Ated(:),   Atem(:,:),   Atee(:)

   !Convection Matrix coefficients
   real(rp), allocatable :: Add(:),   Adm(:,:),   Ade(:)
   real(rp), allocatable :: Amd(:,:), Amm(:,:,:), Ame(:,:)
   real(rp), allocatable :: Aed(:),   Aem(:,:),   Aee(:)

   !Diffusion Matrix coefficients
   real(rp), allocatable :: KGdd(:,:)  !Radial damping diffusion
   real(rp), allocatable :: KGmm(:,:,:,:)
   real(rp), allocatable :: KGem(:,:,:),   KGee(:,:)
   real(rp), allocatable :: Kmm(:,:,:)
   real(rp), allocatable :: Kem(:,:),   Kee(:)

   !Reaction Matrix coefficients
   real(rp), allocatable :: Sdd(:)
   real(rp), allocatable :: Smd(:,:), Smm(:,:,:)
   real(rp), allocatable :: Sed(:),   Sem(:,:),   See(:)

   !Test functions
   !Adjoint Matrix coefficients
   real(rp), allocatable :: LTdd(:),   LTdm(:,:),   LTde(:)
   real(rp), allocatable :: LTmd(:,:), LTmm(:,:,:), LTme(:,:)
   real(rp), allocatable :: LTed(:),   LTem(:,:),   LTee(:)

   !Jacobian Gradient Convection Matrix coefficients
   real(rp), allocatable :: JAdd(:),   JAdm(:,:) , JAde(:)
   real(rp), allocatable :: JAmd(:,:), JAmm(:,:,:), JAme(:,:)
   real(rp), allocatable :: JAed(:),   JAem(:,:),   JAee(:)
   
   !Matrices
   real(rp), allocatable :: elmdq(:,:),   elmmq(:,:,:),   elmeq(:,:)
   real(rp), allocatable :: elmdn(:,:,:), elmmn(:,:,:,:), elmen(:,:,:)
   real(rp), allocatable :: elmdg(:,:),   elmmg(:,:,:),   elmeg(:,:)
   real(rp), allocatable :: elrhd(:),     elrhm(:,:),     elrhe(:)
   
   real(rp) :: dvolt0, dvolt1, dvolt2
   
   !Reference Dtinv
   !This one is for things which have their own integration scheme (taus, dynamic subscales)
   !It is necessary only for CrankNicolson
   real(rp) :: ReferenceDtinv
   
   !Mass Matrix
   real(rp), allocatable :: elmat_mass(:,:,:,:),elrhs_mass(:,:)
   real(rp), allocatable :: elm_mass(:,:)
   
   !#$COMPOSEPROCS 100
#include "COMPOSEPROCS_POINTERS_100.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_100.i90"

  
  
   
   
   
   subroutine SetPointersAndHooksToNULLSUB
   implicit none
      !External Procedures
      procedure() :: NULLSUB
   
   !Pointers
      ProcPointer_nsc_pr_ExternalForces => NULLSUB
      ProcPointer_nsc_pr_ComputeVariables => NULLSUB
      ProcPointer_nsc_pr_ComputeTransientCoefficients => NULLSUB
      ProcPointer_nsc_pr_ComputeConvectionCoefficients => NULLSUB
      ProcPointer_nsc_pr_ComputeDiffusionCoefficients => NULLSUB
      ProcPointer_nsc_pr_ComputeTransportCoefficients => NULLSUB
   
   !Hooks
      ProcHook_nsc_pr_PreLoop => NULLSUB
      ProcHook_nsc_pr_PreAllocate => NULLSUB
      ProcHook_nsc_pr_Initializations => NULLSUB
      ProcHook_nsc_pr_OnIeltyChange => NULLSUB
      ProcHook_nsc_pr_PreGauss => NULLSUB
      ProcHook_nsc_pr_ElmatsToZero => NULLSUB
      ProcHook_nsc_pr_Gathers => NULLSUB
      ProcHook_nsc_pr_InGauss => NULLSUB
      ProcHook_nsc_pr_Interpolates => NULLSUB
      ProcHook_nsc_pr_PostInterpolate => NULLSUB
      ProcHook_nsc_pr_ComputeTaus => NULLSUB
      ProcHook_nsc_pr_ComputeTestf => NULLSUB
      ProcHook_nsc_pr_InGaussElmats => NULLSUB
      ProcHook_nsc_pr_InGaussElmatsAssembly => NULLSUB
      ProcHook_nsc_pr_AssemblyEndite => NULLSUB
      ProcHook_nsc_pr_Finalizations => NULLSUB
      ProcHook_nsc_pr_PostLoop => NULLSUB
      ProcHook_nsc_pr_PhysicalProp => NULLSUB
      ProcHook_nsc_pr_PreDirichlet => NULLSUB
   end subroutine
   
   subroutine AllocateBaseElmopeArrays
      implicit none
      
      type(MemoryMan), pointer :: Memor => NULL()
      
      !Set Time Integrator
      call php_SetTimeIntegrator(a,Integrator,LHSdtinv,nsteps)
      ReferenceDtinv = a%dtinv
      !If CrankNicolson the reference is 1/2 dt, used in dynamic subscales, taus etc
      if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = 2*a%dtinv
      
      !Other arrays alloc
      Memor => a%Memor
      call Memor%alloc(e%mnode,a%ncomp-1,elpre,'elpre','nsc_pr_elmope')
      call Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','nsc_pr_elmope')
      call Memor%alloc(e%mnode,a%ncomp-1,eltem,'eltem','nsc_pr_elmope')
      call Memor%alloc(a%ncomp-1,gppre,'gppre','nsc_pr_elmope')
      call Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','nsc_pr_elmope')
      call Memor%alloc(a%ncomp-1,gptem,'gptem','nsc_pr_elmope')
      call Memor%alloc(e%ndime,grpre,'grpre','nsc_pr_elmope')
      call Memor%alloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_elmope')
      call Memor%alloc(e%ndime,grtem,'grtem','nsc_pr_elmope')
      call Memor%alloc(e%ndime,elexv,'elexv','nsc_pr_elmope')
      call Memor%alloc(e%ndime,eltemv,'eltemv','nsc_pr_elmope')
      call Memor%alloc(e%ndime,gpadv,'gpadv','nsc_pr_elmope')
      call Memor%alloc(e%mnode,AGradV,'AGradV','nsc_pr_elmope')
      call Memor%alloc(e%ndime,vgvel,'vgvel','nsc_pr_elmope')

      !Physical Parameters
      call a%GetPhysicalParameters(acvis,actco,accph,accvh)
      
   end subroutine
   
   subroutine GatherBase
      implicit none
      integer(ip) :: itime
      
      call e%gather(1_ip,elpre(:,1),a%press(:,1))
      call e%gather(1_ip,elpre(:,2),a%press(:,3))
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      call e%gather(1_ip,eltem(:,1),a%tempe(:,1))
      call e%gather(1_ip,eltem(:,2),a%tempe(:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(1_ip,elpre(:,itime),a%press(:,itime+1)) 
         call e%gather(e%ndime,elvel(:,:,itime),a%veloc(:,:,itime+1)) 
         call e%gather(1_ip,eltem(:,itime),a%tempe(:,itime+1)) 
      enddo

   end subroutine  
   
   subroutine InterpolateBase
      implicit none
      integer(ip) :: itime
      
      call e%interpg(1_ip,elpre(:,1),gppre(1))
      call e%interpg(1_ip,elpre(:,2),gppre(2))
      call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1))
      call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2))
      call e%interpg(1_ip,eltem(:,1),gptem(1))
      call e%interpg(1_ip,eltem(:,2),gptem(2))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%interpg(1_ip,elpre(:,itime),gppre(itime))
         call e%interpg(e%ndime,elvel(:,:,itime),gpvel(:,itime))
         call e%interpg(1_ip,eltem(:,itime),gptem(itime))
      enddo
   end subroutine
   
   subroutine ComputeLinearVariables
      implicit none
      
      gpadp = gppre(1) + a%relpre
      gpadv(:) = gpvel(:,1) 
      gpadt = gptem(1) + a%reltem

   end subroutine

   subroutine ComputeBaseGradients
      
      !Compute vel·grad(V)
      call ComputeAGradV(e,gpadv,AGradV)

      ! Compute element variables gradient 
      call e%gradient(1_ip,elpre(:,1),grpre)
      call e%gradient(e%ndime,elvel(:,:,1),grvel)
      call e%gradient(1_ip,eltem(:,1),grtem)

      ! Compute element velocity divergence
      call e%divergence(elvel(:,:,1),divvel)

      ! Compute vel·grad Var
      call nsc_vgvar(e,elpre(:,1),AGradV,vgpre)
      call nsc_vgvec(e,elvel(:,:,1),AGradV,vgvel)
      call nsc_vgvar(e,eltem(:,1),AGradV,vgtem)

   end subroutine

   subroutine DeallocateBaseElmopeArrays
      !Other arrays alloc
      type(MemoryMan), pointer :: Memor => NULL()
      
      Memor => a%Memor
      
      call Memor%dealloc(e%mnode,a%ncomp-1,elpre,'elpre','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','nsc_pr_elmope')
      call Memor%dealloc(e%mnode,a%ncomp-1,eltem,'eltem','nsc_pr_elmope')
      call Memor%dealloc(a%ncomp-1,gppre,'gppre','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','nsc_pr_elmope')
      call Memor%dealloc(a%ncomp-1,gptem,'gptem','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,grpre,'grpre','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,grtem,'grtem','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,elexv,'elexv','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,eltemv,'eltemv','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,gpadv,'gpadv','nsc_pr_elmope')
      call Memor%dealloc(size(AGradV,1),AGradV,'AGradV','nsc_pr_elmope')
      call Memor%dealloc(e%ndime,vgvel,'vgvel','nsc_pr_elmope')

   end subroutine  

end module         









