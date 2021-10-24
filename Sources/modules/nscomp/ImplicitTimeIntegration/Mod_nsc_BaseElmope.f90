module Mod_nsc_BaseElmope
   use Mod_nsc_elmdir_im
   use Mod_Mesh
   use Mod_Memor
   use Mod_Element
   use Mod_php_SetTimeIntegrator
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_NSCompressibleElement
   use Mod_NSCompressibleImplicit
   use Mod_NSCompressibleImplicitElement
   use Mod_NSCompressibleSubroutines
   use Mod_NscExacso   
   use Mod_NscDamping   
   implicit none
   
   class(NSCompressibleImplicitProblem), pointer :: a => NULL()
   
   character(6) :: itask

   !Pointers
      procedure(), pointer :: ProcPointer_nsc_ExternalForces => NULL()
      procedure(), pointer :: ProcPointer_nsc_ComputeVariables => NULL()
      procedure(), pointer :: ProcPointer_nsc_ComputeConvectionCoefficients => NULL()
      procedure(), pointer :: ProcPointer_nsc_ComputeDiffusionCoefficients => NULL()
      procedure(), pointer :: ProcPointer_nsc_ComputeTransportCoefficients => NULL()
      !Additional pointers (for Elmope only)
      !the elemental matrix takes the form
      !************************************************
      !q·  elmbdq  elmbmq  elmbeq  rho =  elmrhd
      !n·  elmbdn  elmbmn  elmben  mom =  elmrhm
      !g·  elmbdg  elmbmg  elmbeg  ene =  elmrhe
      !************************************************
         procedure(), pointer :: ProcPointer_nsc_PreAssembly => NULL()
         procedure(), pointer :: ProcPointer_nsc_PostGaussElmats => NULL()
         procedure(nsc_elmbdq), pointer :: ProcPointer_nsc_elmbdq => NULL()
         procedure(nsc_elmbdn), pointer :: ProcPointer_nsc_elmbdn => NULL()
         procedure(nsc_elmbdg), pointer :: ProcPointer_nsc_elmbdg => NULL()
         procedure(nsc_elmbmq), pointer :: ProcPointer_nsc_elmbmq => NULL()
         procedure(nsc_elmbmn), pointer :: ProcPointer_nsc_elmbmn => NULL()
         procedure(nsc_elmbmg), pointer :: ProcPointer_nsc_elmbmg => NULL()
         procedure(nsc_elmbeq), pointer :: ProcPointer_nsc_elmbeq => NULL()
         procedure(nsc_elmben), pointer :: ProcPointer_nsc_elmben => NULL()
         procedure(nsc_elmbeg), pointer :: ProcPointer_nsc_elmbeg => NULL()
         procedure(nsc_elmrhd), pointer :: ProcPointer_nsc_elmrhd => NULL()
         procedure(nsc_elmrhm), pointer :: ProcPointer_nsc_elmrhm => NULL()
         procedure(nsc_elmrhe), pointer :: ProcPointer_nsc_elmrhe => NULL()
    
      
   
   !Hooks
      procedure(), pointer :: ProcHook_nsc_PreLoop => NULL()
      procedure(), pointer :: ProcHook_nsc_PreAllocate => NULL()
      procedure(), pointer :: ProcHook_nsc_Initializations => NULL()
      procedure(), pointer :: ProcHook_nsc_OnIeltyChange => NULL()
      procedure(), pointer :: ProcHook_nsc_PreGauss => NULL()
      procedure(), pointer :: ProcHook_nsc_ElmatsToZero => NULL()
      procedure(), pointer :: ProcHook_nsc_Gathers => NULL()
      procedure(), pointer :: ProcHook_nsc_InGauss => NULL()
      procedure(), pointer :: ProcHook_nsc_Interpolates => NULL()
      procedure(), pointer :: ProcHook_nsc_PostInterpolate => NULL()
      procedure(), pointer :: ProcHook_nsc_ComputeTaus => NULL()
      procedure(), pointer :: ProcHook_nsc_ComputeTestf => NULL()
      procedure(), pointer :: ProcHook_nsc_InGaussElmats => NULL()
      procedure(), pointer :: ProcHook_nsc_InGaussElmatsAssembly => NULL()
      procedure(), pointer :: ProcHook_nsc_AssemblyEndite => NULL()
      procedure(), pointer :: ProcHook_nsc_Finalizations => NULL()
      procedure(), pointer :: ProcHook_nsc_PostLoop => NULL()
      procedure(), pointer :: ProcHook_nsc_PhysicalProp     => NULL()
      procedure(), pointer :: ProcHook_nsc_PreDirichlet => NULL()
                                       
                                       
   type(NscExacso) :: exacso           
                                       
   type(NscDamping) :: damping           

   class(FiniteElement), pointer :: e => NULL()       
   integer(ip) :: ielem,nelem
  
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv

   integer(ip)           :: ielty0 = 0  !Previous element type
   integer(ip)           :: ipnode0 = 0 !Previous pnode number
      
   real(rp), allocatable :: elden(:,:)
   real(rp), allocatable :: elmom(:,:,:)
   real(rp), allocatable :: elene(:,:)
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: eltemm(:)
   real(rp), allocatable :: elexm(:)
   real(rp), allocatable :: AGradV(:)
   real(rp), allocatable :: gpden(:)  
   real(rp), allocatable :: gpmom(:,:)
   real(rp), allocatable :: gpene(:)  
   real(rp), allocatable :: gpadv(:),gpadm(:)  
   real(rp), allocatable :: vgmom(:)
   real(rp), allocatable :: grden(:), grmom(:,:), grene(:)
   real(rp), allocatable :: hess(:,:,:), lapl(:) 
      
   real(rp)    :: chale(2),timom(3), tidiv
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: gpadd,gpade
   real(rp)    :: gpvno,divvel
   real(rp)    :: gpmno,divmom
   real(rp)    :: vgden,vgene
   real(rp)    :: invcvh,invgpd,sqinvgpd,sqgpmn,sqgpvn
   real(rp)    :: acgamma,actcn,aux,aux_t,aux_d
   real(rp)    :: gppre,gptem,gpspd
   real(rp)    :: eltemd(1),elteme(1)
   real(rp)    :: elexd(1),elexe(1)
   
   integer(ip) :: igaus
   integer(ip) :: kfl_GoIteInGauss
   
   !Convection Matrix coefficients
   real(rp), allocatable :: Add(:),   Adm(:,:)
   real(rp), allocatable :: Amd(:,:), Amm(:,:,:), Ame(:,:)
   real(rp), allocatable :: Aed(:),   Aem(:,:),   Aee(:)

   !Diffusion Matrix coefficients
   real(rp), allocatable :: KGmd(:,:,:), KGmm(:,:,:,:)
   real(rp), allocatable :: KGed(:,:),   KGem(:,:,:),   KGee(:,:)
   real(rp), allocatable :: Kmd(:,:), Kmm(:,:,:)
   real(rp), allocatable :: Ked(:),   Kem(:,:),   Kee(:)

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
   real(rp), allocatable :: JAmd(:,:), JAmm(:,:,:)
   real(rp), allocatable :: JAed(:),   JAem(:,:),   JAee(:)
   
   !Matrices
   real(rp), allocatable :: elmdq(:,:), elmmq(:,:,:), elmeq(:,:)
   real(rp), allocatable :: elmdn(:,:,:), elmmn(:,:,:,:), elmen(:,:,:)
   real(rp), allocatable :: elmdg(:,:), elmmg(:,:,:), elmeg(:,:)
   real(rp), allocatable :: elrhd(:),     elrhm(:,:),     elrhe(:)
   
   real(rp) :: dvolt0, dvolt1, dvolt2
   
   !Reference Dtinv
   !This one is for things which have their own integration scheme (taus, dynamic subscales)
   !It is necessary only for CrankNicolson
   real(rp) :: ReferenceDtinv
   
   !#$COMPOSEPROCS 100
#include "COMPOSEPROCS_POINTERS_100.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_100.i90"

  
  
   
   
   
   subroutine SetPointersAndHooksToNULLSUB
   implicit none
      !External Procedures
      procedure() :: NULLSUB
   
   !Pointers
      ProcPointer_nsc_ExternalForces => NULLSUB
      ProcPointer_nsc_ComputeVariables => NULLSUB
      ProcPointer_nsc_ComputeConvectionCoefficients => NULLSUB
      ProcPointer_nsc_ComputeDiffusionCoefficients => NULLSUB
      ProcPointer_nsc_ComputeTransportCoefficients => NULLSUB
   
   !Hooks
      ProcHook_nsc_PreLoop => NULLSUB
      ProcHook_nsc_PreAllocate => NULLSUB
      ProcHook_nsc_Initializations => NULLSUB
      ProcHook_nsc_OnIeltyChange => NULLSUB
      ProcHook_nsc_PreGauss => NULLSUB
      ProcHook_nsc_ElmatsToZero => NULLSUB
      ProcHook_nsc_Gathers => NULLSUB
      ProcHook_nsc_InGauss => NULLSUB
      ProcHook_nsc_Interpolates => NULLSUB
      ProcHook_nsc_PostInterpolate => NULLSUB
      ProcHook_nsc_ComputeTaus => NULLSUB
      ProcHook_nsc_ComputeTestf => NULLSUB
      ProcHook_nsc_InGaussElmats => NULLSUB
      ProcHook_nsc_InGaussElmatsAssembly => NULLSUB
      ProcHook_nsc_AssemblyEndite => NULLSUB
      ProcHook_nsc_Finalizations => NULLSUB
      ProcHook_nsc_PostLoop => NULLSUB
      ProcHook_nsc_PhysicalProp => NULLSUB
      ProcHook_nsc_PreDirichlet => NULLSUB
   end subroutine
   
   subroutine AllocateBaseElmopeArrays
   
      type(MemoryMan), pointer :: Memor => NULL()
      !Set Time Integrator
      call php_SetTimeIntegrator(a,Integrator,LHSdtinv,nsteps)
      ReferenceDtinv = a%dtinv
      !If CrankNicolson the reference is 1/2 dt, used in dynamic subscales, taus etc
      if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = 2*a%dtinv
      
      Memor => a%Memor
      !Other arrays alloc
      call Memor%alloc(e%ndime,e%mnode,elvel,'elvel','nsc_elmope_im')
      call Memor%alloc(e%mnode,a%ncomp-1,elden,'elden','nsc_elmope_im')
      call Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elmom,'elmom','nsc_elmope_im')
      call Memor%alloc(e%mnode,a%ncomp-1,elene,'elene','nsc_elmope_im')
      call Memor%alloc(a%ncomp-1,gpden,'gpden','nsc_elmope_im')
      call Memor%alloc(e%ndime,a%ncomp-1,gpmom,'gpmom','nsc_elmope_im')
      call Memor%alloc(a%ncomp-1,gpene,'gpene','nsc_elmope_im')
      call Memor%alloc(e%ndime,grden,'grden','nsc_elmope_im')
      call Memor%alloc(e%ndime,e%ndime,grmom,'grmom','nsc_elmope_im')
      call Memor%alloc(e%ndime,grene,'grene','nsc_elmope_im')
      call Memor%alloc(e%ndime,elexm,'elexm','nsc_elmope_im')
      call Memor%alloc(e%ndime,eltemm,'eltemm','nsc_elmope_im')
      call Memor%alloc(e%ndime,gpadv,'gpadv','nsc_elmope_im')
      call Memor%alloc(e%ndime,gpadm,'gpadm','nsc_elmope_im')
      call Memor%alloc(e%mnode,AGradV,'AGradV','nsc_elmope_im')
      call Memor%alloc(e%ndime,vgmom,'vgmom','nsc_elmope_im')

      !Physical Parameters
      call a%GetPhysicalParameters(acvis,actco,accph,accvh)
      
   end subroutine
   
   subroutine GatherBase
      implicit none
      integer(ip) :: itime
      
      call e%gather(1_ip,elden(:,1),a%densf(:,1))
      call e%gather(1_ip,elden(:,2),a%densf(:,3))
      call e%gather(e%ndime,elmom(:,:,1),a%momen(:,:,1))
      call e%gather(e%ndime,elmom(:,:,2),a%momen(:,:,3))
      call e%gather(1_ip,elene(:,1),a%energ(:,1))
      call e%gather(1_ip,elene(:,2),a%energ(:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(1_ip,elden(:,itime),a%densf(:,itime+1)) 
         call e%gather(e%ndime,elmom(:,:,itime),a%momen(:,:,itime+1)) 
         call e%gather(1_ip,elene(:,itime),a%energ(:,itime+1)) 
      enddo

   end subroutine  
   
   subroutine CalculateElmchl
      
      call nsc_ComputeElementVelocity(e,elden(:,1),elmom(:,:,1),elvel)
      call elmchl(e,1_ip,elvel,chale)

   end subroutine

   subroutine InterpolateBase
      implicit none
      integer(ip) :: itime
      
      call e%interpg(1_ip,elden(:,1),gpden(1))
      call e%interpg(1_ip,elden(:,2),gpden(2))
      call e%interpg(e%ndime,elmom(:,:,1),gpmom(:,1))
      call e%interpg(e%ndime,elmom(:,:,2),gpmom(:,2))
      call e%interpg(1_ip,elene(:,1),gpene(1))
      call e%interpg(1_ip,elene(:,2),gpene(2))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%interpg(1_ip,elden(:,itime),gpden(itime))
         call e%interpg(e%ndime,elmom(:,:,itime),gpmom(:,itime))
         call e%interpg(1_ip,elene(:,itime),gpene(itime))
      enddo
   end subroutine
   
   subroutine ComputeLinearVariables
      implicit none
      
      gpadd = gpden(1)
      gpadm(:) = gpmom(:,1)
      gpade = gpene(1)

   end subroutine

   subroutine ComputeBaseGradients
      
      !Compute vel·grad(V)
      call ComputeAGradV(e,gpadv,AGradV)

      ! Compute element variables gradient 
      call e%gradient(1_ip,elden(:,1),grden)
      call e%gradient(e%ndime,elmom(:,:,1),grmom)
      call e%gradient(1_ip,elene(:,1),grene)

      ! Compute element momentum divergence
      call e%divergence(elmom(:,:,1),divmom)

      ! Compute (mom/rho)·grad Var
      call nsc_vgvar(e,elden(:,1),AGradV,vgden)
      call nsc_vgvec(e,elmom(:,:,1),AGradV,vgmom)
      call nsc_vgvar(e,elene(:,1),AGradV,vgene)

   end subroutine

   subroutine DeallocateBaseElmopeArrays
      !Other arrays alloc
      type(MemoryMan), pointer :: Memor => NULL()

      Memor => a%Memor      
      call Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','nsc_elmope_im')
      call Memor%dealloc(e%mnode,a%ncomp-1,elden,'elden','nsc_elmope_im')
      call Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elmom,'elmom','nsc_elmope_im')
      call Memor%dealloc(e%mnode,a%ncomp-1,elene,'elene','nsc_elmope_im')
      call Memor%dealloc(a%ncomp-1,gpden,'gpden','nsc_elmope_im')
      call Memor%dealloc(e%ndime,a%ncomp-1,gpmom,'gpmom','nsc_elmope_im')
      call Memor%dealloc(a%ncomp-1,gpene,'gpene','nsc_elmope_im')
      call Memor%dealloc(e%ndime,grden,'grden','nsc_elmope_im')
      call Memor%dealloc(e%ndime,e%ndime,grmom,'grmom','nsc_elmope_im')
      call Memor%dealloc(e%ndime,grene,'grene','nsc_elmope_im')
      call Memor%dealloc(e%ndime,elexm,'elexm','nsc_elmope_im')
      call Memor%dealloc(e%ndime,eltemm,'eltemm','nsc_elmope_im')
      call Memor%dealloc(e%ndime,gpadv,'gpadv','nsc_elmope_im')
      call Memor%dealloc(e%ndime,gpadm,'gpadm','nsc_elmope_im')
      call Memor%dealloc(size(AGradV,1),AGradV,'AGradV','nsc_elmope_im')
      call Memor%dealloc(e%ndime,vgmom,'vgmom','nsc_elmope_im')
   end subroutine  
end module         









