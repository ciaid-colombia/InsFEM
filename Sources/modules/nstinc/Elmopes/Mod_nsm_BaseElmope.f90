module Mod_nsm_BaseElmope
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_Element
   use Mod_ConvectiveElement
   use Mod_NavierStokesElement
   use Mod_NavierStokes
   use Mod_ThreeField
   use Mod_nsm_elmdir
   use Mod_NsiExacso
   implicit none
   
   class(NavierStokesProblem), pointer :: a   => NULL()
   class(ThreeFieldNSProblem),   pointer :: sup => NULL()
   
   character(6) :: itask

   !Pointers
      procedure(), pointer :: ProcPointer_ExternalForces => NULL()
      procedure(), pointer :: ProcPointer_ComputeAdvectionVelocity => NULL()
      !Additional pointers (for Elmope only)
      procedure(), pointer :: ProcPointer_PreAssembly => NULL()
      procedure(), pointer :: ProcPointer_PostGaussElmats => NULL()
      procedure(nsm_elmbuv_interf), pointer :: ProcPointer_nsm_elmbuv => NULL()
      procedure(nsm_elmrhu_interf), pointer :: ProcPointer_nsm_elmrhu => NULL()
      procedure(nsm_elmrhp_interf), pointer :: ProcPointer_nsm_elmrhp => NULL()
      procedure(nsm_elmbuq_interf), pointer :: ProcPointer_nsm_elmbuq => NULL()
      procedure(nsm_elmbpv_interf), pointer :: ProcPointer_nsm_elmbpv => NULL()
      
   abstract interface    
      subroutine nsm_elmbuv_interf(e,dvolu,denac,dtinv,vtemp,testf,elmat)
         import
         implicit none
         class(FiniteElement) :: e
         real(rp),    intent(in)    :: vtemp(e%pnode)
         real(rp),    intent(in)    :: testf(e%pnode)
         real(rp),    intent(in)    :: dvolu,denac,dtinv
         real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
      end subroutine
      
      subroutine nsm_elmrhu_interf(e,dvolu,testf,elext,eltemp,elrhs)
         import
         implicit none
         class(FiniteElement) :: e
         real(rp),    intent(in)    :: testf(e%pnode),elext(e%ndime),eltemp(e%ndime)
         real(rp),    intent(in)    :: dvolu
         real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      end subroutine
      
       subroutine nsm_elmrhp_interf(e,timom,dvolu,elext,eltemp,elrhs)
         import
            implicit none
            class(FiniteElement) :: e
            real(rp),    intent(in)    :: elext(e%ndime),eltemp(e%ndime)
            real(rp),    intent(in)    :: timom,dvolu
            real(rp),    intent(inout) :: elrhs(1,e%pnode)

      end subroutine nsm_elmrhp_interf
      
      subroutine nsm_elmbuq_interf(e,timom,dvolu,denac,dtinv,vtemp,elmat)
         import
         implicit none
         class(FiniteElement) :: e
         real(rp),    intent(in)    :: vtemp(e%pnode)
         real(rp),    intent(in)    :: dvolu,timom,denac,dtinv
         real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      end subroutine
      
      subroutine nsm_elmbpv_interf(e,dvolu,testf,elmat)
         import
         implicit none
    
         class(FiniteElement) :: e
         real(rp),    intent(in)    :: testf(e%pnode)
         real(rp),    intent(in)    :: dvolu
         real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      end subroutine
   end interface
      
   !Hooks
      procedure(), pointer :: ProcHook_PreLoop          => NULL()
      procedure(), pointer :: ProcHook_PreAllocate      => NULL()
      procedure(), pointer :: ProcHook_Initializations  => NULL()
      procedure(), pointer :: ProcHook_OnIeltyChange    => NULL()
      procedure(), pointer :: ProcHook_PreGauss         => NULL()
      procedure(), pointer :: ProcHook_ElmatsToZero     => NULL()
      procedure(), pointer :: ProcHook_Gathers          => NULL()
      procedure(), pointer :: ProcHook_InGauss          => NULL()
      procedure(), pointer :: ProcHook_Interpolates     => NULL()
      procedure(), pointer :: ProcHook_PostInterpolate  => NULL()
      procedure(), pointer :: ProcHook_ComputeTaus      => NULL()
      procedure(), pointer :: ProcHook_ComputeTestf     => NULL()
      procedure(), pointer :: ProcHook_InGaussElmats    => NULL()
      procedure(), pointer :: ProcHook_AssemblyEndite   => NULL()
      procedure(), pointer :: ProcHook_Finalizations    => NULL()
      procedure(), pointer :: ProcHook_PostLoop         => NULL()
      procedure(), pointer :: ProcHook_PhysicalProp     => NULL()
      procedure(), pointer :: ProcHook_PreDirichlet     => NULL()
      procedure(), pointer :: ProcHook_Assembly         => NULL()
      procedure(), pointer :: ProcHook_AllocateArrays   => NULL()
      procedure(), pointer :: ProcHook_DeallocateArrays => NULL()
      procedure(), pointer :: ProcHook_InGaussElmatsAssembly => NULL()

   type(NsiExacso) :: exacso  
   class(FiniteElement) , pointer :: e => NULL()

   integer(ip) :: ielem,nelem,igaus,idime,ndime
   real(rp)    :: dvol,dvolt0,dvolt1,dvolt2

   integer(ip) :: ielty0 = 0  !Previous element type
   integer(ip) :: ipnode0 = 0 !Previous pnode number
      
   real(rp), allocatable :: elvel(:,:,:)
   real(rp), allocatable :: elpre(:,:)
   real(rp), allocatable :: elext(:),eltemp(:)
   real(rp), allocatable :: AGradV(:)
   real(rp), allocatable :: gpvel(:,:),gpadv(:)  
   real(rp), allocatable :: grpre(:,:),grvel(:,:)
      
   real(rp)    :: chale(2),timom,tidiv
   real(rp)    :: acden,acvis
   real(rp)    :: reyno
   real(rp)    :: gpvno,divvel,gppre(1)
   
   integer(ip) :: kfl_GoIteInGauss
   
   !Test function
   real(rp), allocatable :: testf(:) 
   
   !Matrices
   real(rp), allocatable :: elmat(:,:,:,:),elrhs(:,:)
   real(rp), allocatable :: elmuv(:,:,:,:),elmpq(:,:,:,:),elmpv(:,:,:,:),elmuq(:,:,:,:)
   real(rp), allocatable :: elrhu(:,:),elrhp(:,:)
   real(rp), allocatable :: wrmat1(:,:)
   
   !Reference Dtinv
   !This one is for things which have their own integration scheme (taus, dynamic subscales)
   !It is necessary only for CrankNicolson
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: nsteps
   real(rp)                :: LHSdtinv
   real(rp)                :: ReferenceDtinv
   
   !To define multy materials
   integer(ip) :: imat=1_ip
   
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
      ProcPointer_ExternalForces => NULLSUB
      ProcPointer_ComputeAdvectionVelocity => NULLSUB
   
   !Hooks
      ProcHook_PreLoop          => NULLSUB
      ProcHook_PreAllocate      => NULLSUB
      ProcHook_Initializations  => NULLSUB
      ProcHook_OnIeltyChange    => NULLSUB
      ProcHook_PreGauss         => NULLSUB
      ProcHook_ElmatsToZero     => NULLSUB
      ProcHook_Gathers          => NULLSUB
      ProcHook_InGauss          => NULLSUB
      ProcHook_Interpolates     => NULLSUB
      ProcHook_PostInterpolate  => NULLSUB
      ProcHook_ComputeTaus      => NULLSUB
      ProcHook_ComputeTestf     => NULLSUB
      ProcHook_InGaussElmats    => NULLSUB
      ProcHook_AssemblyEndite   => NULLSUB
      ProcHook_Finalizations    => NULLSUB
      ProcHook_PostLoop         => NULLSUB
      ProcHook_PhysicalProp     => NULLSUB
      ProcHook_PreDirichlet     => NULLSUB
      ProcHook_AllocateArrays   => NULLSUB
      ProcHook_DeallocateArrays => NULLSUB
      ProcHook_Assembly         => LinearSystemAssembly
      ProcHook_InGaussElmatsAssembly => NULLSUB
   end subroutine
   
   subroutine AllocateBaseElmopeArrays
      !Set Time Integrator
      call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
      ReferenceDtinv = a%dtinv
      !If CrankNicolson the reference is 1/2 dt, used in dynamic subscales, taus etc
      if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = 2*a%dtinv
      
      !Other arrays alloc
      call a%Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','nsm_EnditeElmope')
      call a%Memor%alloc(        e%mnode,a%ncomp-1,elpre,'elpre','nsm_EnditeElmope')
      call a%Memor%alloc(e%ndime+1,elext,'elext','nsm_EnditeElmope')
      call a%Memor%alloc(e%ndime,eltemp,'eltemp','nsm_EnditeElmope')
      call a%Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','nsm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpadv,'gpadv','nsm_EnditeElmope')
      call a%Memor%alloc(e%mnode,AGradV,'AGradV','nsm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','nsm_EnditeElmope')
      call a%Memor%alloc(1_ip     ,e%ndime,grpre,'grpre','nsm_EnditeElmope')
      
      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)
      
   end subroutine
   
   subroutine VelocityAndPressureGathers
      implicit none
      integer(ip) :: itime
      
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(e%ndime,elvel(:,:,itime),a%veloc(:,:,itime+1)) 
      enddo
      call e%gather(1_ip   ,elpre(1:e%pnode,1),a%press(1:e%pnode,1))
   end subroutine  
   
   subroutine InterpolateGpVelocities
      implicit none
      integer(ip) :: itime
      
      call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1))
      call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%interpg(e%ndime,elvel(:,:,itime),gpvel(:,itime))
      enddo
   end subroutine
   
   subroutine DeallocateBaseElmopeArrays
      !Other arrays alloc
      call a%Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','nsm_EnditeElmope')
      call a%Memor%dealloc(        e%mnode,a%ncomp-1,elpre,'elpre','nsm_EnditeElmope')   
      call a%Memor%dealloc(e%ndime+1,elext,'elext','nsm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,eltemp,'eltemp','nsm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','nsm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpadv,'gpadv','nsm_EnditeElmope')
      call a%Memor%dealloc(size(AGradV,1),AGradV,'AGradV','nsm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','nsm_EnditeElmope')
      call a%Memor%dealloc(1_ip   ,e%ndime,grpre,'grpre','nsm_EnditeElmope') 
   end subroutine
   
   subroutine LinearSystemAssembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   end subroutine

end module
