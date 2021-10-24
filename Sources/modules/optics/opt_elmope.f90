module mod_opt_elmope
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_Element
   use Mod_Optics
   use Mod_ConvectiveElement
   use Mod_nsm_Viscosity
   implicit none
   class(OpticsProblem), pointer :: a => NULL()
   
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ielem,nelem,ipoin,npoin
   integer(ip) :: igaus
   integer(ip) :: ielty0 = 0

   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: PreGaussElmats => NULL()

   end type
   type(PPointer) :: ProcPointer
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      
      
   end type
   type(PHook) :: ProcHook

   real(rp), allocatable :: elvel(:,:),elpre(:)
   real(rp), allocatable :: eltem(:)
   real(rp), allocatable :: gpvel(:)
   real(rp), allocatable :: grvel(:,:), grtem(:)
   real(rp) :: gptem(1),gppre(1)
   real(rp), allocatable :: elct2(:), elcn2(:), elcn2u53(:)
   
   real(rp), allocatable :: eldisNS(:),eldisTempe(:)
   
   real(rp) :: TGradNorm, VGradNorm, TGradNorm2, VGradNorm2
   real(rp) :: Tdissip(1), Vdissip(1)
   
   real(rp) :: acden,acvis,acsph,actco
   real(rp) :: vista
   
   real(rp) :: dvol,dvolt0
   
   real(rp) :: gpcn2,gpct2,vnor53
   
   real(rp) :: chale(2)
   
   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"
   
   !----------------------------------------------------------
   !SetPointers
   subroutine SetPointers
   use typre
      implicit none

      !External Procedures
      procedure() :: NULLSUB

      integer(ip) :: kfl_nonlinear,nelty
      
      call ResetProcedureComposition
      
      !-----------------------------------------------------------
      !Defaults
      
      !Pointers
      ProcPointer%PreGaussElmats  => NULLSUB
      
      
      !Hooks
      ProcHook%OnIeltyChange   => NULLSUB
      ProcHook%InGauss         => NULLSUB
      
      !-----------------------------------------------------------
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
      endif
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
   end subroutine
   
   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussNonLinear
      implicit none
      
      call e%elmder
   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module


subroutine opt_elmope(OPTProblem)
   use typre
   use Mod_Element
   use Mod_Optics
   use mod_opt_elmope
   use Mod_nsm_Viscosity
   implicit none
   class(OpticsProblem), target :: OPTProblem
   
   integer(ip) :: inode
   
   a=>OPTProblem
   
   !Initialize cn2 and ct2 arrays
   a%cn2(:) = 0.0_rp
   a%cn2u53(:) = 0.0_rp
   a%ct2 = 0.0_rp
   
   !If Dissip is computed here through Smagorinsky, the NSDissipation and TempeDissipation arrays are computed here
   if (a%kfl_dissi /= 1) then
      a%NSDissipation = 0.0_rp
      a%TempeDissipation = 0.0_rp
   endif
   
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Allocations
   !We force a closed rule for smoothing, (we are using vmass)
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'opt_elmope')
   
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','opt_elmope')
   call a%Memor%alloc(e%mnode,elpre,'elpre','opt_elmope')
   call a%Memor%alloc(e%mnode,eltem,'eltem','opt_elmope')
   call a%Memor%alloc(e%ndime,gpvel,'gpvel','opt_gpvel')
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','opt_grvel')
   call a%Memor%alloc(e%ndime,grtem,'grtem','opt_grtem')
   call a%Memor%alloc(e%mnode,elct2,'elct2','opt_grtem')
   call a%Memor%alloc(e%mnode,elcn2,'elcn2','opt_grtem')
   call a%Memor%alloc(e%mnode,elcn2u53,'elcn2u53','opt_grtem')
   call a%Memor%alloc(e%mnode,eldisNS,'eldisNS','opt_elmope')
   call a%Memor%alloc(e%mnode,eldisTempe,'eldisTempe','opt_elmope')

   
   
   !Physical Parameters
   call a%GetPhysicalParameters(acden,acvis,acsph,actco)
   
   !Hook
   !call ProcHook%Initializations  
   
   call a%Mesh%GetNelem(nelem)
   elements: do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !ElmatsToZero
      elct2 = 0.0_rp
      elcn2 = 0.0_rp
      elcn2u53 = 0.0_rp
      
      !Gathers
      call e%gather(e%ndime,elvel,a%veloc)
      call e%gather(1,elpre,a%press)
      call e%gather(1,eltem,a%tempe)
      
      if (a%kfl_dissi /= 1) then
         eldisNS = 0.0_rp
         eldisTempe = 0.0_rp
      
      elseif (a%kfl_dissi == 1) then
         call e%gather(1,eldisNS,a%NSDissipation)
         call e%gather(1,eldisTempe,a%TempeDissipation)
      endif
   
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen
      
      !Compute the characteristic length chale
      call elmchl(e,1_ip,elvel,chale)

      !DvolsToZero
      dvolt0=0.0_rp
      
      call ProcPointer%PreGaussElmats
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Physical Parameters
         call a%GetPhysicalParameters(acden,acvis,acsph,actco)
         
         !Hook
         call ProcHook%InGauss
                 
         !Update dvols
         dvol = e%weigp(e%igaus)*e%detjm
         dvolt0 = dvolt0 + dvol
         
         !Interpolate
         call e%interpg(e%ndime,elvel(:,:),gpvel(:))
         call e%interpg(1,elpre(:),gppre)
         call e%interpg(1,eltem(:),gptem)
         call e%gradient(e%ndime,elvel,grvel)
         call e%gradient(1,eltem,grtem)
         
         !Smagorinsky or WALE (both LES)
         if (a%kfl_dissi == 0 .or. a%kfl_dissi == 2) then
            !Velocity Symmetric gradient norm square
            call vecnor(grvel,e%ndime*e%ndime,VGradNorm,2)
            VGradNorm2 = VGradNorm**2
               
            !Temperature gradient norm square
            call vecnor(grtem,e%ndime,TGradNorm,2)
            TGradNorm2 = TGradNorm**2
            
            !Smagorinsky
            if (a%kfl_dissi == 0) then
               call nsm_smago(e,grvel,acden,a%turbu(1),vista)
            
            !Wale
            elseif (a%kfl_dissi == 2) then
               call nsm_wale(e,grvel,acden,a%turbu(1),vista)
            endif
            
            acvis = acvis + vista
               
            Vdissip = acvis*VGradNorm2
            Tdissip = (actco/acsph+vista/a%Prtur)*TGradNorm2
            
            do inode = 1,e%pnode
               eldisNS(inode) = eldisNS(inode) + e%shape(inode,e%igaus)*Vdissip(1)*dvol
               eldisTempe(inode) = eldisTempe(inode) + e%shape(inode,e%igaus)*Tdissip(1)*dvol
            enddo
           
         !Subscales coming from NS and TEMPE
         elseif (a%kfl_dissi == 1) then
             call e%interpg(1,eldisNS,Vdissip)
             call e%interpg(1,eldisTempe,Tdissip)
         endif
         
         if (abs(Vdissip(1)) < 1e-14) then
            gpct2 = 0.0_rp
         else   
            gpct2 = (a%a_Obukhov**2)*Tdissip(1)*(abs(Vdissip(1))**(-1.0_rp/3.0_rp))
         endif
         
         if (a%units == 0) then
            !Pressure from Pascals to milibars + 1 atmosphere (1000 milibars)
            gppre = (gppre+1e5)/100.0_rp
            !Temperature to Kelvin
            gptem = gptem + 273.15
         endif   
         !Be careful with unit systems
         gpcn2 = (79e-6*gppre(1)/gptem(1)**2)**2*gpct2
         
         !Compute velocity ^5/3
         call vecnor(gpvel,e%ndime,vnor53,2)
         vnor53 = vnor53**(5.0_rp/3.0_rp)
         
         
         !Elemental assembly ct2
         do inode = 1,e%pnode
            elct2(inode) = elct2(inode) + e%shape(inode,e%igaus)*gpct2*dvol
            elcn2(inode) = elcn2(inode) + e%shape(inode,e%igaus)*gpcn2*dvol
            elcn2u53(inode) = elcn2u53(inode) + e%shape(inode,e%igaus)*gpcn2*vnor53*dvol
         enddo
         
      
      enddo gauss_points
   
      !Global Assembly
      a%ct2(e%lnods(1:e%pnode)) = a%ct2(e%lnods(1:e%pnode)) + elct2(1:e%pnode)
      a%cn2(e%lnods(1:e%pnode)) = a%cn2(e%lnods(1:e%pnode)) + elcn2(1:e%pnode)
      a%cn2u53(e%lnods(1:e%pnode)) = a%cn2u53(e%lnods(1:e%pnode)) + elcn2u53(1:e%pnode)
      
      if (a%kfl_dissi /= 1) then
         a%NSDissipation(e%lnods(1:e%pnode)) = a%NSDissipation(e%lnods(1:e%pnode)) + eldisNS(1:e%pnode)
         a%TempeDissipation(e%lnods(1:e%pnode)) = a%TempeDissipation(e%lnods(1:e%pnode)) + eldisTempe(1:e%pnode)
         
      endif
   
   enddo elements
   
   !Smoothings 
   call a%Project(1,a%ct2)
   call a%Project(1,a%cn2(:))
   call a%Project(1,a%cn2u53(:))
   
   !Smagorinsky, computed in Optics module
   if (a%kfl_dissi /= 1) then
      call a%Project(1,a%NSDissipation)
      call a%Project(1,a%TempeDissipation)
   endif
   
   !--------------------------------------------------------
   !Averages
   a%avg_press = a%avg_press + a%press
   call a%Mesh%GetNpoin(npoin)
   forall (ipoin = 1:npoin)
      a%avg_vnorm(ipoin) = a%avg_vnorm(ipoin) + sqrt(dot_product(a%veloc(:,ipoin),a%veloc(:,ipoin)))
   endforall
   a%avg_tempe = a%avg_tempe + a%tempe
   
   
   !Subscales
   a%avg_vdiss = a%avg_vdiss + a%NSDissipation
   a%avg_tdiss = a%avg_tdiss + a%TempeDissipation
   
      
   a%avg_cn2(:,1) = a%avg_cn2(:,1) + a%cn2(:)
   a%avg_cn2u53(:,1) = a%avg_cn2u53(:,1) + a%cn2u53(:)
   
   
   !---------------------------------------------------------
   !Deallocates 
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','opt_elmope')
   call a%Memor%dealloc(e%mnode,elpre,'elpre','opt_elmope')
   call a%Memor%dealloc(e%mnode,eltem,'eltem','opt_elmope')
   call a%Memor%dealloc(e%ndime,gpvel,'gpvel','opt_gpvel')
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','opt_grvel')
   call a%Memor%dealloc(e%ndime,grtem,'grtem','opt_grtem')
   call a%Memor%dealloc(e%mnode,elct2,'elct2','opt_grtem')
   call a%Memor%dealloc(e%mnode,elcn2,'elcn2','opt_grtem')
   call a%Memor%dealloc(e%mnode,elcn2u53,'elcn2u53','opt_grtem')
   call a%Memor%dealloc(e%mnode,eldisNS,'eldisNS','opt_elmope')
   call a%Memor%dealloc(e%mnode,eldisTempe,'eldisTempe','opt_elmope')
   
   call a%Mesh%ElementDealloc(e,a%Memor,a%EndLoopQuadrature,'opt_elmope')
      
end subroutine
