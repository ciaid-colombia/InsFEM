module Mod_supf_elmope_hydro
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_Mesh
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_NavierStokesElement
   use Mod_SUPF_Element   
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_sup_elmdir
   use Mod_nsm_elmdir
   use Mod_php_Elmdir   
   use Mod_ThreeFieldElement
   use Mod_SupExacso    
   use Mod_SupOperations

   implicit none   
   
   class(SUPFractionalStepProblem), pointer :: a
   type(SupExacso) :: exacso     
   real(rp), parameter :: zensi = 0.0_rp   
   
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(), NOPASS, pointer :: FreeSurface => NULL()
      procedure(), NOPASS, pointer :: PreAssembly => NULL()      

   end type
   type(PPointer) :: ProcPointer
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: Gathers => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: Interpolates => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()     
      
    end type
   type(PHook) :: ProcHook 
   
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ielem,nelem
  
   real(rp), allocatable :: elvel(:,:,:)  
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:),elrhu(:,:),elrhp(:,:)
   real(rp)              :: dvol
   
   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv   


   integer(ip)           :: ielty0 = 0 !Previous element type
      

   real(rp), allocatable :: elext(:)
   real(rp), allocatable :: elmuv(:,:,:,:),elmuq(:,:,:,:),elmpq(:,:,:,:),elmpv(:,:,:,:)
   real(rp), allocatable :: wrmat1(:,:)
   real(rp), allocatable :: AGradV(:), testf(:)  
   real(rp), allocatable :: gpadv(:)   

      
   real(rp)    :: chale(2),timom,dvolt0,dvolt1
   integer(ip) :: nmean,auxtens
   real(rp)    :: acden,acvis
   real(rp)    :: reyno
   real(rp)    :: gpvno,gpvno2   

   
   integer(ip) :: idime,itime,igaus,itest
   integer(ip) :: currentbvess,bcstar,auxpba1,auxpba2
   
   !Level Set
   integer(ip)              :: inode,ipoin
   integer(ip)              :: elemStatus,ngauss_minus,ngauss_plus,ngaus_total,poinStatus
   real(rp), allocatable    :: weigp(:),xloc(:,:)   
   
   !todo multy materials
   integer(ip) :: imat=1

   
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
      ProcPointer%PostGaussElmats => PostGaussElmats
      ProcPointer%FreeSurface => NULLSUB
      !PreAssembly to use in free-surface and in enriched elements
      ProcPointer%PreAssembly => NULLSUB      
      !Hooks
      ProcHook%Initializations => NULLSUB
      ProcHook%Gathers         => NULLSUB
      ProcHook%OnIeltyChange   => NULLSUB
      ProcHook%PreGauss        => NULLSUB
      ProcHook%Interpolates    => NULLSUB
      ProcHook%InGauss         => NULLSUB
      ProcHook%InGaussElmats   => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      ProcHook%PhysicalProp    => NULLSUB
      
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      
      !--------------------------------------------------------------------------------
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
      
      if(kfl_nonlinear==1 .or. a%kfl_colev==1)then
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
         call ConcatenateProcedures(ProcHook%InGaussElmats,ProcPointer%PostGaussElmats)
         ProcPointer%PostGaussElmats => NULLSUB      
      end if      
      
      !--------------------------------------------------------------------
      !Level set
      if(a%kfl_colev==1)then
         call ConcatenateProcedures(ProcHook%Initializations,AllocLevelsetTF) 
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocLevelSetTF)
         call PrependProcedure(ProcHook%PreGauss,CutelementsTF)
         call PrependProcedure(ProcHook%PhysicalProp,FluidPropertiesTF) 
         if(a%kfl_fsurf==1)then
            !Always we disconnecting domains
            if(a%kfl_fsurfLapla==1)then
               call PrependProcedure(ProcPointer%PreAssembly,ElmatsToLapla)            
            end if
            call ConcatenateProcedures(ProcPointer%PreAssembly,FreeSurfMatsTF)  
         end if
      end if     

      
      
      
   end subroutine  
   

  !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussNonLinear
      implicit none
      
      call e%elmder
      call e%elmhes
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine  
   
   !----------------------------------------------------------
   !PostGauss Matrices
   subroutine PostGaussElmats
      implicit none  
 
   end subroutine
   
   !-------------------------------------------------------
   !LevelSet
   
   subroutine AllocLevelsetTF
      implicit none   
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_1st')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,xloc,'xloc','supf_elmope_1st') 
   
   end subroutine
   
   subroutine DeallocLevelSetTF
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_1st')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,xloc,'xloc','supf_elmope_1st') 
      
   
   end subroutine
   
   subroutine CutelementsTF
      implicit none
      integer(ip) :: elemStatus
      
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus ==0)then    
         
         ngaus_total=0
     
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         call a%CutMesh%GetWeigpCut(ielem,e%ndime,weigp)
         call a%CutMesh%GetXlocCut(ielem,e%ndime,xloc)
         
         ngaus_total = ngauss_minus+ngauss_plus
         
         !The rutine give the needed shape functions associated 
         if(a%kfl_fsurf==1) weigp(1:ngauss_minus)=0.0_rp         
         
         call e%SetParticularGaussPoints(a%Memor,ngaus_total,xloc,weigp(:))

      end if
   
   end subroutine
   
   subroutine FluidPropertiesTF
      implicit none
      integer(ip) :: elemStatus
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==1)then
         imat=1
      elseif(elemStatus==-1)then
         imat=2
      elseif(elemStatus ==0)then
         
         ngauss_minus=0
         ngauss_plus=0
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         
         if(e%igaus<=ngauss_minus)then
            imat=2
         elseif(e%igaus>ngauss_minus)then  
            imat=1
         end if
         
      end if  

   
   end subroutine 
   
   subroutine FreeSurfMatsTF
      implicit none 
      integer(ip)  :: elemStatus,poinStatus       
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==-1)then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%CutMesh%GetPointType(ipoin,poinStatus)            
            if(poinStatus==-1)then
            
               elmat(:,inode,:,1:e%pnode) = 0.0_rp
               elrhs(:,inode)=0.0_rp              
            
            end if         
         end do        
      end if
   end subroutine    
   
   subroutine ElmatsToLapla
      implicit none
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin,idime,ntens 
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      ntens=(e%ndime-1)*(e%ndime-1)+2
      
      
      if(elemStatus==-1)then
         
         elmat(:,:,:,:) = 0.0_rp        
         elrhs(:,:)=0.0_rp
        
         
         do inode=1,e%pnode
            elmat((e%ndime+1),inode,(e%ndime+1),inode) = 1.0_rp
         end do 
         
         forall (idime = 1:(e%ndime))
            elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
         end forall           

         
      end if      
      
   end subroutine    

   

   
  
end module   

subroutine supf_elmope_hydro(SUPFProblem)
   use Mod_SUPFractionalStep
   use Mod_supf_elmope_hydro
   implicit none
   class(SUPFractionalStepProblem), target :: SUPFProblem   
   
   real(rp) :: eltemp(3) = 0.0_rp
   
   a=>SUPFProblem
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_elmope_hydro')     
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)      

   
   !Matrices Alloc
   call a%Memor%alloc((e%ndime+1),e%mnode,(e%ndime+1),e%mnode,elmat,'elmat','supf_elmope_hydro')
   call a%Memor%alloc((e%ndime+1),e%mnode,elrhs,'elrhs','supf_elmope_hydro')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_hydro')
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supf_elmope_hydro')
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','supf_elmope_hydro')
   call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','supf_elmope_hydro')
   call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','supf_elmope_hydro')
   call a%Memor%alloc(e%ndime,e%mnode,elrhu,'elrhu','supf_elmope_hydro')
   call a%Memor%alloc(1,e%mnode,elrhp,'elrhp','supf_elmope_hydro')
   
   !Other arrays alloc
   call a%Memor%alloc(e%ndime,elext,'elext','supf_elmope_hydro')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','supf_elmope_hydro')
   call a%Memor%alloc(e%mnode,testf,'testf','supf_elmope_hydro')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','supf_elmope_hydro')  
   call a%Memor%alloc(e%ndime,e%mnode,1,elvel,'elvel','supf_elmope_hydro')   

  
   
   !Statistics
   call a%InitStats
   
   !Hook
   call ProcHook%Initializations   
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
   
      !Initializations
      elmat=0.0_rp
      elrhs=0.0_rp
      elmuv=0.0_rp
      elrhu=0.0_rp
      wrmat1=0.0_rp
      elmpq=0.0_rp
      elmpv=0.0_rp
      elmuq=0.0_rp
      elrhp=0.0_rp
   
      !Gathers
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      
      !Hook
      call ProcHook%Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale)
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm                 
         
        !Hook
         call ProcHook%Interpolates                
     
         !Hook
         call ProcHook%PhysicalProp  
         !Physical Parameters
         call a%GetPhysicalParameters(imat,acden,acvis)   
         
         call ProcPointer%FreeSurface
         
         
         gpadv=0.0_rp         
         !Advection velocity norm
         call vecnor(gpadv,e%ndime,gpvno,2)
         !Compute the stability parameters      
         call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)         
         !Stabilization terms : -tau L*v         
         AGradV=0.0_rp         
         call nsm_ComputeTestf(e,acden,timom,AGradV,testf)         

         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                ! w(gp)*detjm*tau1
                  
         !Compute Elext, Temporal Derivatives, Repro...
         elext =0.0_rp   
         !Compute vector of external forces
         call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext)   
         
         !No temporal derivative
         eltemp = 0.0_rp
         
         !We solve only a stationary problem
         LHSDtinv=0.0_rp
         !Ingauss Elmats
         !Hook
         call ProcHook%InGaussElmats         
         ! Viscosity terms : we only consider mu*(grad v, grad u)         
         call elmvis(e,dvol,acvis,wrmat1)  
         !Compute contributions to elemental matrix : Block U,Q
         call nsm_elmbuq(e,timom,dvol,acden,LHSdtinv,AGradV,elmuq)           
         !Compute contributions to RHS : Block U
         call nsm_elmrhu(e,dvol,testf,elext,eltemp,elrhu)                       
         !Compute contributions to elemental matrix : Block V,P
         call nsm_elmbpv(e,dvol,testf,elmpv)
         !BLOCK P,Q : tau1*(graq q, grad p)
         call nsm_elmbpq(e,dvol*timom,elmpq)           
         !Compute contributions to RHS : Block P
         call nsm_elmrhp(e,timom,dvol,elext,eltemp,elrhp)    
        
         
         !Statistics
         call a%InGaussStats(acden,acvis,gpvno,chale,timom)
         
      enddo gauss_points
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer%PostGaussElmats  
      
      
      !Matrix composition
      ! Assembly wrmat1  in elmuv
      forall (idime = 1:e%ndime)
         elmuv(idime,1:e%pnode,idime,1:e%pnode) = elmuv(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      end forall  
      
      ! Momentum equation
      ! Assembly elmuv to elmat
      elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) + elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      ! Assembly elmpv to elmat
      elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) + elmpv(1:e%ndime,1:e%pnode,1,1:e%pnode)
      ! Assembly elrhu to elrhs
      elrhs(1:e%ndime,1:e%pnode) = elrhs(1:e%ndime,1:e%pnode) + elrhu(1:e%ndime,1:e%pnode)

      ! Continuity Equation
      ! Assembly elmuq to elmat
      elmat(e%ndime+1,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,1:e%ndime,1:e%pnode) + elmuq(1,1:e%pnode,1:e%ndime,1:e%pnode)
      ! Assembly elmpq to elmat
      elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) + elmpq(1,1:e%pnode,1,1:e%pnode)
      ! Assembly  elrhp to elrhs
      elrhs(e%ndime+1,1:e%pnode) = elrhs(e%ndime+1,1:e%pnode) + elrhp(1,1:e%pnode)     
      
      
      
      !Pre Assembly Modifications
      !Pointer
      call ProcPointer%PreAssembly       
      
      call nsm_rotdir(a,e,(e%ndime+1),elmat,elrhs)      
      !Dirichlet Boundary Conditions
      call supf_elmdir_press(a,e,elmat((e%ndime+1),:,(e%ndime+1),:),elrhs((e%ndime+1),:))
      !inside the subroutine
      ! php_elmdir(a,e,ndofn,ndofbc,ndofbcstart,currentbvess,elmat,elrhs)
      currentbvess=auxtens+1 
      bcstar=0_ip
      call php_elmdir(a,e,(e%ndime+1),e%ndime,bcstar,currentbvess,elmat,elrhs) 
   
      
      !Assembly
      call a%LinearSystemST%Assembly(e,elmat,elrhs)
      
   enddo elements
   
   call a%FinalizeStats
   
   !Hook
   call ProcHook%Finalizations
   
   !Matrices Dealloc
   call a%Memor%dealloc((e%ndime+1),e%mnode,(e%ndime+1),e%mnode,elmat,'elmat','supf_elmope_hydro')
   call a%Memor%dealloc((e%ndime+1),e%mnode,elrhs,'elrhs','supf_elmope_hydro')
   call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_hydro')
   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supf_elmope_hydro')
   call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','supf_elmope_hydro')
   call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','supf_elmope_hydro')
   call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','supf_elmope_hydro')
   call a%Memor%dealloc(e%ndime,e%mnode,elrhu,'elrhu','supf_elmope_hydro')
   call a%Memor%dealloc(1,e%mnode,elrhp,'elrhp','supf_elmope_hydro')
   
   !Other arrays Dealloc
   call a%Memor%dealloc(e%ndime,elext,'elext','supf_elmope_hydro')    
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','supf_elmope_hydro')
   call a%Memor%dealloc(e%mnode,testf,'testf','supf_elmope_hydro')  
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','supf_elmope_hydro')  
   call a%Memor%dealloc(e%ndime,e%mnode,1,elvel,'elvel','supf_elmope_hydro')   
   
   

   !ElementDealloc
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','supf_elmope_hydro')
   
   
end subroutine
