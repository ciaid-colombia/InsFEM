module Mod_supf_EndElmope_s
   use Mod_Mesh
   use Mod_Memor
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_NavierStokesElement
   use Mod_ThreeFieldElement   
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_nsm_elmdir
   use Mod_SupExacso 
   use Mod_SupOperations
   implicit none
   
   class(SUPFractionalStepProblem), pointer :: a
   type(SupExacso) :: exacso    
   
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: ComputeAdvectionVelocity => NULL()
      procedure(), NOPASS, pointer :: ExternalForces    => NULL()
      procedure(), NOPASS, pointer :: FEResPro => NULL()
      procedure(), NOPASS, pointer :: TemporalDerivatives => NULL()  
      procedure(), NOPASS, pointer :: PreAssembly => NULL()     
      
   end type
   type(PPointer) :: ProcPointer   
   
   
   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: PreLoop => NULL()
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: ElmatsToZero => NULL()
      procedure(), NOPASS, pointer :: Gathers => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: Interpolates => NULL()
      procedure(), NOPASS, pointer :: Elext => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: AssemblyEndite => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PostLoop => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()     
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
      
   real(rp), allocatable :: elvel(:,:,:),elsig(:,:,:)
   real(rp), allocatable :: elpre(:,:)
   real(rp), allocatable :: elext(:),elextC(:),elextS(:)
   real(rp), allocatable :: AGradV(:), testf(:) 
   real(rp), allocatable :: gpvel(:,:),gpadv(:),gpsig(:,:)  
      
   real(rp)    :: chale(2),timom, tidiv,dvolt0,dvolt1,dvolt2,gprhs(3)
   integer(ip) :: nmean,auxiter,aux12,auxSU                       ! Stabilization
   real(rp)    :: acden,acvis,vista,beta,auxVE,auxG,lambda
   real(rp)    :: reyno
   real(rp)    :: gpvno,gppre(1,1),divvel,dummr
   
   !Residual Projections   
   real(rp), allocatable :: elrep(:,:),gprep(:)
   real(rp), allocatable :: grpre(:,:),grvel(:,:),grsig(:,:)
   integer(ip)           :: ipoin,ibopo,npoin
   real(rp), pointer     :: exnor(:,:)
   
   !Oss Split
   real(rp), allocatable :: gpdivs(:),elrepdivs(:,:)   
   !Discontinuity Capturing
   real(rp), allocatable :: gpgrad(:),elrepGrad(:,:)   
   !Gradient Projection
   real(rp), allocatable :: gpsgrad(:),elrepSGrad(:,:)   
   
   integer(ip) :: idime,itime,igaus,auxtens,auxGrad,auxSGrad,auxdim,auxrepro
   
   
   !Level Set
   integer(ip)              :: inode
   integer(ip)              :: elemStatus,ngauss_minus,ngauss_plus,ngaus_total,poinStatus
   real(rp), allocatable    :: weigp(:),xloc(:,:)     
   
   
   !todo multy materials
   integer(ip) :: imat=1
   
   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"
   
   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty
      
      !External Procedures
      procedure() :: NULLSUB
      
      call ResetProcedureComposition
      !Deafults
      !Pointers
      ProcPointer%ComputeAdvectionVelocity => ComputeAdvectionVelocity
      ProcPointer%ExternalForces  => NULLSUB  
      ProcPointer%FEResPro        => FEResPro2dVE 
      !Viscoelastic Part
      ProcPointer%TemporalDerivatives => sup_temporalpart     
      !PreAssembly to use in free-surface and in enriched elements
      ProcPointer%PreAssembly => NULLSUB       
      
      !Hooks
      ProcHook%PreLoop            => NULLSUB
      ProcHook%Initializations    => NULLSUB
      ProcHook%OnIeltyChange      => NULLSUB
      ProcHook%PreGauss           => NULLSUB
      ProcHook%ElmatsToZero       => NULLSUB
      ProcHook%Gathers            => NULLSUB
      ProcHook%InGauss            => NULLSUB
      ProcHook%Interpolates       => NULLSUB
      ProcHook%Elext              => NULLSUB
      ProcHook%InGaussElmats      => NULLSUB
      ProcHook%AssemblyEndite     => NULLSUB
      ProcHook%Finalizations      => NULLSUB
      ProcHook%PostLoop           => NULLSUB
      ProcHook%PhysicalProp       => NULLSUB      

      
      !MaterialProperties is called Always
      call ConcatenateProcedures(ProcHook%PhysicalProp,MaterialProperties)    
      
      !Advection velocity
      if (a%kfl_advec == 0) then
         ProcPointer%ComputeAdvectionVelocity => NULLSUB
      endif
      
      !-----------------------------------------------------------
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (a%MatProp(imat)%lawvi<0 .and. kfl_nonlinear == 1) then
           call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
      endif
      
      auxdim=e%ndime

      if(e%ndime==3)then
         ProcPointer%FEResPro => FEResPro3dVE
      end if      
      
      !Residual Projection
      if (a%kfl_repro == 1 .or. a%kfl_repro == 2) then
         call ConcatenateProcedures(ProcHook%PreLoop,PreLoopRep)
         call ConcatenateProcedures(ProcHook%Initializations,InitializationsRep)
         call ConcatenateProcedures(ProcHook%ElmatsToZero,ElmatsToZeroRep)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGradients)
         
         call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeRHSRep)
         if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeRHSNonLinearRep)
         call ConcatenateProcedures(ProcHook%InGaussElmats,GaussPointAssemblyRep)         
         call ConcatenateProcedures(ProcHook%AssemblyEndite,AssemblyResidual)     
         call ConcatenateProcedures(ProcHook%Finalizations,FinalizationsRep)
         call ConcatenateProcedures(ProcHook%PostLoop,ResidualBoundary)
         call ConcatenateProcedures(ProcHook%PostLoop,ProjectResidual)
      endif

      !Discontinuity Capturing
      if (a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 1) then
         call ConcatenateProcedures(ProcHook%PreLoop,PreLoopDiscCapt)
         call ConcatenateProcedures(ProcHook%Initializations,InitDiscontinuity)
         call ConcatenateProcedures(ProcHook%ElmatsToZero,ElmatsToZeroRepGrad)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGradients)
         call ConcatenateProcedures(ProcHook%InGaussElmats,GradientToVector)
         call ConcatenateProcedures(ProcHook%InGaussElmats,GaussPointAssemblyRepGrad)         
         call ConcatenateProcedures(ProcHook%AssemblyEndite,AssemblyGradient)
         call ConcatenateProcedures(ProcHook%Finalizations,FinDiscontinuity)
         call ConcatenateProcedures(ProcHook%PostLoop,SmoothGradient)
      endif        
      
      auxrepro=a%kfl_repro

      !Split-Oss
      if (a%MatProp(imat)%lawvi<0) then
         if(a%kfl_repro == 2 .or. a%kfl_repro == 3)then
            call ConcatenateProcedures(ProcHook%PreLoop,PreLoopOssSplit)
            call ConcatenateProcedures(ProcHook%Initializations,InitOssSplit)
            call ConcatenateProcedures(ProcHook%ElmatsToZero,ElmatsToZeroOssSplit)
            call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGradients)
            call ConcatenateProcedures(ProcHook%InGaussElmats,termToVector)
            call ConcatenateProcedures(ProcHook%InGaussElmats,GaussPointAssemblyOssSplit)
            call ConcatenateProcedures(ProcHook%AssemblyEndite,AssemblyOssSplit)         
            call ConcatenateProcedures(ProcHook%Finalizations,FinOssSplit)
            call ConcatenateProcedures(ProcHook%PostLoop,ResidualBoundaryOssSplit)         
            call ConcatenateProcedures(ProcHook%PostLoop,SmoothOssSplit) 
         end if         
      endif        
      
       !ExactSol
      if (a%kfl_exacs/=0) then  
         ProcPointer%ExternalForces => ExactSolutionForces
      endif      
      
       !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
      

      !--------------------------------------------------------------------
      !Level set
      if(a%kfl_colev==1)then
         call ConcatenateProcedures(ProcHook%Initializations,AllocLevelsetTF) 
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocLevelSetTF)
         call PrependProcedure(ProcHook%PreGauss,CutelementsTF)
         call PrependProcedure(ProcHook%PhysicalProp,FluidPropertiesTF) 
!          if(a%kfl_fsurf==1)then
! !             call ConcatenateProcedures(ProcPointer%PreAssembly,FreeSurfMatsTF)
!          end if 
         
      end if 
      
      
      
   end subroutine
   
   !-------------------------------------------------------------------
   !AdvectionVelocity
   subroutine ComputeAdvectionVelocity
      implicit none
      
      gpadv = gpvel(:,1)
   end subroutine
   
   !Temporal Derivatives
   subroutine sup_temporalpart
      implicit none
      
      call sup_TimeIntegrationToElext(e,auxtens,integrator,auxVE,a%dtinv,gpsig,elextS)  
          
   end subroutine   
   
   !-------------------------------------------------------------------
   !FOR NON-LINEAR ELEMENTS
   subroutine InGaussNonLinear
      implicit none
      
      call e%elmder
      call e%elmhes
   end subroutine
   
   !-------------------------------------------------------------------
   !Physical Properties
   subroutine MaterialProperties 
      implicit none 
      !-----------------------------------------------------------
      !Physical Properties
      acvis  = a%MatProp(imat)%LawViParam(1) 
      acden  = a%MatProp(imat)%densi
      beta   = a%MatProp(imat)%LawViParam(2)
         
      !incremental scheme for Weissenberg number
      call a%IncrementalLambda(imat,a%MatProp(imat)%LawViParam(3), lambda)  
      auxG   = a%MatProp(imat)%LawViParam(4)/((1.0_rp-beta) + 0.00001_rp)    

   end subroutine    
   
   
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   !-------------------------------------------------------------------
   !FOR RESIDUAL PROJECTION
   subroutine PreLoopRep
      implicit none
      
      a%repro = 0.0_rp
   end subroutine
   
   subroutine InitializationsRep
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(a%ndofn,e%mnode,elrep,'elrep','supm_EnditeElmope')
      call a%Memor%alloc(a%ndofn,gprep,'gprep','supm_EnditeElmope')
   end subroutine

   subroutine FinalizationsRep
      implicit none
        
      !Matrices dealloc
      call a%Memor%dealloc(a%ndofn,e%mnode,elrep,'elrep','supm_EnditeElmope')
      call a%Memor%dealloc(a%ndofn,gprep,'gprep','supm_EnditeElmope')
   end subroutine

   subroutine ElmatsToZeroRep
      implicit none
      
      elrep = 0.0_rp
   end subroutine
   
   subroutine InterpolateGradients
      implicit none
      call e%gradient(1   ,elpre(:,1)  ,grpre)       !Press. gradient
      call e%gradient(e%ndime,elvel(:,:,1),grvel)    !Vel. gradient
      call e%gradient(auxtens,elsig(:,:,1),grsig)    !Sig. gradient  
   end subroutine
   

   
   !Force term   
   !Non Exact Solution Case
   subroutine  nsiForces
      implicit none    
      !Compute vector of external forces
      call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext)  
   end subroutine   
        
   subroutine  ExactSolutionForces    
      use Mod_Mesh    
      implicit none      
      real(rp)    :: gpcod(e%ndime)
      !Interpolate
      call e%interpg(e%ndime,e%elcod,gpcod)
      !Compute vector of external forces      
      call exacso%sup_ComputeSolution(e%ndime,gpcod,a%ctime,a%LogFormulation,a)
      !call exacso%sup_GetForce(e%ndime,a%LogFormulation,a%kfl_LCR,elext,elextC,elextS,a)        
   end subroutine    
      

   subroutine ComputeRHSRep
      implicit none    
         call ProcPointer%FEResPro        
   end subroutine   
   
   subroutine FEResPro2dVE
      implicit none      
      !Compute   

      
      call supm_elmrfeVE_oto(e,LHSdtinv,acden,acvis,beta,lambda,auxVE,auxG,a%PTT_model,gpadv, &
        gpvel,gpsig,grpre,grvel,grsig,auxtens,elext,elextC,elextS,gprep)      
        
        
   end subroutine   
   

   subroutine FEResPro3dVE
      implicit none      
      !Compute
     call supm_elmrfeVE_oto(e,LHSdtinv,acden,acvis,beta,lambda,auxVE,auxG,a%PTT_model,gpadv, &
        gpvel,gpsig,grpre,grvel,grsig,auxtens,elext,elextC,elextS,gprep)                
        
   end subroutine 

   
   subroutine GaussPointAssemblyRep
      implicit none  
      integer(ip) :: elemStatus
      
      !this is crucial when we have two fluids and enriched pressure
      if(a%kfl_colev==1)then
         call a%CutMesh%GetElementType(ielem,elemStatus)
         
         if(elemStatus ==0)then 
            gprep=0.0_rp
         end if         
      
      end if       
      
      call supm_elmrep(e,dvol,a%ndofn,gprep,elrep)
  
   end subroutine
   
   subroutine ComputeRHSNonLinearRep
      implicit none   
      
      call supm_elmrfe_oto_nonlinear(e,acvis,auxtens,beta,elvel,gprep)

   end subroutine
   
   subroutine AssemblyResidual
      implicit none
     
      a%repro(:,e%lnods(1:e%pnode)) = a%repro(:,e%lnods(1:e%pnode)) + elrep(:,1:e%pnode)
           
   end subroutine
   
   subroutine ResidualBoundary
   
   call a%Mesh%Getnpoin(npoin)      
      
      do ipoin=1,npoin      
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if(ibopo>0)then 

            do idime=1,auxtens+e%ndime+1
               a%repro(idime,ipoin) = 0.0_rp
            end do           
           
         end if
      end do    
         
   end subroutine

   subroutine ProjectResidual
      implicit none
      
      call a%Project(a%ndofn,a%repro) 
      
   end subroutine
   
   !------------------------------------------------------
   !OSS Split

   subroutine PreLoopOssSplit
      implicit none
      
      a%reproSDiv   = 0.0_rp !Div(S)
      
   end subroutine   
   
   subroutine InitOssSplit
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%ndime,gpdivs, 'gpdivs','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrepdivs, 'elrepdivs','supm_EnditeElmope')
   
   end subroutine
   
   subroutine ElmatsToZeroOssSplit
      implicit none
      
      elrepdivs = 0.0_rp   
      
   end subroutine   

   subroutine FinOssSplit
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(e%ndime,gpdivs, 'gpdivs','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepdivs, 'elrepdivs','supm_EnditeElmope')    
   end subroutine
      

   subroutine termToVector
      implicit none
      
      call sup_GetDivS(e%ndime,auxtens,grsig,gpdivs)
      
   end subroutine
   
   subroutine GaussPointAssemblyOssSplit
      implicit none  
      integer(ip) :: elemStatus
      
      !this is crucial when we have two fluids and enriched pressure
      if(a%kfl_colev==1)then
         call a%CutMesh%GetElementType(ielem,elemStatus)
         
         if(elemStatus ==0)then 
            gpdivs=0.0_rp
         end if         
      
      end if      
      
      call supm_elmrep(e,dvol,e%ndime,gpdivs,elrepdivs)
      
   end subroutine   
   
   subroutine AssemblyOssSplit
      implicit none    
      
      a%reproSDiv(:,e%lnods(1:e%pnode)) = a%reproSDiv(:,e%lnods(1:e%pnode)) + elrepdivs(:,1:e%pnode)

   end subroutine
   
   subroutine ResidualBoundaryOssSplit
   
   
         
   end subroutine
      
   
   subroutine SmoothOssSplit
      implicit none
      
      call a%Project(e%ndime,a%reproSDiv)  
  
   end subroutine   

   !------------------------------------------------------
   !Discontinuity Capturing
   subroutine PreLoopDiscCapt
      implicit none
      
      a%reproGrad = 0.0_rp
      
   end subroutine   
   
   subroutine InitDiscontinuity
      implicit none
      
      auxGrad=e%ndime*auxtens
      !Matrices alloc
      call a%Memor%alloc(auxGrad,gpgrad,'gpgrad','supm_EnditeElmope')
      call a%Memor%alloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supm_EnditeElmope')         
   
   end subroutine
   
   subroutine ElmatsToZeroRepGrad
      implicit none
      
      elrepGrad = 0.0_rp
      
   end subroutine   

   subroutine FinDiscontinuity
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(auxGrad,gpgrad,'gpgrad','supm_EnditeElmope')
      call a%Memor%dealloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supm_EnditeElmope')               
   end subroutine
      

   subroutine GradientToVector
      implicit none
      
      call sup_GradStress(e%ndime,auxtens,auxGrad,grsig,gpgrad)
      
   end subroutine
   
   subroutine GaussPointAssemblyRepGrad
      implicit none   
      
      call supm_elmrep(e,dvol,auxGrad,gpgrad,elrepGrad)
      
   end subroutine   
   
   subroutine AssemblyGradient
      implicit none
      
      a%reproGrad(:,e%lnods(1:e%pnode)) = a%reproGrad(:,e%lnods(1:e%pnode)) + elrepGrad(:,1:e%pnode)
      
   end subroutine
   
   subroutine SmoothGradient
      implicit none
      
      call a%Project(auxGrad,a%reproGrad) 
   
   end subroutine    
   
   !-------------------------------------------------------
   !LevelSet
   
   subroutine AllocLevelsetTF
      implicit none   
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,xloc,'xloc','nsm_elmope') 
   
   end subroutine
   
   subroutine DeallocLevelSetTF
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,xloc,'xloc','nsm_elmope') 
      
   
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
      integer(ip)  :: elemStatus      
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==1)then
         imat=1
         acden=a%MatProp(imat)%densi
         acvis=a%MatProp(imat)%visco           
      elseif(elemStatus==-1)then
         imat=2
         acden=a%MatProp(imat)%densi
         acvis=a%MatProp(imat)%visco           
      elseif(elemStatus ==0)then
         
         ngauss_minus=0
         ngauss_plus=0
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
        
         if(e%igaus<=ngauss_minus)then
            imat=2
            acden=a%MatProp(imat)%densi
            acvis=a%MatProp(imat)%visco               
         elseif(e%igaus>ngauss_minus)then
            imat=1
            acden=a%MatProp(imat)%densi
            acvis=a%MatProp(imat)%visco               
         end if
         
      end if
   
   end subroutine
   
   subroutine FreeSurfMatsTF
      implicit none 
      integer(ip) :: elemStatus
      
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==-1)then
         
         !Split OSS
         if (a%kfl_repro == 1 .or. a%kfl_repro == 2) then 
            elrep(:,:)     = 0.0_rp 
         end if
         
         !Slipt total and Split OSS
         elrepdivs(:,:) = 0.0_rp
         
         !monolithic scheme works without the following sentence 
         
!          !Discontinuity Capturing
!          if (a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 1) then
!             elrepGrad(:,:) = 0.0_rp
!          end if
      end if
   end subroutine        
   
   
   
   
   
   
end module 



subroutine supf_EndElmope_s(SUPFProblem)
   use Mod_supf_EndElmope_s
   implicit none
   class(SUPFractionalStepProblem), target :: SUPFProblem
   a=>SUPFProblem
   
   !Return if there is nothing to be done
   if (a%kfl_repro == 0 .and. a%kfl_shock==0) return
   
   !Element allocation
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'supm_EnditeElmope')    
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call ProcHook%PreLoop
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   call a%Mesh%GetNelem(nelem)
   
   
   auxtens=(e%ndime-1)*(e%ndime-1)+2
   !Other arrays alloc
   call a%Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','supm_EnditeElmope')
   call a%Memor%alloc(      e%mnode,a%ncomp-1,elpre,'elpre','supm_EnditeElmope')
   call a%Memor%alloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supm_EnditeElmope')   
   call a%Memor%alloc(e%ndime,elext,'elext','supm_EnditeElmope')
   call a%Memor%alloc(1,elextC,'elextC','supm_EnditeElmope')  
   call a%Memor%alloc(auxtens,elextS,'elextS','supm_EnditeElmope')   
   call a%Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supm_EnditeElmope')
   call a%Memor%alloc(auxtens,a%ncomp-1,gpsig,'gpsig','supm_EnditeElmope')   
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','supm_EnditeElmope')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','supm_EnditeElmope')
   call a%Memor%alloc(e%mnode,testf,'testf','supm_EnditeElmope')
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','supm_EnditeElmope')
   call a%Memor%alloc(1_ip     ,e%ndime,grpre,'grpre','supm_EnditeElmope')
   call a%Memor%alloc(auxtens,e%ndime,grsig,'grsig','supm_EnditeElmope')
   
   
   !Hook
   call ProcHook%Initializations
   
   !do itest = 1,100
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)    
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
      
      !Elmats to Zero
      call ProcHook%ElmatsToZero
      
      !Gathers
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(auxtens,elsig(:,:,1),a%sigma(:,:,1))      
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      !viscoelastic case
      call e%gather(auxtens,elsig(:,:,2),a%sigma(:,:,3))      
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(e%ndime,elvel(:,:,itime),a%veloc(:,:,itime+1))
         call e%gather(auxtens,elsig(:,:,itime),a%sigma(:,:,itime+1))          
      enddo
      call e%gather(1_ip   ,elpre(1:e%pnode,1),a%press(1:e%pnode,1))
      
      !Hook
      call ProcHook%Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale)
      
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1))
         call e%interpg(auxtens,elsig(:,:,1),gpsig(:,1))         
         call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2))
         call e%interpg(auxtens,elsig(:,:,2),gpsig(:,2)) !i-1 iteration           
         do itime = 3,nsteps ! Time bdf2 and others
            call e%interpg(e%ndime,elvel(:,:,itime),gpvel(:,itime))
            call e%interpg(auxtens,elsig(:,:,itime),gpsig(:,itime))            
         enddo
         !Hook
         call ProcHook%Interpolates
         
         !Physical Parameters
         call a%GetPhysicalParameters(imat,acden,acvis)          
     
         !Hook
         call ProcHook%PhysicalProp     
         
         !Advection velocity      
         call ProcPointer%ComputeAdvectionVelocity
         
         !Advection velocity norm
         call vecnor(gpadv,e%ndime,gpvno,2)
         
         !Compute a·grad(V)
         call ComputeAGradV(e,gpadv,AGradV)
         
         !Compute Elext, Temporal Derivatives, Repro...
         
         !Compute Elext, Temporal Derivatives
         elext = 0.0_rp
         elextC= 0.0_rp
         elextS= 0.0_rp         

         call ProcPointer%TemporalDerivatives  
         
         !Compute vector of external forces
         call ProcPointer%ExternalForces   
                         
         !InGaussElmats
         !Hook
         call ProcHook%InGaussElmats
      enddo gauss_points
      
      
      !Pre Assembly Modifications
      !Pointer
      call ProcPointer%PreAssembly            
      
      
      !Assembly Endite
      !Hook
      call ProcHook%AssemblyEndite
      
      
   enddo elements
   !enddo
   
   !Hook
   call ProcHook%Finalizations
   
   !Other arrays alloc
   call a%Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','supm_EnditeElmope')
   call a%Memor%dealloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supm_EnditeElmope')   
   call a%Memor%dealloc(      e%mnode,a%ncomp-1,elpre,'elpre','supm_EnditeElmope')   
   call a%Memor%dealloc(e%ndime,elext,'elext','supm_EnditeElmope')
   call a%Memor%dealloc(1,elextC,'elextC','supm_EnditeElmope')  
   call a%Memor%dealloc(auxtens,elextS,'elextS','supm_EnditeElmope')    
   call a%Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supm_EnditeElmope')
   call a%Memor%dealloc(auxtens,a%ncomp-1,gpsig,'gpsig','supm_EnditeElmope')   
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','supm_EnditeElmope')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','supm_EnditeElmope')
   call a%Memor%dealloc(e%mnode,testf,'testf','supm_EnditeElmope')
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','supm_EnditeElmope')
   call a%Memor%dealloc(1_ip   ,e%ndime,grpre,'grpre','supm_EnditeElmope') 
   call a%Memor%dealloc(auxtens,e%ndime,grsig,'grsig','supm_EnditeElmope')   
   
   !Operations to be done after the Elemental Loop
   !Hook
   call ProcHook%PostLoop
   
   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,a%EndLoopQuadrature,'supm_EnditeElmope')

   
end subroutine











