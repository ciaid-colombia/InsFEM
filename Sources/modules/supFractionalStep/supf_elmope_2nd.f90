module Mod_supf_elmope_2nd
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_Mesh
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_NavierStokesElement
   use Mod_ThreeFieldElement   
   use Mod_SUPF_Element   
   use Mod_ConvectiveElement
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_sup_elmdir   
   use Mod_nsm_elmdir
   use Mod_php_Elmdir
   use Mod_SupExacso    
   use Mod_SupOperations

   implicit none   
   
   class(SUPFractionalStepProblem), pointer :: a
   type(SupExacso) :: exacso     
   real(rp), parameter :: zensi = 0.0_rp   
   
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(), NOPASS, pointer :: ExternalForces => NULL()
      procedure(), NOPASS, pointer :: ComputeAdvectionVelocity  => NULL()
      procedure(), NOPASS, pointer :: ResProRHS => NULL()
      procedure(), NOPASS, pointer :: SplitOSSComponents => NULL()      
      procedure(), NOPASS, pointer :: ExtrapolationTerms => NULL()
      procedure(), NOPASS, pointer :: TauConstitutive => NULL()
      !Viscoelastic pointers      
      procedure(), NOPASS, pointer :: ViscoGalerkinMatrix => NULL()
      procedure(), NOPASS, pointer :: ViscoEstabMatrix => NULL()    
      !Discontinuity Capturing
      procedure(), NOPASS, pointer :: DiscontinuityCapturing => NULL()
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
      procedure(), NOPASS, pointer :: Testf => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp  => NULL()
      
      
    end type
   type(PHook) :: ProcHook 
   
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ielem,nelem
  
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:),elrhs2(:,:)
   real(rp)              :: dvol

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv

   integer(ip)           :: ielty0 = 0 !Previous element type
      
   real(rp), allocatable :: elvel(:,:,:)
   real(rp), allocatable :: elsig(:,:,:)
   real(rp), allocatable :: elext(:),elextS(:)
   real(rp), allocatable :: elmst(:,:,:,:)
   real(rp), allocatable :: AGradV(:), testf(:) 
   real(rp), allocatable :: gpvel(:,:),gpadv(:)
   real(rp), allocatable :: gpsig(:,:),grsig(:,:)
   !Residual Projections   
   real(rp), allocatable :: elrep(:,:),gprep(:)
   real(rp), allocatable :: grpre(:,:),grvel(:,:)  
   !Split Oss
   real(rp), allocatable :: gpdivs(:),elrepdivs(:,:) 
   integer(ip) :: idime,itime,igaus,itest,auxntens,auxndim,auxGrad,kdime,ldime,aux12
   !Discontinuity Capturing
   real(rp), allocatable :: elrepGrad(:,:),gpgrad(:),grsigRP(:,:),grsigRPO(:,:) 
   real(rp)    :: cshock(2) !shock parameters
   
   real(rp)    :: chale(2),timom,tisig,dvolt0,dvolt1,dvolt3,gprhs(3)
   integer(ip) :: nmean,auxtens,auxPTT
   real(rp)    :: acden,acvis,beta,lambda,auxVE,auxG,elextC(1)
   real(rp)    :: reyno
   real(rp)    :: gpvno,gpvno2,divvel,dummr,auxp1,auxp2,auxp3,auxp4,facdisc,kdisc  
   integer(ip) :: currentbvess,bcstar,auxconv,auxtrac1,auxtrac2,auxGie
   
   !Level Set
   integer(ip)              :: inode,ipoin
   integer(ip)              :: elemStatus,ngauss_minus,ngauss_plus,ngaus_total,poinStatus
   real(rp), allocatable    :: weigp(:),xloc(:,:)
   !Laplacian Term
   real(rp), allocatable :: wrmat1(:,:)   
   
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
      ProcPointer%ComputeAdvectionVelocity => ComputeAdvectionVelocity
      ProcPointer%ExternalForces  => NULLSUB
      ProcPointer%ExtrapolationTerms  => NULLSUB
      ProcPointer%ResProRHS       => ResProRHS2d      

      !Viscoelastic Part
      ProcPointer%TauConstitutive => TauElastic      
      ProcPointer%ViscoGalerkinMatrix => elmatGalerk
      ProcPointer%ViscoEstabMatrix => elmatEstab 
      ProcPointer%SplitOSSComponents => ResProSOSS       

      
      !Discontinuity Capturing
      ProcPointer%DiscontinuityCapturing => discontinuity2d  
      
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
      ProcHook%Testf           => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      ProcHook%PhysicalProp    => NULLSUB 
      
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_elmope_2nd')     
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      
      
      !MaterialProperties is called Always
      call ConcatenateProcedures(ProcHook%PhysicalProp,MaterialProperties)      
      
      
      !Advection velocity
      if (a%kfl_advec == 0) then
         ProcPointer%ComputeAdvectionVelocity => NULLSUB
      endif
            
      !-----------------------------------------------------------
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
      endif
      
      !-----------------------------------------------------------
      !Three-dimensional case
      if(e%ndime==3)then    
         ProcPointer%ViscoGalerkinMatrix   => elmatGalerk3d
         ProcPointer%ViscoEstabMatrix      => elmatEstab3d         
         ProcPointer%DiscontinuityCapturing => discontinuity3d 
         ProcPointer%SplitOSSComponents     => ResProSOSS3d
      endif      
 
      
      !-----------------------------------------------------------
      !ResidualProjection
      !OSS
      if (a%kfl_repro == 1 .or. a%kfl_repro==0) then
         call runend('ASGS and OSS in Viscoelastic Constitutivefractional dont work yet')
       endif
       
      !Split OSS terms  
      if(a%kfl_repro == 2 .or. a%kfl_repro == 3)then           
         call ConcatenateProcedures(ProcHook%Initializations,AllocOssSplit)
         call ConcatenateProcedures(ProcHook%Gathers,GatherOssSplit)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateOssSplit)
         call ConcatenateProcedures(ProcHook%InGaussElmats,IngausSplitOss)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocOssSplit)  
      end if       
      
      !-----------------------------------------------------------      
      !OSS constitutive
      if(a%kfl_repro == 2)then                   
         if(e%ndime == 3)then
            ProcPointer%ResProRHS       => ResProRHS3d
         end if        
         call ConcatenateProcedures(ProcHook%Initializations,AllocRep)
         call ConcatenateProcedures(ProcHook%Gathers,GatherRep)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateRep)
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsRep)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocRep)             
      end if
      
      !-----------------------------------------------------------
      !Discontinuity Capturing
      if(a%kfl_shock == 1)then
         call ConcatenateProcedures(ProcHook%Initializations,AllocRepGrad)
         call ConcatenateProcedures(ProcHook%Gathers,GatherRepGrad)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateRepGrad)
         call ConcatenateProcedures(ProcHook%Interpolates,GradRespro)         
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsRepGrad)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocRepGrad)
      end if      
      
      !-----------------------------------------------------------
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
         if(a%kfl_fsurf==1)then
            if(a%kfl_fsurfLapla==1 .or. a%kfl_fsurfLapla==2)then
               call PrependProcedure(ProcPointer%PreAssembly,ElmatsToLapla)
            end if
            
            if(a%kfl_fsurfLapla==2)then
               call ConcatenateProcedures(ProcHook%Initializations,AllocStokes) 
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocStokes)
               call ConcatenateProcedures(ProcHook%PreGauss,InitRHSStokes)
               call ConcatenateProcedures(ProcHook%InGaussElmats,RHSStokes)
               call ConcatenateProcedures(ProcPointer%PreAssembly,ElmatsToStokes)
            end if
            
            call ConcatenateProcedures(ProcPointer%PreAssembly,FreeSurfMatsTF)            
         end if 
         
      end if   




      
   end subroutine   
   
   !AdvectionVelocity
   subroutine ComputeAdvectionVelocity
      implicit none
      
      gpadv = gpvel(:,1)
   end subroutine 
   

   !----------------------------------------------------------
   !PostGauss Matrices
   subroutine PostGaussElmats
      implicit none      

      
   end subroutine   
   
   !-------------------------------------------------------------------
   !Physical Properties
   subroutine MaterialProperties 
      implicit none 
      !----------------------------------------------------------------
      !Incremental strategy (Continuation method in terms of lambda)
      acvis  = a%MatProp(imat)%LawViParam(1) 
      beta   = a%MatProp(imat)%LawViParam(2)
         
      !incremental scheme for Weissenberg number
      call a%IncrementalLambda(imat,a%MatProp(imat)%LawViParam(3), lambda)                
         
      auxPTT = 1_ip
      auxG   = a%MatProp(imat)%LawViParam(4)/((1.0_rp-beta) + 0.00001_rp)
      auxVE  = lambda/(2.0_rp*acvis)  
      
      if(a%MatProp(imat)%lawvi==-3)then
         auxPTT = 0_ip
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
   
   !-------------------------------------------------------------------
   !FOR RESIDUAL PROJECTION
   subroutine AllocRep
      implicit none

      auxndim=a%ndofn
      call a%Memor%alloc(a%ndofn,e%mnode,elrep,'elrep','supf_elmope_2nd')
      call a%Memor%alloc(a%ndofn,gprep,'gprep','supf_elmope_2nd')
      call a%Memor%alloc(1_ip     ,e%ndime,grpre,'grpre','supf_elmope_2nd')
   end subroutine

   subroutine DeallocRep
      implicit none

      call a%Memor%dealloc(a%ndofn,e%mnode,elrep,'elrep','supf_elmope_2nd')
      call a%Memor%dealloc(a%ndofn,gprep,'gprep','supf_elmope_2nd')
      call a%Memor%dealloc(1_ip     ,e%ndime,grpre,'grpre','supf_elmope_2nd')
   end subroutine

   subroutine GatherRep
      implicit none
      
      call e%gather(a%ndofn,elrep,a%repro)
   end subroutine

   subroutine InterpolateRep
      implicit none
      
      !Interpolate
      call e%interpg(a%ndofn,elrep,gprep)
   end subroutine

   subroutine InGaussElmatsRep
      implicit none

     call ProcPointer%ResProRHS
 
   end subroutine
   
   
   subroutine ResProRHS2d
      implicit none

      !Compute contributions to RHS : Block C
      call supf_elmrhcVES_oss(e,auxGie,auxVE,auxG,auxPTT,tisig,acvis,dvol,auxtens,AGradV,grvel,gpsig,gprep,elrhs) 

   end subroutine  
  
   subroutine ResProRHS3d
      implicit none

      !Compute contributions to RHS : Block C
      call supf_elmrhcVES_oss3d(e,auxGie,auxVE,auxG,auxPTT,tisig,acvis,dvol,auxtens,AGradV,grvel,gpsig,gprep,elrhs) 
   
   end subroutine
   
   !------------------------------------------------------------------------------------------
   !Discontinuity Capturing
   
   subroutine AllocRepGrad
      implicit none

      cshock(1)=a%shock(1)
      cshock(2)=a%shock(2)
      
      auxGrad=e%ndime*auxtens
      
      call a%Memor%alloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supf_elmope_2nd')
      call a%Memor%alloc(auxGrad,gpgrad,'gpgrad','supf_elmope_2nd')
      call a%Memor%alloc(auxtens,e%ndime,grsigRP,'grsigRP','supf_elmope_2nd')
      call a%Memor%alloc(auxtens,e%ndime,grsigRPO,'grsigRPO','supf_elmope_2nd')   
      
   end subroutine

   subroutine DeallocRepGrad
      implicit none

      call a%Memor%dealloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supf_elmope_2nd')
      call a%Memor%dealloc(auxGrad,gpgrad,'gpgrad','supf_elmope_2nd')
      call a%Memor%dealloc(auxtens,e%ndime,grsigRP,'grsigRP','supf_elmope_2nd')      
      call a%Memor%dealloc(auxtens,e%ndime,grsigRPO,'grsigRPO','supf_elmope_2nd')       
      
   end subroutine
   
   subroutine GatherRepGrad
      implicit none
      
      call e%gather(auxGrad,elrepGrad,a%reproGrad)    
      
   end subroutine

   subroutine InterpolateRepGrad
      implicit none
      
      !Interpolate
      call e%interpg(auxGrad,elrepGrad,gpgrad)
   end subroutine
   
   subroutine GradRespro
      implicit none
      
      !Interpolate
      call sup_grsigRP(e%ndime,auxtens,auxGrad,gpgrad,grsigRP)      
            
      grsigRPO=0.0_rp
      !Orthogonal Projection
      do kdime=1,auxtens
         do ldime=1,e%ndime
            grsigRPO(kdime,ldime)= grsigRPO(kdime,ldime) + (grsig(kdime,ldime)-grsigRP(kdime,ldime))
         end do
      enddo
      
   end subroutine
   
   subroutine InGaussElmatsRepGrad
      use  Mod_sup_Kdiscont
      implicit none
      
      call sup_Kdiscont(e,cshock,auxtens,acvis,lambda,chale,grsig,grsigRPO,grsigRP,grvel,gpvno2,facdisc,kdisc)    
     
      if(a%kfl_colev==0)then
         a%viscarray(ielem)%a(e%igaus) = Kdisc
      end if 
      !Discontinuity capturing
      call ProcPointer%DiscontinuityCapturing

   end subroutine
   
   subroutine discontinuity2d
      implicit none   
      
      call supm_elmbstDC(e,auxtens,kdisc,dvol,elmat)   
   
   end subroutine
   
   subroutine discontinuity3d
      implicit none   
      call supm_elmbstDC(e,auxtens,kdisc,dvol,elmat)
   end subroutine
      

   !------------------------------------------------------------------------------------------------------------
   !Split-OSS
   
   subroutine AllocOssSplit
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%ndime,gpdivs, 'gpdivs','supf_elmope_2nd')
      call a%Memor%alloc(e%ndime,e%mnode,elrepdivs, 'elrepdivs','supf_elmope_2nd')  
      
   end subroutine

   subroutine DeallocOssSplit
      implicit none

      !Matrices alloc
      call a%Memor%dealloc(e%ndime,gpdivs, 'gpdivs','supf_elmope_2nd')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepdivs, 'elrepdivs','supf_elmope_2nd')    
      
   end subroutine
   
   subroutine GatherOssSplit
      implicit none
      
        call e%gather(e%ndime,elrepdivs,a%reproSDiv)
      
   end subroutine

   subroutine InterpolateOssSplit
      implicit none
      
      !Interpolate
      call e%interpg(e%ndime,elrepdivs,gpdivs)
      
   end subroutine   
   
   subroutine ResProSOSS
      implicit none
      
      call supf_elmbst_splitoss(e,dvol,timom,auxtens,beta,elmat)
      call supf_elmrhc_splitoss(e,auxtens,beta,timom,dvol,gpdivs,elrhs)
      
   end subroutine  
   
   subroutine ResProSOSS3d
      implicit none
      
      call supf_elmbst_splitoss3d(e,dvol,timom,auxtens,beta,elmat)
      call supf_elmrhc_splitoss3d(e,auxtens,beta,timom,dvol,gpdivs,elrhs)
      
   end subroutine 
   
  subroutine IngausSplitOss
   implicit none
   
      call ProcPointer%SplitOSSComponents
  
  end subroutine
   
   
   !--------------------------------------------------------------------------------
   ! Stabilization Parameters
   subroutine TauElastic
      use Mod_sup_ComputeTisig
      implicit none
      
      !Compute the stability parameters      
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)      
      !Viscoelastic Stabilization Constituive parameter 
      call sup_ComputeTiSig(e,lambda,acvis,gpvno2,grvel,a%staco,chale,auxG,a%PTT_model,gpsig,tisig)        
 
   
   end subroutine    

   !-------------------------------------------------------------------
   !Physical Properties
   subroutine ViscosityLaw 
      implicit none   

   end subroutine
   
   !------------------------------------------------------------------
   !Laplacian elemental subroutines
   subroutine InGaussElmatsNonLinear
      implicit none 
        
      
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
   
   !-------------------------------------------------------------------------------
   !extrapolation terms
   
   subroutine sup_ExtraExtrapoaltion2d
   implicit none   
  
   
   end subroutine
   
   subroutine sup_ExtraExtrapoaltion3d
   implicit none   
      
   end subroutine
   
   
   subroutine elmatGalerk
  
      !Galerkin
      call supf_elmbstGal(e,auxconv,auxtrac1,acvis,auxVE,auxG,auxPTT,LHSdtinv,dvol,auxtens,AGradV,grvel,gpsig,elmat)          
      call supf_elmrhcGal(e,auxVE,beta,auxG,auxPTT,acvis,dvol,auxtens,elextS,gpvel,grvel,gpsig,elrhs)

   
   end subroutine
   
   subroutine elmatGalerk3d
  
      !Galerkin
      call supf_elmbstGal3d(e,auxconv,auxtrac1,acvis,auxVE,auxG,auxPTT,LHSdtinv,dvol,auxtens,AGradV,grvel,gpsig,elmat)    
      call supf_elmrhcGal3d(e,auxVE,beta,auxG,auxPTT,acvis,dvol,auxtens,elextS,gpvel,grvel,gpsig,elrhs)
      
   end subroutine 
   
   subroutine elmatEstab

      !Constitutive
      call supf_elmbstEst(e,auxGie,auxconv,auxtrac1,auxtrac2,auxVE,auxG,auxPTT,acvis,tisig,LHSdtinv,AGradV,gpsig,dvol,auxtens,grvel,elmat) 
      call supf_elmrhcEst(e,auxGie,auxconv,auxtrac1,auxtrac2,auxVE,beta,auxG,auxPTT,tisig,acvis,gpsig,gpvel,AGradV,grvel,dvol,auxtens,elextS,elrhs) 

   end subroutine
   
   subroutine elmatEstab_Skip
   
      !Constitutive
      call supf_elmbstEst(e,auxGie,auxconv,auxtrac1,auxtrac2,auxVE,auxG,auxPTT,acvis,tisig,0.0_rp,AGradV,gpsig,dvol,auxtens,grvel,elmat) 
      call supf_elmrhcEst(e,auxGie,auxconv,auxtrac1,auxtrac2,auxVE,beta,auxG,auxPTT,tisig,acvis,gpsig,gpvel,AGradV,grvel,dvol,auxtens,elextS,elrhs) 
   
   
   end subroutine
   
   subroutine elmatEstab3d

      !Constitutive
      call supf_elmbstEst3d(e,auxconv,auxtrac1,auxtrac2,auxVE,auxG,auxPTT,acvis,tisig,LHSdtinv,AGradV,gpsig,dvol,auxtens,grvel,elmat) 
      call supf_elmrhcEst3d(e,auxconv,auxtrac1,auxtrac2,auxVE,beta,auxG,auxPTT,tisig,acvis,gpsig,gpvel,AGradV,grvel,dvol,auxtens,elextS,elrhs) 
 
   end subroutine  
   
   subroutine elmatEstab3d_Skip

      !Constitutive
      call supf_elmbstEst3d(e,auxconv,auxtrac1,auxtrac2,auxVE,auxG,auxPTT,acvis,tisig,0.0_rp,AGradV,gpsig,dvol,auxtens,grvel,elmat) 
      call supf_elmrhcEst3d(e,auxconv,auxtrac1,auxtrac2,auxVE,beta,auxG,auxPTT,tisig,acvis,gpsig,gpvel,AGradV,grvel,dvol,auxtens,elextS,elrhs) 
 
   end subroutine     
  
  
   
   !-------------------------------------------------------
   !LevelSet
   
   subroutine AllocLevelsetTF
      implicit none   
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_2nd')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,xloc,'xloc','supf_elmope_2nd') 
      call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_2nd')
      
   
   end subroutine
   
   subroutine DeallocLevelSetTF
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_2nd')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,xloc,'xloc','supf_elmope_2nd') 
      call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_2nd')      
      
   
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
      integer(ip) :: elemStatus,poinStatus
      
      
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
         wrmat1=0.0_rp
         
         
         ! Viscosity terms : we only consider mu*(grad v, grad u)         
         call elmvis(e,dvolt0,acvis,wrmat1)  

         forall (idime = 1:ntens)
            elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
         end forall          
         
      end if      
      
   end subroutine    
   
   subroutine AllocStokes
      implicit none
      integer(ip) :: ntens
      
      ntens = (e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(ntens,e%mnode,elrhs2,'elrhs2','supf_elmope_1st')      
   
   end subroutine
   
   subroutine DeallocStokes
      implicit none
      integer(ip) :: ntens
      
      ntens = (e%ndime-1)*(e%ndime-1)+2      
      
      call a%Memor%dealloc(ntens,e%mnode,elrhs2,'elrhs2','supf_elmope_1st')         
      
   end subroutine
   
   subroutine InitRHSStokes
      implicit none
      
      elrhs2 = 0.0_rp
      
   end subroutine
   
   subroutine RHSStokes
      implicit none
      integer(ip)  :: elemStatus
      
      call a%CutMesh%GetElementType(ielem,elemStatus)      
      
      if(elemStatus==-1)then
         
         elextS=0.0_rp
         
         if(e%ndime==2)then
            
            call supf_elmrhcGal(e,0.0_rp,beta,auxG,auxPTT,acvis,dvol,auxtens,elextS,gpvel,grvel,gpsig,elrhs2)         
         
         elseif(e%ndime==3)then
         
            call supf_elmrhcGal3d(e,0.0_rp,beta,auxG,auxPTT,acvis,dvol,auxtens,elextS,gpvel,grvel,gpsig,elrhs2)          
         
         end if
      end if
      
   end subroutine
   
   
   subroutine ElmatsToStokes
      implicit none 
      integer(ip)  :: elemStatus
      integer(ip) :: ntens
      
      ntens = (e%ndime-1)*(e%ndime-1)+2          
      
      call a%CutMesh%GetElementType(ielem,elemStatus)     
      
      if(elemStatus==-1)then
      
         elmat(:,:,:,:) = 0.0_rp     
         
         do inode=1,e%pnode
            do idime=1,ntens
               elmat(idime,inode,idime,inode) = elmat(idime,inode,idime,inode) + 1.0_rp/(2.0_rp*acvis)
            end do
         end do        
      
                 
         elrhs  = elrhs2       
         
      end if       
      
   end subroutine
   
   
   
        
  
  
  
  
end module   

subroutine supf_elmope_2nd(SUPFProblem)
   use Mod_SUPFractionalStep
   use Mod_supf_elmope_2nd
   implicit none
   class(SUPFractionalStepProblem), target :: SUPFProblem 

  
   
   a=>SUPFProblem
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   !Matrices Alloc
   call a%Memor%alloc(auxtens,e%mnode,auxtens,e%mnode,elmat,'elmat','supf_elmope_2nd')
   call a%Memor%alloc(auxtens,e%mnode,elrhs,'elrhs','supf_elmope_2nd')
   
   !Other arrays alloc
   call a%Memor%alloc(e%ndime,e%mnode,2,elvel,'elvel','supf_elmope_2nd')
   call a%Memor%alloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supf_elmope_2nd')
   call a%Memor%alloc(e%ndime,elext,'elext','supf_elmope_2nd')
   call a%Memor%alloc(auxtens,elextS,'elextS','supf_elmope_2nd') 
   call a%Memor%alloc(e%ndime,2,gpvel,'gpvel','supf_elmope_2nd')
   call a%Memor%alloc(auxtens,a%ncomp-1,gpsig,'gpsig','supf_elmope_2nd')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','supf_elmope_2nd')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','supf_elmope_2nd')
   call a%Memor%alloc(e%mnode,testf,'testf','supf_elmope_2nd')
   
   !gradients
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','supf_elmope_2nd') 
   call a%Memor%alloc(auxtens,e%ndime,grsig,'grsig','supm_elmope_2nd')  
   
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
      
      !constitutive parameters
      auxconv=1_ip
      auxtrac1=1_ip
      auxtrac2=1_ip
      auxGie = 1_ip
   
      !Gathers
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      call e%gather(auxtens,elsig(:,:,1),a%sigma(:,:,1))
      !viscoelastic case
      call e%gather(auxtens,elsig(:,:,2),a%sigma(:,:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(auxtens,elsig(:,:,itime),a%sigma(:,:,itime+1))  
      enddo

      
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
         
         !Interpolate         
         call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1)) !i-1 iteration
         call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2)) !j-1 time
         call e%interpg(auxtens,elsig(:,:,1),gpsig(:,1)) !i-1 iteration
         call e%interpg(auxtens,elsig(:,:,2),gpsig(:,2)) !i-1 iteration         
         do itime = 3,nsteps ! Time bdf2 and others
            call e%interpg(auxtens,elsig(:,:,itime),gpsig(:,itime))
         enddo
         
         !velocity and pressure gradient 
         call e%gradient(e%ndime,elvel,grvel)   
         !Stress gradient 
         call e%gradient(auxtens,elsig,grsig) 
         
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
         !Advection velocity norm used in alfa3 and kdisc
         call vecnor(gpvel(:,1),e%ndime,gpvno2,2)
         !Compute a·grad(V)
         call ComputeAGradV(e,gpvel(:,1),AGradV)
         
         !Taus
         call ProcPointer%TauConstitutive         
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         call nsm_ComputeTestf(e,acden,timom,AGradV,testf)
         !Hook
         call ProcHook%Testf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                ! w(gp)*detjm*tau1
                  
         !Compute Elext, Temporal Derivatives, Repro...
         elext =0.0_rp
         elextS=0.0_rp
         elextC=0.0_rp
         
         !Viscoelastic time integration
         call sup_TimeIntegrationToElext(e,auxtens,integrator,auxVE,a%dtinv,gpsig,elextS)
   
         !Compute vector of external forces
         call ProcPointer%ExternalForces        
         !----------------------------------------------------------------------
         ! Elemental matrix elmat and RHS         

         !Compute contributions to elmat and RHS
         !Galerkin
         call ProcPointer%ViscoGalerkinMatrix
         !Residual based stabalized terms
         call ProcPointer%ViscoEstabMatrix    
         
         !Terms involving the extrapolate pressure and stress to the RHS         
         call ProcPointer%ExtrapolationTerms
         
         !Hook
         call ProcHook%InGaussElmats
         
         !Statistics
         call a%InGaussStats(acden,acvis,gpvno,chale,timom)
         
      enddo gauss_points
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer%PostGaussElmats
      
      
      !Pre Assembly Modifications
      !Pointer
      call ProcPointer%PreAssembly          

      !Dirichlet Boundary Conditions
      call nsm_rotdir(a,e,auxtens,elmat,elrhs)
      currentbvess=1
      bcstar=0_ip
      
      if(a%kfl_exacs/=0.and.a%kfl_timei==1)then
         call sup_Exaelmdir(a,e,auxtens,auxtens,bcstar,currentbvess,elmat,elrhs)
      else
         call php_elmdir(a,e,auxtens,auxtens,bcstar,currentbvess,elmat,elrhs)
      end if     
      
      !Assembly
      call a%LinearSystemS%Assembly(e,elmat,elrhs)
      
   enddo elements
   
   call a%FinalizeStats
   
   !Hook
   call ProcHook%Finalizations
   
   !Matrices deAlloc
   call a%Memor%dealloc(auxtens,e%mnode,auxtens,e%mnode,elmat,'elmat','supf_elmope_2nd')
   call a%Memor%dealloc(auxtens,e%mnode,elrhs,'elrhs','supf_elmope_2nd')
   
   !Other arrays dealloc
   call a%Memor%dealloc(e%ndime,e%mnode,2,elvel,'elvel','supf_elmope_2nd')
   call a%Memor%dealloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supf_elmope_2nd')
   call a%Memor%dealloc(e%ndime,elext,'elext','supf_elmope_2nd')
   call a%Memor%dealloc(auxtens,elextS,'elextS','supf_elmope_2nd')   
   call a%Memor%dealloc(e%ndime,2,gpvel,'gpvel','supf_elmope_2nd')
   call a%Memor%dealloc(auxtens,a%ncomp-1,gpsig,'gpsig','supf_elmope_2nd')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','supf_elmope_2nd')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','supf_elmope_2nd')
   call a%Memor%dealloc(e%mnode,testf,'testf','supf_elmope_2nd')
   
   !gradients dealloc
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','supf_elmope_2nd')
   call a%Memor%dealloc(auxtens,e%ndime,grsig,'grsig','supm_elmope_2nd')

   !ElementDealloc
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','supf_elmope_2nd')
   
end subroutine
