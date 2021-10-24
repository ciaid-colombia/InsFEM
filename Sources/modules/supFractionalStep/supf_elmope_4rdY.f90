module Mod_supf_elmope_4rdY
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
      procedure(), NOPASS, pointer :: ComputeAdvectionVelocity  => NULL() 
      procedure(), NOPASS, pointer :: ResProRHS => NULL()
      procedure(), NOPASS, pointer :: ExtrapolationTerms => NULL()
      procedure(), NOPASS, pointer :: TauConstitutive  => NULL()
      procedure(), NOPASS, pointer :: ElmatUV => NULL()
      procedure(), NOPASS, pointer :: RHSUV => NULL()
      procedure(), NOPASS, pointer :: ConvectiveTerms => NULL()
      procedure(), NOPASS, pointer :: PreAssembly   => NULL()    

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
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv

   integer(ip)           :: ielty0 = 0 !Previous element type
      
   real(rp), allocatable :: elvel(:,:,:)
   real(rp), allocatable :: elsig(:,:,:)
   real(rp), allocatable :: elpre(:,:),elext(:)
   real(rp), allocatable :: elmuv(:,:,:,:)
   real(rp), allocatable :: wrmat1(:,:)
   real(rp), allocatable :: AGradV(:), testf(:) 
   real(rp), allocatable :: gpvel(:,:),gpadv(:),gpadv2(:)
   real(rp), allocatable :: gpsig(:,:)
   
   real(rp)             ::  val_aux   
      
   real(rp)    :: chale(2),timom,tidiv,tisig,dvolt0,dvolt1,dvolt2,dvolt3,gprhs(3)
   integer(ip) :: nmean,auxtens
   real(rp)    :: acden,acvis,beta,auxvis,lambda,auxG
   real(rp)    :: reyno
   real(rp)    :: gpvno,gpvno2,gppre(1),gppre2(1),divvel,dummr,auxp1,auxp2,auxp3,auxp4,gppredif(1)   
   !Residual Projections   
   real(rp), allocatable :: elrep(:,:),gprep(:)
   real(rp), allocatable :: grpre(:,:),grvel(:,:),grvel2(:,:)  
   !Split Oss
   real(rp), allocatable :: gpconv(:),elrepconv(:,:) 
   real(rp), allocatable :: gplapl(:),elreplapl(:,:)
   real(rp)              :: gpdiv(1)
   real(rp), allocatable :: elrepdiv(:,:)
   
   integer(ip) :: idime,itime,igaus,itest,auxntens,auxndim,auxGrad,kdime,ldime,aux12
   integer(ip) :: currentbvess,bcstar
   
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
      ProcPointer%ComputeAdvectionVelocity => ComputeAdvectionVelocity
      ProcPointer%ExtrapolationTerms  => sup_ExtraExtrapoaltion2d
      ProcPointer%ResProRHS       => ResProRHS2d  
      ProcPointer%ElmatUV         => elmat2d
      ProcPointer%RHSUV           => rhs2d
      ProcPointer%ConvectiveTerms => elmatConvective 
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
      !Viscoelastic Part
      ProcPointer%TauConstitutive => TauElastic 
      
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      
      !MaterialProperties is called Always
      call ConcatenateProcedures(ProcHook%PhysicalProp,MaterialProperties)      
      
      
      !Advection velocity
      if(e%ndime==3)then
         ProcPointer%ConvectiveTerms => elmatConvective3d
      endif      
      
      !Advection velocity
      if (a%kfl_advec == 0) then
         ProcPointer%ComputeAdvectionVelocity => NULLSUB
         ProcPointer%ConvectiveTerms => NULLSUB
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
         ProcPointer%ElmatUV         => elmat3d
         ProcPointer%RHSUV           => rhs3d       
      endif      
 
      
      !-----------------------------------------------------------
      !ResidualProjection
      !OSS
      if (a%kfl_repro == 1 .or. a%kfl_repro==0) then
         call runend('ASGS and OSS in Viscoelastic fractional dont work yet')
       endif
      
      !-----------------------------------------------------------
      
      !Split-Oss
      if(a%kfl_repro == 2 .or. a%kfl_repro == 3)then                   
         call ConcatenateProcedures(ProcHook%Initializations,AllocOssSplit)
         call ConcatenateProcedures(ProcHook%Gathers,GatherOssSplit)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateOssSplit)
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsOssSplit)
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
      
      
      !--------------------------------------------------------------------------------
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
            call ConcatenateProcedures(ProcPointer%PreAssembly,FreeSurfMatsTF)            
         end if 
         
      end if       
      
      
   end subroutine   
   
   !AdvectionVelocity
   subroutine ComputeAdvectionVelocity
      implicit none
      
      gpadv  = gpvel(:,1)
      gpadv2 = gpvel(:,2)
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
      acden  = a%MatProp(imat)%densi
      beta   = a%MatProp(imat)%LawViParam(2)
         
      !incremental scheme for Weissenberg number
      call a%IncrementalLambda(imat,a%MatProp(imat)%LawViParam(3), lambda)                
         
      auxG   = a%MatProp(imat)%LawViParam(4)/((1.0_rp-beta) + 0.00001_rp)      
      

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
      dvolt2=0.0_rp
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
      call a%Memor%alloc(a%ndofn,e%mnode,elrep,'elrep','supf_elmope_1st')
      call a%Memor%alloc(a%ndofn,gprep,'gprep','supf_elmope_1st')

   end subroutine

   subroutine DeallocRep
      implicit none

      call a%Memor%dealloc(a%ndofn,e%mnode,elrep,'elrep','supf_elmope_1st')
      call a%Memor%dealloc(a%ndofn,gprep,'gprep','supf_elmope_1st')

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
      call supf_elmrhu_oss(e,auxtens,tisig,dvol,gprep,beta,elrhs)      

   end subroutine  
  
   subroutine ResProRHS3d
      implicit none

      !Compute contributions to RHS : Block C
      call supf_elmrhu_oss3d(e,auxtens,tisig,dvol,gprep,beta,elrhs)        
   
   end subroutine   

   !------------------------------------------------------------------------------------------------------------
   !Split-OSS
   
   subroutine AllocOssSplit
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%ndime,gpconv, 'gpconv','supf_elmope_1st')
      call a%Memor%alloc(e%ndime,e%mnode,elrepconv, 'elrepconv','supf_elmope_1st')   

      call a%Memor%alloc(1,e%mnode,elrepdiv, 'elrepdiv','supf_elmope_1st')
      
   end subroutine

   subroutine DeallocOssSplit
      implicit none

      !Matrices alloc
      call a%Memor%dealloc(e%ndime,gpconv, 'gpconv','supf_elmope_1st')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepconv, 'elrepconv','supf_elmope_1st')     

      call a%Memor%dealloc(1,e%mnode,elrepdiv, 'elrepdiv','supf_elmope_1st')
      
   end subroutine
   
   subroutine GatherOssSplit
      implicit none
      
      call e%gather(e%ndime,elrepconv,a%reproUGradU)
      call e%gather(1_ip,elrepdiv,a%reproDivU)
      
   end subroutine

   subroutine InterpolateOssSplit
      implicit none
      
      !Interpolate
      call e%interpg(e%ndime,elrepconv,gpconv)
      call e%interpg(1_ip,elrepdiv,gpdiv)
      
   end subroutine   
   
   subroutine InGaussElmatsOssSplit
      implicit none      
      
      call supf_rhu_soss_step1(e,dvol,acden,tidiv,gpdiv,elrhs)
!       call supf_rhu_soss_conv(e,dvol,acden,timom,AGradV,gpconv,elrhs)
      
   end subroutine  

   !--------------------------------------------------------------------------------
   ! Stabilization Parameters
   subroutine TauElastic
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none
      
      !Compute the stability parameters
      
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)      
      call sup_ComputeTidivVE(e,timom,a%staco,chale,tidiv)
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
   
   !-------------------------------------------------------------------------------
   !extrapolation terms
   
   subroutine sup_ExtraExtrapoaltion2d
   implicit none
      
      !Compute contributions to RHS : Block U
      call supf_elmrhu_corr(e,dvol,acden,LHSdtinv,gpvel(:,2),gppre(1),gpsig(:,1),elrhs)
        
      !Second order fractional step
      !Extra terms from the extrapolation initial values
      call supf_elmrhu_corr_extra(e,dvol,gppre2(1),gpsig(:,2),elrhs)
             
   
   end subroutine 
 
   
   !matrices
   
   subroutine elmat2d   
  
     !Constitutive stabilization term
     call supf_elmbuv_estS(e,dvol,acvis,acden,LHSdtinv,tidiv,tisig,beta,grvel,elmat)
   
   end subroutine
   
   subroutine elmat3d
     
     !Constitutive stabilization term
     call supf_elmbuv_estS3d(e,dvol,acvis,acden,LHSdtinv,tidiv,tisig,beta,grvel,elmat)   
   
   
   end subroutine
   
   
   subroutine rhs2d
   
      !Compute contributions to RHS : Block U
      call supf_elmrhuY(e,acden,timom,tisig,beta,acvis,dvol,gpvel(:,1),gpvel(:,2),grvel,grvel2,AGradV,elrhs)   
   
   end subroutine
   
   subroutine rhs3d
   
      !Compute contributions to RHS : Block U
      call supf_elmrhuY3d(e,acden,timom,tisig,beta,acvis,dvol,gpvel(:,1),gpvel(:,2),grvel,grvel2,AGradV,elrhs)
   
   end subroutine
   
   
   subroutine elmatConvective
      
      !Compute contributions to elemental matrix : Block U,V
      call supf_elmbuv1(e,dvol,acvis,beta,acden,LHSdtinv,aGradV,timom,gpvel,grvel,elmat)   
      call supf_rhu_soss_conv(e,dvol,acden,timom,AGradV,gpconv,elrhs)
      !Compute contributions to RHS : Block U
      call supf_elmrhuYConv(e,acden,timom,tisig,beta,acvis,dvol,gpvel(:,1),gpvel(:,2),grvel,grvel2,AGradV,elrhs)        
   end subroutine
   
   
   subroutine elmatConvective3d
   
      !Compute contributions to elemental matrix : Block U,V
      call supf_elmbuv13d(e,dvol,acvis,beta,acden,LHSdtinv,aGradV,timom,gpvel,grvel,elmat)    
      call supf_rhu_soss_conv(e,dvol,acden,timom,AGradV,gpconv,elrhs) 
      !Compute contributions to RHS : Block U
      call supf_elmrhuYConv3d(e,acden,timom,tisig,beta,acvis,dvol,gpvel(:,1),gpvel(:,2),grvel,grvel2,AGradV,elrhs)
  
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
      integer(ip)   :: elemStatus
      
      call a%CutMesh%GetElementType(ielem,elemStatus)        
      if(elemStatus==-1)then
      
         elmat = 0.0_rp
         elrhs = 0.0_rp
         
         do inode=1,e%pnode
            do idime=1,e%ndime
               elmat(idime,inode,idime,inode) = elmat(idime,inode,idime,inode) + 1.0_rp
               
            end do
         end do           
         
         
         elrhs(1:e%ndime,:) = elvel(1:e%ndime,:,1) 
               
         
      end if    
      
   end subroutine   
   
   
   
   
  
end module   

subroutine supf_elmope_4rdY(SUPFProblem)
   use Mod_SUPFractionalStep
   use Mod_supf_elmope_4rdY
   implicit none
   class(SUPFractionalStepProblem), target :: SUPFProblem   
   
   a=>SUPFProblem
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_elmope_1st')     
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   !Matrices Alloc
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','supf_elmope_1st')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','supf_elmope_1st')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_1st')
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supf_elmope_1st')
   
   !Other arrays alloc
   call a%Memor%alloc(e%ndime,e%mnode,3,elvel,'elvel','supf_elmope_1st')
   call a%Memor%alloc(auxtens,e%mnode,2,elsig,'elsig','supf_elmope_1st')
   call a%Memor%alloc(      e%mnode,2,elpre,'elpre','supf_elmope_1st')
   call a%Memor%alloc(e%ndime,3,gpvel,'gpvel','supf_elmope_1st')
   call a%Memor%alloc(auxtens,2,gpsig,'gpsig','supf_elmope_1st')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','supf_elmope_1st')
   call a%Memor%alloc(e%ndime,gpadv2,'gpadv2','supf_elmope_1st') 
   call a%Memor%alloc(e%ndime,elext,'elext','supf_elmope_1st')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','supf_elmope_1st')
   call a%Memor%alloc(e%mnode,testf,'testf','supf_elmope_1st')
   
   !gradients
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','supf_elmope_1st')
   call a%Memor%alloc(e%ndime,e%ndime,grvel2,'grvel2','supf_elmope_1st')  
   
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
      wrmat1 = 0.0_rp
   
      !Gathers
      !iteration i-1
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      !intermediate values
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,2))
      call e%gather(auxtens,elsig(:,:,1),a%sigma(:,:,2))
      call e%gather(1,elpre(:,1),a%press(:,1))
      !time step n
      call e%gather(e%ndime,elvel(:,:,3),a%veloc(:,:,3))
      call e%gather(auxtens,elsig(:,:,2),a%sigma(:,:,3))
      call e%gather(1,elpre(:,2),a%press(:,3))
      
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
      dvolt2=0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm

         !Interpolate
         !iteration i-1
         call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1)) !i-1 iteration
         !intermediate values
         call e%interpg(e%ndime,elvel(:,:,2),gpvel(:,2)) !j-1 time
         call e%interpg(auxtens,elsig(:,:,1),gpsig(:,1)) !i-1 iteration
         call e%interpg(1,elpre(:,1),gppre(1))         
         !time step n 
         call e%interpg(e%ndime,elvel(:,:,3),gpvel(:,3)) !i-1 iteration         
         call e%interpg(auxtens,elsig(:,:,2),gpsig(:,2)) !i-1 iteration
         call e%interpg(1,elpre(:,2),gppre2(1))
         
         !velocity and pressure gradient 
         call e%gradient(e%ndime,elvel(:,:,1),grvel)
         
         !velocity and pressure gradient 
         call e%gradient(e%ndime,elvel(:,:,2),grvel2)         
     
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
         !Advection velocity norm used in alfa3
         call vecnor(gpvel(:,1),e%ndime,gpvno2,2)
         
         !Compute a·grad(V)
         call ComputeAGradV(e,gpadv,AGradV)
         
         call ProcPointer%TauConstitutive 
         
         elext =0.0_rp
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         call nsm_ComputeTestf(e,acden,timom,AGradV,testf)
         !Hook
         call ProcHook%Testf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                ! w(gp)*detjm*tau1
         dvolt2=dvol*tidiv + dvolt2                ! w(gp)*detjm*tau2  
         
         !----------------------------------------------------------------------
         ! Elemental matrix elmat and RHS         
         
         call ProcPointer%ElmatUV
         call ProcPointer%ConvectiveTerms           
         call ProcPointer%RHSUV

         !--------------------------------------------------------------------
         ! Second order fractional step
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
      
      
      !Matrix composition
      do idime = 1,e%ndime
         elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      enddo      
      
      !Pre Assembly Modifications
      !Pointer
      call ProcPointer%PreAssembly        
      
     
      !ACA TENGO DUDAS CON LAS CONDICIONES DE BORDE
      !Dirichlet Boundary Conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
      !inside the subroutine
      ! php_elmdir(a,e,ndofn,ndofbc,ndofbcstart,currentbvess,elmat,elrhs)
      currentbvess=auxtens+1
      bcstar=0_ip
      
      if(a%kfl_confi==1.and.a%MatProp(imat)%lawvi<0)then
         call sup_Exaelmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)
      elseif(a%kfl_confi==1.and.a%kfl_bc_number>0)then
         call sup_Exaelmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)         
      else
         call php_elmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)
      end if        
     
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo elements
   
   call a%FinalizeStats
   
   !Hook
   call ProcHook%Finalizations
   
   !Matrices deAlloc
   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','supf_elmope_1st')
   call a%Memor%dealloc(e%ndime,e%mnode,elrhs,'elrhs','supf_elmope_1st')
   call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_1st')
   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supf_elmope_1st')
   
   !Other arrays dealloc
   call a%Memor%dealloc(e%ndime,e%mnode,3,elvel,'elvel','supf_elmope_1st')
   call a%Memor%dealloc(auxtens,e%mnode,2,elsig,'elsig','supf_elmope_1st')
   call a%Memor%dealloc(      e%mnode,2,elpre,'elpre','supf_elmope_1st')  
   call a%Memor%dealloc(e%ndime,3,gpvel,'gpvel','supf_elmope_1st')
   call a%Memor%dealloc(auxtens,2,gpsig,'gpsig','supf_elmope_1st')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','supf_elmope_1st')
   call a%Memor%dealloc(e%ndime,gpadv2,'gpadv2','supf_elmope_1st')
   call a%Memor%dealloc(e%ndime,elext,'elext','supf_elmope_1st')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','supf_elmope_1st')
   call a%Memor%dealloc(e%mnode,testf,'testf','supf_elmope_1st')
   
   !gradients dealloc
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','supf_elmope_1st')        
   call a%Memor%dealloc(e%ndime,e%ndime,grvel2,'grvel2','supf_elmope_1st')  

   !ElementDealloc
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','supf_elmope_1st')
   
end subroutine
