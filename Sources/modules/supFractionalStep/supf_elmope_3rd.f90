module Mod_supf_elmope_3rd
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
      procedure(), NOPASS, pointer :: ExternalForces => NULL()
      procedure(), NOPASS, pointer :: ComputeAdvectionVelocity   => NULL()
      procedure(), NOPASS, pointer :: ResProRHS  => NULL()
      procedure(), NOPASS, pointer :: ExtrapolationTerms => NULL()
      procedure(), NOPASS, pointer :: TauConstitutive => NULL()
      procedure(), NOPASS, pointer :: PenaltyTerm => NULL()
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
   real(rp)              :: dvol,dvolt02

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv

   integer(ip)           :: ielty0 = 0 !Previous element type
      
   real(rp), allocatable :: elvel(:,:,:)
   real(rp), allocatable :: elsig(:,:,:)
   real(rp), allocatable :: elpre(:,:)
   real(rp), allocatable :: elext(:),elextS(:),elext2(:)
   real(rp), allocatable :: AGradV(:), testf(:) 
   real(rp), allocatable :: gpvel(:,:),gpadv(:)
   real(rp), allocatable :: gpsig(:,:)
   real(rp), allocatable :: divsig1(:),divsig2(:)   
   
      
   real(rp)    :: chale(2),timom,dvolt0,dvolt1,gprhs(3),AverageTimom
   integer(ip) :: nmean,auxtens,jpoin,jnode
   real(rp)    :: acden,acvis,elextC(1)
   real(rp)    :: reyno
   real(rp)    :: gpvno,gppre(1,1),divvel,dummr,Residuo(3)   
   !Gradients   
   real(rp), allocatable :: grpre(:,:),grvel(:,:)  
   !Split Oss
   real(rp), allocatable :: gpgrap(:),elrepgrap(:,:)     
   
   integer(ip) :: idime,itime,igaus,itest,auxntens,auxndim,auxGrad,kdime,ldime
   
   !Level Set
   integer(ip)              :: inode,ipoin
   integer(ip)              :: elemStatus,ngauss_minus,ngauss_plus,ngaus_total,poinStatus
   real(rp), allocatable    :: weigp(:)
   real(rp), pointer        :: xloc(:,:) => NULL()
   !Laplacian Term
   real(rp), allocatable :: wrmat1(:,:)    
   !To define the pressure boundary condition over free surface
   real(rp), allocatable        :: points(:,:)
   !Nodal coordinates
   real(rp), pointer       :: coord(:)  => NULL()  
   real(rp)                :: MaxResiduo(3)
      
   
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
      ProcPointer%ExtrapolationTerms  => sup_ExtraExtrapolation2d
      !Penalty
      ProcPointer%PenaltyTerm  => NULLSUB      
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
      
      
      !First order Fractional step method
      if(a%kfl_tsche_1st_datafile == 'BDF1 ')then 
         ProcPointer%ExtrapolationTerms  => NULLSUB   
      end if
      
      !MaterialProperties is called Always
      call ConcatenateProcedures(ProcHook%PhysicalProp,MaterialProperties)         
      
      !Advection velocity
      if (a%kfl_advec == 0) then
         ProcPointer%ComputeAdvectionVelocity => NULLSUB
      endif
      
      !Penalty term
      if(a%kfl_penal ==1)then
         ProcPointer%PenaltyTerm => elmatpenalty
      end if      
            
      !-----------------------------------------------------------
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
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
            !Always after FreeSurfMatsTF
            call ConcatenateProcedures(ProcPointer%PreAssembly,PressureZeroFurf)
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

!       ! laplacin term: (grad p, grad q)
!       call elmvis(e,dvolt0,1.0_rp,elmat)     

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
   
   !-------------------------------------------------------------------
   !Physical Properties
   subroutine MaterialProperties 
      implicit none 
      !----------------------------------------------------------------
      !Incremental strategy (Continuation method in terms of lambda)
      acvis  = a%MatProp(imat)%LawViParam(1) 
      acden  = a%MatProp(imat)%densi      

   end subroutine     

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
   !------------------------------------------------------------------------------------------------------------------------------------------------------
   !Penalty term       
   subroutine elmatpenalty
      
      call supm_elmbpqpena(e,acvis,a%penal,dvol,elmat) 
   
   end subroutine   

   !------------------------------------------------------------------------------------------------------------
   !Split-OSS
   
   subroutine AllocOssSplit
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%ndime,gpgrap, 'gpgrap','supf_elmope_3rd')
      call a%Memor%alloc(e%ndime,e%mnode,elrepgrap, 'elrepgrap','supf_elmope_3rd')   
      
   end subroutine

   subroutine DeallocOssSplit
      implicit none

      !Matrices alloc
      call a%Memor%dealloc(e%ndime,gpgrap, 'gpgrap','supf_elmope_3rd')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepgrap, 'elrepgrap','supf_elmope_3rd')     
      
   end subroutine
   
   subroutine GatherOssSplit
      implicit none
      
      call e%gather(e%ndime,elrepgrap,a%reproGradp)
      
   end subroutine

   subroutine InterpolateOssSplit
      implicit none
      
      !Interpolate
      call e%interpg(e%ndime,elrepgrap,gpgrap)
      
   end subroutine   
   
   subroutine InGaussElmatsOssSplit
      implicit none
      
      call supf_elmrhp_splitoss(e,timom,dvol,gpgrap,elrhs)
      
   end subroutine  
  
   !--------------------------------------------------------------------------------
   ! Stabilization Parameters
   subroutine TauElastic
      implicit none
      integer(ip)  :: elemStatus 
      
      !Compute the stability parameters      
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)        
      
      
!       if(a%kfl_colev==1)then
!          timom=0.0_rp
!       end if

      
!       if(a%kfl_colev==1)then
!       
!          call a%CutMesh%GetElementType(ielem,elemStatus)      
!          
!          if(elemStatus==0)then
!             timom=0.0_rp
!          end if
!          
!       end if
      
   
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
   
   subroutine sup_ExtraExtrapolation2d
   implicit none
   
     call supf_elmrhpextra(e,dvol,acden,LHSDtinv,divsig2,grpre,elrhs)
     
   end subroutine 
   
   !-------------------------------------------------------
   !LevelSet
   
   subroutine AllocLevelsetTF
      implicit none   
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_3rd')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,points,'points','supf_elmope_3rd')     
      call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_3rd')
      

   
   end subroutine
   
   subroutine DeallocLevelSetTF
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_3rd')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,points,'points','supf_elmope_3rd')   
      call a%Memor%dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_3rd')
      
      
   
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
   
   
   subroutine PressureZeroFurf
      implicit none
      integer(ip) :: elemStatus,poinStatus,JpoinStatus,ndone,ntodo
      integer(ip) :: nipoin,ipoin,jpoin,jnode
      real(rp)    :: pdsurf,toleralce
      real(rp)    :: diference(3),hipo,tmp
      real(rp)    :: localx1(3),localx2(3),centro(3),checkvector(3)
      !3d case
      real(rp)    :: diference2(3),hipo2,localx3(3),hipo3,projCheck
      integer(ip) :: ialone !node alone
      real(rp)    :: Residuo(3)

      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      toleralce=1e-10
      
      if(elemStatus==0)then  
      
         !To integrate over the surface or line intersection
         call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
         call a%CutMesh%GetRealIntersectionPoints(ielem,e,points)
         call a%CutMesh%GetSurfaceIntersection(ielem,e,pdsurf)
         call a%CutMesh%GetNinters(ielem,nipoin)
         
               
         !Compute shape functions at the intersection points
         weigp(1:nipoin) = 1.0_rp/nipoin
         
            
         !The rutine give the needed shape functions associated 
         call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)    
               
         
         !Initializations
         diference=0.0_rp
         hipo=0.0_rp
         diference2=0.0_rp
         hipo2=0.0_rp
         hipo3=0.0_rp
         localx1=0.0_rp
         localx2=0.0_rp
         localx3=0.0_rp
         centro=0.0_rp
         checkvector=0.0_rp
         projCheck=0.0_rp
         !using as a reference points(:,1)
         
         do idime=1,e%ndime
            diference(idime) = diference(idime) + (points(idime,1) - points(idime,2))
            hipo = hipo + diference(idime)*diference(idime)
         end do
         
         hipo=sqrt(hipo)
         !known direction
         if(hipo>toleralce)then
            localx1   = diference/hipo   
         elseif(hipo<toleralce)then
            localx1 = 0.0_rp
         end if
         
         
         if(e%ndime==2)then
            !in 2d the orthogonal vector is directly
            localx2(1)=-localx1(2)
            localx2(2)=localx1(1)
         elseif(e%ndime==3)then   
         
            do idime=1,e%ndime
               diference2(idime) = diference2(idime) + (points(idime,1) - points(idime,3))
               hipo2 = hipo2 + diference2(idime)*diference2(idime)
            end do           
            hipo2=sqrt(hipo2)
            
            if(hipo2>toleralce)then
               localx3=diference2/hipo2
            elseif(hipo2<toleralce)then
               localx3=0.0_rp
            endif
            
            !We use the classical vectorial product rule
            localx2(1)= localx1(2)*localx3(3)-localx1(3)*localx3(2)
            localx2(2)= -(localx1(1)*localx3(3)-localx1(3)*localx3(1))
            localx2(3)= localx1(1)*localx3(2)-localx1(2)*localx3(1)
            
            do idime=1,e%ndime
               hipo3=hipo3+localx2(idime)*localx2(idime)               
            end do
            hipo3=sqrt(hipo3)
            
            if(hipo3>toleralce)then
               localx2=localx2/hipo3
            elseif(hipo3<toleralce)then            
               localx2=0.0_rp
            end if
         end if  
         
         !we need check if the normal is outer or inner. We need the outer  
         
         do ipoin=1,nipoin
            do idime=1,e%ndime
               centro(idime)=centro(idime) + points(idime,ipoin)
            end do
         end do
         !Geometrical gravity center 
         centro=centro/nipoin
         
         call a%CutMesh%GetIalone(ielem,ialone)
         if(ialone/=0)then
            ipoin=e%lnods(ialone)   
            call a%Mesh%GetPointCoord(ipoin,coord)      
            
            do idime=1,e%ndime
               checkvector(idime)= checkvector(idime) + (coord(idime)-centro(idime))
            end do
            
            projCheck   = dot_product(checkvector,localx2)
            
            call a%CutMesh%GetPointType(ipoin,poinStatus)
            
            if(projCheck>0.0_rp .and. poinStatus>0_ip) localx2=-localx2
            if(projCheck<0.0_rp .and. poinStatus<=0_ip) localx2=-localx2
         elseif(ialone==0)then 
         
            do inode=1,e%pnode
               jpoin=e%lnods(inode)
               if(poinStatus<=0.0_rp) ipoin=jpoin
            end do
            
            call a%Mesh%GetPointCoord(ipoin,coord)      
            
            do idime=1,e%ndime
               checkvector(idime)= checkvector(idime) + (coord(idime)-centro(idime))
            end do
            
            projCheck   = dot_product(checkvector,localx2)
            
            if(projCheck<0.0_rp) localx2=-localx2         
         
         end if    
         
         if(pdsurf>toleralce)then
         
            !Approximate imposition of boundary conditions in immersed boundary methods
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               if(poinStatus==-1)then  
                  elmat(:,inode,:,1:e%pnode) = 0.0_rp
                  elrhs(:,inode) = 0.0_rp
                  do jnode = 1,e%pnode
                     do igaus = 1,e%pgaus
                        e%igaus = igaus                    
                        tmp=1.0_rp/(a%MatProp(1)%densi*LHSdtinv)
      
                        elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
                        - e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*(tmp/chale(1))*weigp(e%igaus)*pdsurf

                     enddo  
                  enddo                
               
               end if         
            end do           
            
            !terms from the integrations by parts
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               if(poinStatus==1)then
                  do jnode = 1,e%pnode
                     do idime=1,e%ndime                  
                        do igaus = 1,e%pgaus
                           e%igaus = igaus
                           tmp=1.0_rp/(a%MatProp(1)%densi*LHSdtinv)
                           
                           elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
                              + e%shape(inode,e%igaus)*e%cartd(idime,jnode)*localx2(idime)*(tmp)*(weigp(e%igaus)*pdsurf)     
               
                           
                        end do
                     enddo                  
                  enddo                
               
               end if         
            end do  
            
            !terms from the integrations by parts
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               if(poinStatus==1)then
                  do idime=1,e%ndime                  
                     do igaus = 1,e%pgaus
                        e%igaus = igaus
                        tmp=1.0_rp/(a%MatProp(1)%densi*LHSdtinv)               

                              
                        !First order Fractional step method
                        if(a%kfl_tsche_1st_datafile /= 'BDF1 ')then   
                           
                           elrhs(1,inode) = elrhs(1,inode) &
                              + e%shape(inode,e%igaus)*grpre(1,idime)*tmp*localx2(idime)*(weigp(e%igaus)*pdsurf)                             
                              
                        end if
                              
                           
                     end do
                  enddo                  
               
               end if         
            end do           
         
         elseif(pdsurf<=toleralce)then
         
            !The node is very close to the surface so the node is fixed to zero value
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               if(poinStatus==-1)then  
                  elmat(:,inode,:,1:e%pnode) = 0.0_rp
                  elrhs(:,inode) = 0.0_rp
                     
                  elmat(1,inode,1,inode) = elmat(1,inode,1,inode) + 1.0_rp
                  elrhs(:,inode) = 0.0_rp

               end if         
            end do              
         
         end if
         
!          Residuo=0.0_rp
!          do inode=1,e%pnode
!             do jnode=1,e%pnode
!             jpoin=e%lnods(jnode)
!             
!             Residuo(inode)=Residuo(inode)+(elmat(1,inode,1,jnode)*a%press(jpoin,1))
!             
!             end do
!          end do
!          
!          Residuo(1:e%pnode)=Residuo(1:e%pnode)-elrhs(1,1:e%pnode)         
      
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
         
!          ! Viscosity terms : we only consider mu*(grad v, grad u)         
!          call elmvis(e,dvolt0,acvis,wrmat1)  
! 
!          elmat(1,1:e%pnode,1,1:e%pnode) = elmat(1,1:e%pnode,1,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)    
         
         do inode=1,e%pnode         
            elmat(1,inode,1,inode) =  1.0_rp
         end do
         
         
      end if      
      
   end subroutine   
   
   subroutine AllocStokes
      implicit none
      
      call a%Memor%alloc(1,e%mnode,elrhs2,'elrhs2','supf_elmope_3rd')      
   
   end subroutine
   
   subroutine DeallocStokes
      implicit none
      
      call a%Memor%dealloc(1,e%mnode,elrhs2,'elrhs2','supf_elmope_3rd')         
      
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
         divsig1 = 0.0_rp
         
         call supf_elmrhp(e,dvol,acden,timom,LHSDtinv,divvel,divsig1,elrhs2)
      end if   
   end subroutine
   
   
   subroutine ElmatsToStokes
      implicit none 
      integer(ip)  :: elemStatus   
      
      call a%CutMesh%GetElementType(ielem,elemStatus)        
      if(elemStatus==-1)then
      
         elmat(:,:,:,:) = 0.0_rp         
         call nsm_elmbpq(e,dvolt1,elmat)
                 
         elrhs  = elrhs2       
         
      end if       
      
   end subroutine
   
   
   
   
   
  
end module   

subroutine supf_elmope_3rd(SUPFProblem)
   use Mod_SUPFractionalStep
   use Mod_supf_elmope_3rd
   implicit none
   class(SUPFractionalStepProblem), target :: SUPFProblem   
   
   a=>SUPFProblem
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_elmope_3rd')     
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   !Matrices Alloc
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmat,'elmat','supf_elmope_3rd')
   call a%Memor%alloc(1,e%mnode,elrhs,'elrhs','supf_elmope_3rd')
   
   !Other arrays alloc
   call a%Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','supf_elmope_3rd')
   call a%Memor%alloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supf_elmope_3rd')
   call a%Memor%alloc(      e%mnode,a%ncomp-1,elpre,'elpre','supf_elmope_3rd')
   call a%Memor%alloc(e%ndime,elext,'elext','supf_elmope_3rd')
   call a%Memor%alloc(e%ndime,elext2,'elext2','supf_elmope_3rd')   
   call a%Memor%alloc(auxtens,elextS,'elextS','supf_elmope_3rd') 
   call a%Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supf_elmope_3rd')
   call a%Memor%alloc(auxtens,a%ncomp-1,gpsig,'gpsig','supf_elmope_3rd')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','supf_elmope_3rd')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','supf_elmope_3rd')
   call a%Memor%alloc(e%mnode,testf,'testf','supf_elmope_3rd')
   
   !gradients
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','supf_elmope_3rd')        
   call a%Memor%alloc(1_ip     ,e%ndime,grpre,'grpre','supf_elmope_3rd')
   !Divsig
   call a%Memor%alloc(e%ndime,divsig1,'divsig1','supf_elmope_3rd')
   call a%Memor%alloc(e%ndime,divsig2,'divsig2','supf_elmope_3rd') 
   
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
      AverageTimom=0.0_rp
      
 
   
      !Gathers
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1))
      call e%gather(e%ndime,elvel(:,:,2),a%veloc(:,:,3))
      call e%gather(auxtens,elsig(:,:,1),a%sigma(:,:,1))
      !viscoelastic case
      call e%gather(auxtens,elsig(:,:,2),a%sigma(:,:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(e%ndime,elvel(:,:,itime),a%veloc(:,:,itime+1))
         call e%gather(auxtens,elsig(:,:,itime),a%sigma(:,:,itime+1))  
      enddo
      !Pressure Gather
      call e%gather(1_ip   ,elpre(1:e%pnode,1),a%press(1:e%pnode,1))
      
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
      dvolt02=0.0_rp
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
            call e%interpg(e%ndime,elvel(:,:,itime),gpvel(:,itime))
            call e%interpg(auxtens,elsig(:,:,itime),gpsig(:,itime))
         enddo
         
         call e%interpg(1,elpre(:,1),gppre)
         !velocity and pressure gradient 
         call e%gradient(1,elpre(:,1),grpre)
         call e%gradient(e%ndime,elvel,grvel)         
         
         call e%divergence(elvel(:,:,1),divvel)    
         
!          if(a%istep==1) write(*,*) 'divvel',divvel,ielem
!          if(a%istep==1) write(*,*) 'grpre',grpre(1,1),grpre(1,2),ielem
!          if(a%istep==1) write(*,*) 'gpgrap',gpgrap(1),gpgrap(2),ielem
         
         
        !Hook
         call supf_stressdivergence(e,auxtens,elsig(:,:,2),divsig2)
         call supf_stressdivergence(e,auxtens,elsig(:,:,1),divsig1)
         
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
         
         !----------------------------------------------------------------------
         ! Elemental matrix elmat and RHS        
         dvolt02 = dvolt02  + (1.0_rp/(acden*LHSdtinv))*dvol  
         !Compute contributions to RHS : Block U
         !InGaussElmats         
         
         call supf_embpq_step3(e,dvol,acden,timom,LHSDtinv,elmat)
         
         call supf_elmrhp(e,dvol,acden,timom,LHSDtinv,divvel,divsig1,elrhs)
         
         !only in splitt OSS with gravity         
         elext2=0.0_rp      
         !Compute vector of external forces
         call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext2)  
         call supm_elmrhpForce(e,0_ip,(-1.0_rp*timom),dvol,elext2,elrhs)         
         
         !-----------------------------------------------------------------------
         !Second order Fractional error     
         !Terms involving the extrapolate pressure and stress to the RHS         
         call ProcPointer%ExtrapolationTerms
         
         call ProcPointer%PenaltyTerm
         
         !Hook
         call ProcHook%InGaussElmats
         
         !Statistics
         call a%InGaussStats(acden,acvis,gpvno,chale,timom)
         
      enddo gauss_points
      
      
      AverageTimom=dvolt1/dvolt0
      
      
      !Post Gauss Matrices
      !Pointer
      call ProcPointer%PostGaussElmats
      
      !Pre Assembly Modifications
      !Pointer
      call ProcPointer%PreAssembly  
      
     
      !Dirichlet Boundary conditions
      if(a%kfl_exacs/=0.and.a%kfl_timei==1)then
         call supf_ExaPreElmdir(a,e,elmat,elrhs)
      else
         call supf_elmdir_press(a,e,elmat,elrhs)
      end if
    

!       if(e%ndime==2)then
!          Residuo=0.0_rp
!          do inode=1,e%pnode
!             do jnode=1,e%pnode
!             jpoin=e%lnods(jnode)
!             
!             Residuo(inode)=Residuo(inode)+(elmat(1,inode,1,jnode)*a%press(jpoin,1))
!             
!             end do
!          end do
!          
!          Residuo(1:e%pnode)=Residuo(1:e%pnode)-elrhs(1,1:e%pnode)         
!            
!       end if  
     
!       if(e%ndime==2)then
!          do inode=1,e%pnode
!             if(Residuo(inode)>1e-14) write(*,*) 'Residuo(inode)',ielem,Residuo(inode)
!          end do
!       end if
      
      !Assembly
      call a%LinearSystemC%Assembly(e,elmat,elrhs)
      
   enddo elements
   
   call a%FinalizeStats
   
   !Hook
   call ProcHook%Finalizations
   
   !Matrices deAlloc
   call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmat,'elmat','supf_elmope_3rd')
   call a%Memor%dealloc(1,e%mnode,elrhs,'elrhs','supf_elmope_3rd')
   
   !Other arrays dealloc
   call a%Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elvel,'elvel','supf_elmope_3rd')
   call a%Memor%dealloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supf_elmope_3rd')
   call a%Memor%dealloc(      e%mnode,a%ncomp-1,elpre,'elpre','supf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,elext,'elext','supf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,elext2,'elext2','supf_elmope_3rd')   
   call a%Memor%dealloc(auxtens,elextS,'elextS','supf_elmope_3rd')   
   call a%Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supf_elmope_3rd')
   call a%Memor%dealloc(auxtens,a%ncomp-1,gpsig,'gpsig','supf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','supf_elmope_3rd')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','supf_elmope_3rd')
   call a%Memor%dealloc(e%mnode,testf,'testf','supf_elmope_3rd')
   
   !gradients dealloc
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','supf_elmope_3rd')        
   call a%Memor%dealloc(1_ip     ,e%ndime,grpre,'grpre','supf_elmope_3rd') 
   
   !Divsig
   call a%Memor%dealloc(e%ndime,divsig1,'divsig1','supf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,divsig2,'divsig2','supf_elmope_3rd')   

   !ElementDealloc
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','supf_elmope_3rd')
   
end subroutine
