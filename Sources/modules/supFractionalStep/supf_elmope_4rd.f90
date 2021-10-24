Module Mod_supf_elmope_4rd
   use typre
   use Mod_SUPFractionalStep
   use Mod_TimeIntegrator
   use Mod_Memor
   use Mod_Element
   use Mod_SUPF_Element   
   use Mod_ConvectiveElement
   use Mod_php_SetTimeIntegrator
   use Mod_sup_elmdir
   use Mod_nsm_elmdir
   use Mod_php_elmdir

   implicit none   
   
   class(SUPFractionalStepProblem), pointer :: a   
   real(rp), parameter :: zensi = 0.0_rp   
   
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: PostGaussElmats => NULL()
      procedure(), NOPASS, pointer :: ComputeAdvectionVelocity  => NULL() 
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
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()      
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()
    end type
   type(PHook) :: ProcHook       
      
   class(FiniteElement), pointer :: e => NULL()

   real(rp), allocatable    :: elmat(:,:,:,:)
   real(rp), allocatable    :: elrhs(:,:)
   real(rp), allocatable    :: wrmat1(:,:),elmuv(:,:,:,:)
   real(rp)                 :: acden,acvis,beta

   real(rp), allocatable    :: elvel(:,:,:),elpre(:,:)
   real(rp), allocatable    :: elsig(:,:,:),gpsig(:,:)
   real(rp), allocatable    :: elpreInc(:,:)  !Pressure Increment
   real(rp), allocatable    :: gpvel(:)
   real(rp)                 :: gppre1(1),gppre2(1)
   
   type(TimeIntegratorDt1):: Integrator
   real(rp)             :: LHSdtinv
   integer(ip)          :: nsteps,auxtens


   integer(ip)          :: ielty0
   integer(ip)          :: idime,igaus,nelem,ielem
   
   real(rp)             :: dvol,dvolt0
   real(rp)             :: val_aux
   integer(ip)          :: currentbvess,bcstar
   
   !Level Set
   integer(ip)              :: inode,ipoin
   integer(ip)              :: elemStatus,ngauss_minus,ngauss_plus,ngaus_total,poinStatus
   real(rp), allocatable    :: weigp(:),xloc(:,:)    
   
   !Todo Multy Materials
   integer(ip) :: imat=1 

   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"

   !----------------------------------------------------------  

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
      ProcPointer%PostGaussElmats => NULLSUB
      ProcPointer%ComputeAdvectionVelocity => NULLSUB     
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
      
      !MaterialProperties is called Always
      call ConcatenateProcedures(ProcHook%PhysicalProp,MaterialProperties)
      
      
      call ConcatenateProcedures(ProcHook%InGaussElmats,IGElmats)
      

      !--------------------------------------------------------------------
      !Level set
      if(a%kfl_colev==1)then
         call ConcatenateProcedures(ProcHook%Initializations,AllocLevelsetTF) 
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocLevelSetTF)
         call PrependProcedure(ProcHook%PreGauss,CutelementsTF)
         call PrependProcedure(ProcHook%PhysicalProp,FluidPropertiesTF) 
         if(a%kfl_fsurf==1)then
            if(a%kfl_fsurfLapla==1 .or. a%kfl_fsurfLapla==2)then
               call PrependProcedure(ProcPointer%PreAssembly,ElmatsToStokes)
            end if
            
            call ConcatenateProcedures(ProcPointer%PreAssembly,FreeSurfMatsTF)            
         end if 
         
      end if
      
      
      
   
   end subroutine
   
   !-------------------------------------------------------------------
   !Physical Properties
   subroutine MaterialProperties 
      implicit none 
      !----------------------------------------------------------------
      !Incremental strategy (Continuation method in terms of lambda)
      acvis  = a%MatProp(imat)%LawViParam(1) 
      acden  = a%MatProp(imat)%densi
      beta   = a%MatProp(imat)%Lawviparam(2)       

   end subroutine   
   
   subroutine IGElmats
      implicit none
      
      !Compute contributions to RHS : Block U
      call supf_elmrhu_corr(e,dvol,acden,LHSdtinv,gpvel,gppre1(1),gpsig(:,1),elrhs)
         
      !-------------------------------------------------------------------------------
      !Second order fractional step
      !Extra terms from the extrapolation initial values
      if(a%kfl_tsche_1st_datafile /= 'BDF1 ')then
         call supf_elmrhu_corr_extra(e,dvol,gppre2(1),gpsig(:,2),elrhs)
      end if      
      
   end subroutine
   
   
   
   !-------------------------------------------------------
   !LevelSet
   
   subroutine AllocLevelsetTF
      implicit none   
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_4rd')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,xloc,'xloc','supf_elmope_4rd') 
   
   end subroutine
   
   subroutine DeallocLevelSetTF
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','supf_elmope_4rd')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,xloc,'xloc','supf_elmope_4rd') 
      
   
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
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin,idime
      
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
   
   
   subroutine ElmatsToStokes
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

subroutine supf_elmope_4rd(SUPFProblem)
   use Mod_SUPFractionalStep
   use Mod_supf_elmope_4rd
   implicit none
   class(SUPFractionalStepProblem), target :: SUPFProblem   
   
   a=>SUPFProblem
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_elmope_4rd')
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers   
  
   !Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)

   !Memory Allocations
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','supf_elmope_4rd')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_4rd')
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supf_elmope_4rd')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','supf_elmope_4rd')
   
   call a%Memor%alloc(e%ndime,e%mnode,1,elvel,'elvel','supf_elmope_4rd')
   call a%Memor%alloc(e%mnode,2,elpre,'elpre','supf_elmope_4rd')
   call a%Memor%alloc(e%mnode,1,elpreInc,'elpreInc','supf_elmope_4rd')
   call a%Memor%alloc(e%ndime,gpvel,'gpvel','supf_elmope_4rd')
   
   auxtens=(e%ndime-1)*(e%ndime-1) + 2
   call a%Memor%alloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supf_elmope_4rd')      
   call a%Memor%alloc(auxtens,a%ncomp-1,gpsig,'gpsig','supf_elmope_4rd')      
   
   
   !Hook
   call ProcHook%Initializations   

   !Loop over elements
   call a%Mesh%GetNelem(nelem)
   elements: do ielem=1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
      
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      wrmat1=0.0_rp

      !Gathering operations
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1)) ! int. u_n+1
      call e%gather(1,elpre(:,1),a%press(:,1)) ! p_n+1
      call e%gather(1,elpre(:,2),a%press(:,3)) ! p_n
      call e%gather(auxtens,elsig(:,:,1),a%sigma(:,:,1))
      call e%gather(auxtens,elsig(:,:,2),a%sigma(:,:,3))
      
      !Hook
      call ProcHook%Gathers      
      
      !Cartesian derivatives and Jacobian for linear elements
      call e%elmdcg
      
      !DvolsToZero
      dvolt0=0.0_rp      

      !Loop on gauss points
      gauss_points: do igaus=1,e%pgaus
         e%igaus = igaus
         if (e%linea == 0) call e%elmder
         
         !Hook
         call ProcHook%InGauss         
         
         dvol = e%weigp(e%igaus)*e%detjm

         ! Interpolation (Gauss point values)
         call e%interpg(e%ndime,elvel(:,:,1),gpvel) ! u_int
         call e%interpg(1,elpre(:,1),gppre1)
         call e%interpg(1,elpre(:,2),gppre2)         
         
         call e%interpg(auxtens,elsig(:,:,1),gpsig(:,1)) !i-1 iteration
         call e%interpg(auxtens,elsig(:,:,2),gpsig(:,2)) !i-1 time  
         
        !Hook
         call ProcHook%Interpolates
         
         !Physical Parameters
         call a%GetPhysicalParameters(imat,acden,acvis)          
     
         !Hook
         call ProcHook%PhysicalProp  
         
         !Advection velocity      
         call ProcPointer%ComputeAdvectionVelocity    
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm       

         
         !Compute contributions to elemental matrix : Block U,V
         val_aux = LHSdtinv*acden*dvol
         call elmmas(e,val_aux,wrmat1)
         
         !Hook
         call ProcHook%InGaussElmats         
         
         

      end do gauss_points
      
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
      
      
      !Prescribe Dirichlet boundary conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
      !inside the subroutine
      ! php_elmdir(a,e,ndofn,ndofbc,ndofbcstart,currentbvess,elmat,elrhs)
      currentbvess=auxtens+1
      bcstar=0_ip
      if(a%kfl_exacs/=0.and.a%kfl_timei==1)then
         call sup_Exaelmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)
      elseif(a%kfl_confi==1.and.a%kfl_bc_number>0)then
         call sup_Exaelmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)
      else
         call php_elmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)
      end if  
   
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   end do elements
   
   !Hook
   call ProcHook%Finalizations   
   
   !Memory Deallocations
   call a%Memor%Dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','supf_elmope_4rd')
   call a%Memor%Dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','supf_elmope_4rd')
   call a%Memor%Dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','supf_elmope_4rd')   
   call a%Memor%Dealloc(e%ndime,e%mnode,elrhs,'elrhs','supf_elmope_4rd')
   
   call a%Memor%Dealloc(e%ndime,e%mnode,1,elvel,'elvel','supf_elmope_4rd')
   call a%Memor%Dealloc(e%mnode,2,elpre,'elpre','supf_elmope_4rd')
   call a%Memor%Dealloc(e%mnode,1,elpreInc,'elpreInc','supf_elmope_4rd')
   call a%Memor%Dealloc(e%ndime,gpvel,'gpvel','supf_elmope_4rd')
   
   call a%Memor%dealloc(auxtens,e%mnode,a%ncomp-1,elsig,'elsig','supf_elmope_4rd')      
   call a%Memor%dealloc(auxtens,a%ncomp-1,gpsig,'gpsig','supf_elmope_4rd')     
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supf_elmope_4rd')
end subroutine supf_elmope_4rd



