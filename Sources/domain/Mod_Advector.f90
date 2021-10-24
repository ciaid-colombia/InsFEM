module Mod_Advector
   use typre
   use Mod_ParallelLibraryInterface
   use Mod_ParallelSystemInterface
   use Mod_MPIObject
   use Mod_Memor
   use Mod_TimeIntegrator
   use Mod_Element
   use Mod_Mesh
   use Mod_ConvectiveElement
   use Mod_TemperatureElement
   implicit none
   private
   public Advector
   
   real(rp), parameter :: staco(4) = (/400.0_rp, 200.0_rp, 100.0_rp, 0.0_rp/)

   type, extends(MPIObject) :: Advector
      type(FemMesh), pointer           :: Mesh => NULL()
      type(MemoryMan), pointer :: Memor2 => NULL()
      class(ParallelLibraryInterface), pointer   :: ParallelLibrary => NULL()
      class(ParallelSystemInterface), pointer   :: LinearSystem => NULL()
      
      real(rp), pointer :: AdvectionVelocity(:,:) => NULL()
      real(rp) :: dtime

      
contains
      procedure :: SetMesh
      procedure :: SetParallelLibrary
      procedure :: SetMemor
      procedure :: Initialize
      procedure :: Finalize
      
      
      procedure :: BuildMatrix
      procedure :: BuildRHS
      procedure :: Advect
      
   end type


   contains

   subroutine SetMesh(a,Mesh)
      class(Advector) :: a
      type(FemMesh), target :: Mesh
      a%Mesh => Mesh
   end subroutine
   
   subroutine SetParallelLibrary(a,ParallelLibrary)
      class(Advector) :: a
      class(ParallelLibraryInterface), pointer :: ParallelLibrary
      a%ParallelLibrary => ParallelLibrary
   end subroutine
   
   subroutine SetMemor(a,Memor2)
      class(Advector) :: a
      type(MemoryMan), target :: Memor2
      a%Memor2 => Memor2
   end subroutine
   
   subroutine Initialize(a)
      class(Advector) :: a
      
      character(150) :: auxstring
      character(150) :: auxstring2
      
      auxstring = trim(a%Mesh%InputFolder)//'/'//adjustl(trim(a%Mesh%namda))//'.dom.sol '
      auxstring2 = 'dom '
  
      call a%ParallelLibrary%CreateSystem(a%LinearSystem, a%Memor2)
      call a%Mesh%InitializeSystem(1_ip,0_ip,a%LinearSystem,a%Memor2,auxstring,auxstring2,a%Mesh%lun_solve_dom)
   end subroutine
   
   subroutine Finalize(a)
      class(Advector) :: a
      
      call a%LinearSystem%Deallocate
      call a%ParallelLibrary%DeallocateSystem(a%LinearSystem, a%Memor2) 
   end subroutine
   
   

   subroutine BuildMatrix(a,dt,AdvectionVelocity)
      implicit none
      class(Advector) :: a
      real(rp) :: dt
      real(rp), target :: AdvectionVelocity(:,:)
      
      
      
      class(FiniteElement), pointer :: e => NULL()
      type(TimeIntegratorDt1) :: Integrator
      integer(ip)           :: nsteps
      real(rp)              :: LHSdtinv
      real(rp)              :: dsurf
      real(rp), allocatable :: elmat(:,:,:,:),elrhs(:,:)
      real(rp)              :: dvol
      
      real(rp), allocatable :: elvel(:,:),grvel(:,:)
      real(rp), allocatable :: AGradV(:),gpvel(:)
      real(rp)              :: gpvno
      real(rp), allocatable :: testf(:)
      
      real(rp) :: chale(2), timom
      integer(ip) :: igaus,ielem,nelem
      real(rp) :: dtinv
      integer(ip) :: kfl_HangingNodes

      a%AdvectionVelocity => AdvectionVelocity
      a%dtime = dt
      
      
      !Linear system to zero
      call a%LinearSystem%ToZeroMatrix
      
      !Time integrator
      dtinv = 1.0_rp/dt
      call Integrator%Init('BDF1 ')
      call Integrator%GetLHSDtinv(dtinv,LHSdtinv)
      
      
      
      call a%Mesh%ElementAlloc(e,a%Memor2,'DefaultRule','nsm_EnditeElmope')
      call a%Memor2%alloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsm_elmope')
      call a%Memor2%alloc(1_ip,e%mnode,elrhs,'elrhs','nsm_elmope')
      
      call a%Memor2%alloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
      call a%Memor2%alloc(e%ndime,gpvel,'gpvel','tem_elmope')
      call a%Memor2%alloc(e%mnode,AGradV,'AGradV','tem_elmope')
      call a%Memor2%alloc(e%mnode,testf,'testf','tem_Elmope')
      
      call a%Mesh%GetHanging(kfl_HangingNodes)
      
      call a%Mesh%GetNelem(nelem)
      elements : do ielem = 1,nelem
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)  
         
         !ElmatsToZero
         elmat=0.0_rp
      
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         
         !Element length at center of gravity
         call e%elmlen
         
         call e%gather(e%ndime,elvel(:,:),AdvectionVelocity)

         ! Compute the characteristic length chale
         call elmchl(e,1_ip,elvel,chale)
         
         !Gauss Point Loop
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            
            call e%elmder
         
            dvol = e%weigp(e%igaus)*e%detjm
            
            !Interpolate velocity
            call e%interpg(e%ndime,elvel(:,:),gpvel(:))
            
            !Advection velocity norm
            call vecnor(gpvel,e%ndime,gpvno,2)
         
            !Compute a·grad(V)
            call ComputeAGradV(e,gpvel,AGradV)
            
            !Compute the stability parameters
            call ComputeTauCDR(e,1.0_rp,0.0_rp,0.0_rp,gpvno,staco,chale,timom)
            
            !Adjoint Test Function
            !Hook
            call tem_ComputeTestf(e,1.0_rp,timom,AGradV,0.0_rp,testf)
                  
            !Compute contributions to elemental matrix : Block U,V
            call tem_elmbuv(e,dvol,1.0_rp,0.0_rp,LHSdtinv,AGradV,testf,elmat(1,:,1,:))
         enddo gauss_points
         
         if (kfl_HangingNodes == 1) call a%Mesh%PrepareHangingMatrices(e,elmat,elrhs)  
         call a%LinearSystem%AssemblyElmat(e,elmat)
      enddo elements
      
      if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(1_ip,a%LinearSystem,a%Memor2)
 
 
      
      call a%Memor2%Dealloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
      call a%Memor2%Dealloc(e%ndime,gpvel,'gpvel','tem_elmope')
      call a%Memor2%Dealloc(e%mnode,AGradV,'AGradV','tem_elmope')
      call a%Memor2%dealloc(e%mnode,testf,'testf','tem_Elmope')
      call a%Memor2%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsm_elmope')
      call a%Memor2%dealloc(1_ip,e%mnode,elrhs,'elrhs','nsm_elmope')
      call a%Mesh%ElementDeAlloc(e,a%Memor2,'DefaultRule','nsm_EnditeElmope')
      
   end subroutine
   
   
   subroutine Advect(a,ndofn,vectorToAdvect,AdvectedVector)
      class(Advector) :: a
      integer(ip) :: ndofn
      real(rp) :: VectorToAdvect(:,:), AdvectedVector(:,:)
      
      integer(ip) :: idofn
      integer(ip) :: kfl_HangingNodes,kfl_perio
      
      call a%Mesh%GetHanging(kfl_HangingNodes)
      do idofn = 1,ndofn
         call a%BuildRHS(VectorToAdvect(idofn,:))
         call a%LinearSystem%Solve(AdvectedVector(idofn,:))
      enddo
      
      !Ghostcommunicate
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndofn,AdvectedVector)

      !HangingNodes
      call a%Mesh%GetHanging(kfl_HangingNodes)
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(ndofn,AdvectedVector)

      !Periodic boundary conditions
      call a%Mesh%GetPerio(kfl_perio)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(ndofn,AdvectedVector)

   end subroutine
      
   subroutine BuildRHS(a,vectorToAdvect)
      implicit none
      class(Advector) :: a
      real(rp) :: dt
      !real(rp) :: AdvectionVelocity(:,:)
      real(rp) :: vectorToAdvect(:)

      class(FiniteElement), pointer :: e => NULL()
      type(TimeIntegratorDt1) :: Integrator
      integer(ip)           :: nsteps
      real(rp)              :: LHSdtinv
      real(rp)              :: dsurf
      real(rp), allocatable :: elmat(:,:,:,:),elrhs(:,:)
      real(rp)              :: dvol
      
      real(rp), allocatable :: elvel(:,:),grvel(:,:)
      real(rp), allocatable :: AGradV(:),gpvel(:)
      real(rp)              :: gpvno
      real(rp), allocatable :: testf(:)
      
      real(rp) :: chale(2), timom
      integer(ip) :: igaus,ielem,nelem
      real(rp) :: dtinv
      real(rp), allocatable :: eltem(:,:),gptem(:)
      real(rp) :: elext
      integer(ip) :: kfl_HangingNodes
      
      !Linear system to zero
      call a%LinearSystem%ToZeroRHS
      
      !Time integrator
      dtinv = 1.0_rp/a%dtime
      call Integrator%Init('BDF1 ')
      call Integrator%GetLHSDtinv(dtinv,LHSdtinv)
      
      
      
      call a%Mesh%ElementAlloc(e,a%Memor2,'DefaultRule','nsm_EnditeElmope')
      call a%Memor2%alloc(1_ip,e%mnode,elrhs,'elrhs','nsm_elmope')
      call a%Memor2%alloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsm_elmope')
      
      call a%Memor2%alloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
      call a%Memor2%alloc(e%ndime,gpvel,'gpvel','tem_elmope')
      call a%Memor2%alloc(e%mnode,AGradV,'AGradV','tem_elmope')
      call a%Memor2%alloc(e%mnode,testf,'testf','tem_Elmope')
      
      call a%Memor%alloc(e%mnode,3,eltem,'eltem','tem_elmope')
      call a%Memor%alloc(3,gptem,'gptem','tem_elmope')
      
      call a%Mesh%GetHanging(kfl_HangingNodes)
      
      call a%Mesh%GetNelem(nelem)
      elements : do ielem = 1,nelem
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)  
         
         !ElmatsToZero
         elrhs=0.0_rp
      
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         
         !Element length at center of gravity
         call e%elmlen
         
         call e%gather(e%ndime,elvel(:,:),a%AdvectionVelocity)
         !Gathers
         call e%gather(1,eltem(:,1),vectorToAdvect(:))
         !call e%gather(1,eltem(:,2),vectorToAdvect(:,2))
         !call e%gather(1,eltem(:,3),vectorToAdvect(:,3))
      

         ! Compute the characteristic length chale
         call elmchl(e,1_ip,elvel,chale)
         
         !Gauss Point Loop
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            
            call e%elmder
         
            dvol = e%weigp(e%igaus)*e%detjm
            
            !Interpolate velocity
            call e%interpg(e%ndime,elvel(:,:),gpvel(:))
            
            !Advection velocity norm
            call vecnor(gpvel,e%ndime,gpvno,2)
         
            !Compute a·grad(V)
            call ComputeAGradV(e,gpvel,AGradV)
            
            !Compute the stability parameters
            call ComputeTauCDR(e,1.0_rp,0.0_rp,0.0_rp,gpvno,staco,chale,timom)
            
            !Adjoint Test Function
            !Hook
            call tem_ComputeTestf(e,1.0_rp,timom,AGradV,0.0_rp,testf)
            
            call e%interpg(1,eltem(:,1),gptem(2))
            !call e%interpg(1,eltem(:,2),gptem(2))
            !call e%interpg(1,eltem(:,3),gptem(3))
            
            !Compute Elext, Temporal Derivatives
            elext=0.0_rp
            !Time integration
            call tem_TimeIntegrationToElext(e,Integrator,1.0_rp,dtinv,gptem,elext)
                    
            !Compute contributions to elemental matrix : Block U,V
            call tem_elmrhu(e,dvol,testf,elext,elrhs)
         enddo gauss_points
         
         if (kfl_HangingNodes == 1) call a%Mesh%PrepareHangingMatrices(e,elmat,elrhs)
            
         call a%LinearSystem%AssemblyElrhs(e,elrhs)
      enddo elements
        
      call a%Memor2%Dealloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
      call a%Memor2%Dealloc(e%ndime,gpvel,'gpvel','tem_elmope')
      call a%Memor2%Dealloc(e%mnode,AGradV,'AGradV','tem_elmope')
      call a%Memor2%dealloc(e%mnode,testf,'testf','tem_Elmope')
      call a%Memor2%dealloc(1_ip,e%mnode,elrhs,'elrhs','nsm_elmope')
      call a%Memor2%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsm_elmope')
      
      call a%Mesh%ElementDeAlloc(e,a%Memor2,'DefaultRule','nsm_EnditeElmope')
      
   end subroutine

end module
