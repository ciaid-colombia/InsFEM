module Mod_nsm_LevelSetCoupling
   use typre
   use Mod_PointerSetter
   use Mod_CutMesh
   use Mod_nsm_BaseElmope
   implicit none
   private
   public :: SetPointersLevelSetCoupling
   
   type, extends(PointerSetter) :: SPLevelSetCoupling
contains
      procedure :: SpecificSet => SpecificSetLevelSetCoupling
   end type
   type(SPLevelSetCoupling) :: SetPointersLevelSetCoupling
 
   real(rp)              :: gplev(1)
   integer(ip)           :: elemStatus(1),ngauss_minus,ngauss_plus,ngaus_total
   real(rp), allocatable :: weigp(:),xloc(:,:)
   real(rp), allocatable :: elCutGradient(:,:)
   
contains

   subroutine SpecificSetLevelSetCoupling(d)
      implicit none
      class(SPLevelSetCoupling) :: d
         
      if (a%kfl_colev == 1) then            
         !Level Set Coupling
         call ConcatenateProcedures(ProcHook_Initializations,AllocLev)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocLev)
         call ConcatenateProcedures(ProcHook_PreGauss,Cutelements)
         call PrependProcedure(ProcHook_PhysicalProp,FluidProperties) 
      endif
   end subroutine   
   
   !-------------------------------------------------------------------
   !LevelSet Coupling, Computation Subroutines
   subroutine AllocLev
      implicit none
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
      call a%Memor%alloc(e%ndime,auxdim*auxdim1,xloc,'xloc','nsm_elmope')
      if (a%kfl_SurfaceTension == 1) then
         call a%Memor%alloc(e%ndime,e%mnode,elCutGradient,'elCutGradient','AllocLev')
      endif
   end subroutine
   
   subroutine Cutelements
      implicit none
      integer(ip)  :: elemStatus,iauxgauss
      
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
   
   subroutine FluidProperties
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
    
   subroutine SurfaceTensionGathers
      implicit none
      integer(ip) :: inode
      real(rp)    :: vnor
      
      call e%gather(e%ndime,elCutGradient,a%CutGradient)
      do inode = 1,e%pnode
         call vecnor(elCutGradient(:,inode),e%ndime,vnor,2)
         if(vnor/=0.0_rp)then
            elCutGradient(:,inode) = elCutGradient(:,inode)/vnor
         end if
      enddo
   end subroutine
   
   subroutine SurfaceTensionPostGaussElmats
      implicit none
      real(rp), pointer :: ipoints(:,:) => NULL()
      real(rp), pointer :: xloc(:,:) => NULL()
      integer(ip) :: nipoin
      integer(ip) :: elemStatus
      integer(ip) :: kpoin,igaus,inode,idime
      real(rp)    :: dsurf
      real(rp)    :: gamma_surface = 0.0_rp
      real(rp)    :: CutHessian(e%ndime,e%ndime),Kappa,vnor
      real(rp)    :: gpnormal(e%ndime)
      real(rp)    :: weigp(4)
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus ==0)then 
         call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
         call a%CutMesh%GetSurfaceIntersection(ielem,e,dsurf)
         call a%CutMesh%GetNinters(ielem,nipoin)
               
         !Compute shape functions at the intersection points
         do kpoin = 1, nipoin         
            weigp(kpoin) = 1.0_rp/nipoin
         end do      
            
         !The rutine give the needed shape functions associated 
         call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp) 
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg 
         
         do igaus = 1,e%pgaus
            e%igaus = igaus
            
            call e%elmder
            call e%gradient(e%ndime,elCutGradient,CutHessian)
            !Compute Curvature
            Kappa = 0.0_rp
            do idime = 1,e%ndime
               Kappa = Kappa + CutHessian(idime,idime)
            enddo
            Kappa = -Kappa/2_rp
            call e%interpg(e%ndime,elCutGradient,gpnormal)
            call vecnor(gpnormal,e%ndime,vnor,2)
           
            if(vnor/=0.0_rp) gpnormal = gpnormal/vnor
            !Compute Surface Tension
            do inode = 1,e%pnode
               do idime = 1,e%ndime
                  elrhu(idime,inode) = elrhu(idime,inode) +  gamma_surface*Kappa*e%shape(inode,e%igaus)*gpnormal(idime)*dsurf
               enddo
            enddo
         enddo
      endif
   end subroutine   
   
   subroutine DeallocLev
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
      call a%Memor%dealloc(e%ndime,auxdim*auxdim1,xloc,'xloc','nsm_elmope') 
      if (a%kfl_SurfaceTension == 1) then
         call a%Memor%dealloc(e%ndime,e%mnode,elCutGradient,'elCutGradient','AllocLev')
      endif
   end subroutine

end module
