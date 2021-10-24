module Mod_nsm_FreeSurfaceMatrices
   use typre
   use Mod_PointerSetter
   use Mod_CutMesh
   use Mod_nsm_BaseElmope
   use Mod_nsm_HangingNodes, only : SetPointersHangingNodes
   use Mod_nsm_InterpolateGradients
   implicit none
   private
   public SetPointersFreeSurfaceMatrices
   
   type, extends(PointerSetter) :: SPFreeSurfaceMatrices
contains
      procedure :: SpecificSet => SpecificSetFreeSurfaceMatrices
   end type
   type(SPFreeSurfaceMatrices) :: SetPointersFreeSurfaceMatrices
 
   real(rp)    :: gplev(1)
   integer(ip) :: elemStatus(1),ngauss_minus,ngauss_plus,ngaus_total
   integer(ip) :: kfl_isHanging
   integer(ip), allocatable :: HangingStatus(:)
   
   !DirichletBoundaryConditions
   real(rp), allocatable :: elDirichletVeloc(:,:)
   real(rp) :: gpDirichletVeloc(3),gpDirNor
   class(FiniteElement), pointer :: DirichletElement => NULL()
   
contains

   subroutine SpecificSetFreeSurfaceMatrices(d)
      implicit none
      class(SPFreeSurfaceMatrices) :: d
         
      if (a%kfl_colev == 1) then        
         if(a%kfl_fsurf==1)then
            !I need to make sure that, if hanging nodes, then the hanging nodes modification is done  BEFORE I modify the matrices
            call a%Mesh%GetHanging(kfl_isHanging)
            if (kfl_isHanging == 1) then
               call SetPointersHangingNodes%Set
            endif
            
            call ConcatenateProcedures(ProcHook_PreDirichlet,FreeSurfMats)
            if(a%kfl_fsurfLapla==1)then
               call PrependProcedure(ProcHook_PreDirichlet,ElmatsToLapla)                     
            elseif (a%kfl_fsurfLapla == 2) then
               call runend('Stokes not ready for kfl_fsurflapla, two fileds incompressible navier stokes')
            elseif (a%kfl_fsurfLapla == 3) then
               call PrependProcedure(ProcHook_PreDirichlet,ElmatsToDiagonalBlock)                     
            end if
            
            if (a%kfl_fsurfDirichlet > 0) then
               call ConcatenateProcedures(ProcHook_Initializations,AllocDirichlet)
               call ConcatenateProcedures(ProcHook_Finalizations,DeAllocDirichlet)
               !call ConcatenateProcedures(ProcHook_PreGauss,DirichletPreGauss)
               call ConcatenateProcedures(ProcHook_PreDirichlet,DirichletPreGauss)
            endif
            
            call a%Mesh%GetHanging(kfl_isHanging)
            if (kfl_isHanging == 1) then
               call ConcatenateProcedures(ProcHook_Initializations,AllocHanging)
               call ConcatenateProcedures(ProcHook_Finalizations,DeAllocHanging)
            endif
         end if      
      endif

   end subroutine   
   
   !Actual Computations
   subroutine FreeSurfMats
      implicit none 
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==-1)then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%CutMesh%GetPointType(ipoin,poinStatus)            
            
            !For HangingNodes case only
            if (kfl_isHanging == 1) poinStatus = HangingStatus(ipoin)
            
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
         
         elmat(:,:,:,:) = 0.0_rp
         elrhs(:,:)=0.0_rp
         
         wrmat1=0.0_rp                  
         elmat=0.0_rp
         
         call a%GetPhysicalParameters(1_ip,acden,acvis)
         
         ! Viscosity terms : we only consider mu*(grad v, grad u)         
         call elmvis(e,dvolt0,1.0_rp,wrmat1)         
         
         forall (idime = 1:e%ndime)
            elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + acvis*wrmat1(1:e%pnode,1:e%pnode)
         end forall   
      
         !Do not do it in fractional step 1st stage
         if (size(elmat,1) > e%ndime) then
            elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) + timom*wrmat1(1:e%pnode,1:e%pnode)   
            
            do inode = 1,e%pnode
               !Rhs
               elrhs(e%ndime+1,inode) = elrhs(e%ndime+1,inode) + timom*dot_product(e%cartd(:,inode),acden*a%grnor*a%gravi(1:e%ndime))*dvolt0 
            enddo 
            
         endif
      end if      
      
    end subroutine

   subroutine ElmatsToDiagonalBlock
      implicit none
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin,idime 
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      
      if(elemStatus==-1)then
         !Eliminate off diagonal blocks for the fake fluid
         elmat(1:e%ndime,:,e%ndime+1,:) = 0.0_rp
         elmat(e%ndime+1,:,1:e%ndime,:) = 0.0_rp
      end if      
      
   end subroutine
   
   !HangingNodes
   subroutine AllocHanging
      implicit none
      
      integer(ip) :: npoin
      integer(ip) :: ipoin,phang
      integer(ip), pointer :: lhang(:) => NULL()
      
      integer(ip) :: poinStatus,inode
      
      
      interface
         subroutine HangingGetParents(a,ipoin,phang,lhang)
            use typre
            use Mod_Mesh
            implicit none
            class(FemMesh) :: a
            integer(ip) :: ipoin,phang
            integer(ip), pointer :: lhang(:)
         end subroutine
      end interface
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetHanging(kfl_isHanging)
      
      
      if (kfl_isHanging == 1) then
         call a%Memor%alloc(npoin,HangingStatus,'HangingStatus','AllocHanging')
         do ipoin = 1,npoin
            call a%CutMesh%GetPointType(ipoin,poinStatus)            
            if(poinStatus==-1)then
               call HangingGetParents(a%Mesh,ipoin,phang,lhang)
               HangingStatus(lhang) = -1
               HangingStatus(ipoin) = -1
            endif
         enddo
      endif
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,HangingStatus)
      
   end subroutine
      
   subroutine DeallocHanging
      implicit none
      
      integer(ip) :: npoin
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetHanging(kfl_isHanging)
      
      
      if (kfl_isHanging == 1) then
         call a%Memor%dealloc(npoin,HangingStatus,'HangingStatus','AllocHanging')
      endif
   end subroutine
   
   !Dirichlet Boundary Conditions
   subroutine AllocDirichlet
      implicit none
      
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,e%mnode,elDirichletVeloc,'elDirichletVeloc','AllocDirichlet')
      
      !allocate
      call a%Mesh%ElementAlloc(DirichletElement,a%Memor,'DefaultRule','lev_ReinitLevel')
      
   end subroutine
   
   subroutine DeAllocDirichlet
      implicit none
      
      call a%Memor%dealloc(ndime,e%mnode,elDirichletVeloc,'elDirichletVeloc','AllocDirichlet')
      call a%Mesh%ElementDeAlloc(DirichletElement,a%Memor,'DefaultRule','lev_ReinitLevel')
   end subroutine
   
      
   subroutine DirichletPreGauss
      implicit none
      
      real(rp), pointer :: xloc(:,:) => NULL()
      real(rp) :: weigp(3),UnitOrthogonalVector(3)
      
      integer(ip) :: idime,inode,jnode,kpoin,nipoin,kdime
      real(rp) :: pdsurf,stabterm
      
      real(rp) :: xmuit,grvelnorm
      integer(ip) :: ipoin,pstatus
      
      call a%CutMesh%GetElementType(ielem,elemStatus(1))
      
      
      if(elemStatus(1) ==0)then
         if (a%kfl_fsurfDirichletMethod > 2) then
            call runend('method not implemented')
         endif
      
         call a%Mesh%ElementLoad(ielem,DirichletElement)  
         if (a%kfl_fsurfDirichlet==1) then
            call DirichletElement%gather(ndime,elDirichletVeloc,a%DirichletVelocity)  
         elseif (a%kfl_fsurfDirichlet==2) then
            elDirichletVeloc = 0.0_rp
         endif
         
         !Cartesian derivatives and Jacobian at center of gravity
         call DirichletElement%elmdcg  
         
         call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
         call a%CutMesh%GetSurfaceIntersection(ielem,DirichletElement,pdsurf)
         call a%CutMesh%GetNinters(ielem,nipoin)
         
         !Compute shape functions at the intersection points
         do kpoin = 1, nipoin         
            weigp(kpoin) = 1.0_rp/nipoin
         end do      
            
         !The rutine give the needed shape functions associated 
         call DirichletElement%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)
         
         !Physical Properties
         !First material
         call a%GetPhysicalParameters(1_ip,acden,acvis)
         
         
         !Assembly the tractions in the boundary
         call a%CutMesh%GetExternalNormal(ielem,DirichletElement,UnitOrthogonalVector)
         
         call DirichletElement%elmdel
         
         do igaus = 1,DirichletElement%pgaus
            DirichletElement%igaus = igaus    
               
            call DirichletElement%elmder
            
            do inode= 1,DirichletElement%pnode
               do jnode = 1,DirichletElement%pnode
                  do idime = 1,ndime
                     !pressure
                     elmat(idime,inode,ndime+1,jnode) = elmat(idime,inode,ndime+1,jnode) &
                        +DirichletElement%shape(inode,DirichletElement%igaus)*DirichletElement%shape(jnode,DirichletElement%igaus)*weigp(DirichletElement%igaus)*pdsurf*UnitOrthogonalVector(idime)
                  enddo
               enddo
            enddo
            
            !Viscous terms
            do inode=1,DirichletElement%pnode
               xmuit=DirichletElement%shape(inode,DirichletElement%igaus)*acvis
               do idime=1,DirichletElement%ndime
                  do jnode=1,DirichletElement%pnode
                     elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode)  &
                        -xmuit*dot_product(DirichletElement%cartd(1:DirichletElement%ndime,jnode),UnitOrthogonalVector(1:DirichletElement%ndime))*weigp(DirichletElement%igaus)*pdsurf
                  end do
                  if (a%fvins>0.0_rp) then
                     do jnode=1,DirichletElement%pnode
                        do kdime=1,DirichletElement%ndime
                           elmat(idime,inode,kdime,jnode)=elmat(idime,inode,kdime,jnode)&
                                 -xmuit*DirichletElement%cartd(idime,jnode)*UnitOrthogonalVector(kdime)*weigp(DirichletElement%igaus)*pdsurf
                        end do
                     end do
                  end if
               end do
            end do
            
         enddo
         
         !If we are using Strong method, delete corresponding lines
         if (a%kfl_fsurfDirichletMethod == 2) then
            do inode = 1,e%pnode
               ipoin = e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,pStatus)
               if (pstatus == -1) then
                  elmat(1:e%ndime,inode,:,:) = 0.0_rp
                  elrhs(1:e%ndime,inode) = 0.0_rp
               endif
            enddo
         endif
         
         
        
         pstatus = -1
         !Simple penalty term
         do igaus = 1,DirichletElement%pgaus
            DirichletElement%igaus = igaus  
            call DirichletElement%elmder
            
            call DirichletElement%gradient(e%ndime,elvel,grvel) 
            call vecnor(grvel,e%ndime*e%ndime,grvelnorm,2)
            
            call DirichletElement%interpg(ndime,elDirichletVeloc,gpDirichletVeloc)
            call vecnor(gpDirichletVeloc,ndime,gpDirNor,2)
            
            !stabterm = acvis/e%hleng(1)+gpdirNor
            stabterm = acvis/e%hleng(1)+1000*grvelnorm*e%hleng(2)
            !stabterm = acvis/e%hleng(1)+1e12
            !stabterm = stabterm*200
            
            do inode= 1,DirichletElement%pnode
               !Selective if strong method
               if (a%kfl_fsurfDirichletMethod == 2) then
                  ipoin = e%lnods(inode)
                  call a%CutMesh%GetPointType(ipoin,pStatus)
               endif
               if (pstatus == -1) then
                  do jnode = 1,DirichletElement%pnode
                     do idime = 1,ndime
                        elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) &
                        + stabterm*DirichletElement%shape(inode,DirichletElement%igaus)*DirichletElement%shape(jnode,DirichletElement%igaus)*weigp(DirichletElement%igaus)*pdsurf
                     enddo   
                  enddo
                  do idime = 1,ndime      
                     elrhs(idime,inode) = elrhs(idime,inode) + stabterm*DirichletElement%shape(inode,DirichletElement%igaus)*gpDirichletVeloc(idime)*weigp(DirichletElement%igaus)*pdsurf
                  enddo   
               endif
            enddo
         enddo
         
         if (size(elmat,1) == ndime) call runend('wrong elmat dimensions, (not ready for fractional step)')
         
         
         
         
      endif
   end subroutine
   
   

end module
