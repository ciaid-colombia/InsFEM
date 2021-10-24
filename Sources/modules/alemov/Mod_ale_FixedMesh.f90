module Mod_ale_FixedMesh
   use Mod_PointerSetter
   use Mod_Alemov
   use Mod_ale_BaseElmope
   implicit none
   private
   public SetPointersFixedMesh
   
   type, extends(PointerSetter) :: SPFixedMesh
contains
      procedure :: SpecificSet => SpecificSetFixedMesh
   end type
   
   type(SPFixedMesh) :: SetPointersFixedMesh
   
   integer(ip) :: elemStatus
   real(rp) :: pdsurf, weigp(10)
   integer(ip) :: nipoin
   real(rp), allocatable :: eldisp(:,:)
   real(rp) :: gpdis(1)
   real(rp), pointer :: xloc(:,:)
   integer(ip) :: poinStatus,ndime,igaus,jnode,kpoin
   
   integer(ip) :: method = 2    !0: Fix surrounding nodes
                                !1: Art: Approximate Imposition of Boundary Conditions in Immersed Boundary Methods
                                !2: Nitsche's method (penalty, no need for flux integrals)
   
contains   
   subroutine SpecificSetFixedMesh(d)
      implicit none
      class(SPFixedMesh) :: d
      
      !Set the pointers
      if (a%kfl_IsFixedMesh == 1) then
         call ConcatenateProcedures(ProcHook%Initializations,Initializations)
         call ConcatenateProcedures(ProcHook%PostGaussElmats,PostGaussElmats)
         call ConcatenateProcedures(ProcHook%Finalizations,Finalizations)
      endif
   end subroutine
   

   subroutine Initializations
      implicit none
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      
      call a%Memor%alloc(a%ndofn,e%mnode,eldisp,'eldisp','ale_elmope')
      
      !Fixe surrounding nodes
      if (method == 0) then
         !For Fixed-Mesh ALE we fix the velocity of the first layer of nodes 
         !to be Lagrangian
         if (a%kfl_IsFixedMesh == 1) then
            do ipoin = 1,npoin
               call a%CutMesh%GetPointType(ipoin,poinStatus) 
               if (poinStatus == 1 .or. poinStatus == -1) then
                  a%kfl_fixno(currentbvess,ipoin) = 1
               endif
            enddo
         endif
      endif
      
   end subroutine
   
   subroutine PostGaussElmats
      implicit none
      
      !Nothing to do if prescribing only adjacent nodes
      if (method == 0) return
   
      !Art Approximate imposition of boundary conditions in immersed boundary methods
      if (method == 1) then
         !Fixed Mesh ALE, we need to fix the velocity of the nodes at the interface
         !Do not contribute to nodes dedicated to the boundary conditions, on all elements
         call ExternalNodesToZero
      endif
      
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      !Cut elements, boundary condition
      if(elemStatus ==0)then   
      
         call ComputeCutData
         
      
         do inode = 1,e%pnode
            ipoin = e%lnods(inode)
            
            
            call a%CutMesh%GetPointType(ipoin,poinStatus)
            
            !External nodes if method 1, or all the nodes if method 2
            if(poinStatus == -1 .or. method == 2)then
            
               call AddPenalty
            
            endif
         enddo
      endif
      
      
      
      
contains
      subroutine ComputeCutData
         implicit none
         
         !Gather prescribed displacements at the nodes
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            eldisp(1_ip,inode) = a%bvess(currentbvess,ipoin,1)
         end do
         
         call e%elmdcg
         call e%elmlen
         
         call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
         call a%CutMesh%GetSurfaceIntersection(ielem,e,pdsurf)
         call a%CutMesh%GetNinters(ielem,nipoin)
      
         !Compute shape functions at the intersection points
         do kpoin = 1, nipoin         
            weigp(kpoin) = 1.0_rp/nipoin
         end do      
            
         !The rutine gives the needed shape functions associated 
         call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)
      end subroutine
      
      subroutine ExternalNodesToZero
         implicit none
         
         do inode = 1,e%pnode
            ipoin = e%lnods(inode)
            call a%CutMesh%GetPointType(ipoin,poinStatus)
            
            if(poinStatus == -1)then
               elmat(:,inode,:,:) = 0.0_rp
               elrhs(:,inode) = 0.0_rp
            endif
         enddo
      end subroutine
      
      subroutine AddPenalty
         implicit none
         
         
         
         do igaus = 1,e%pgaus
            e%igaus = igaus  
            do jnode = 1,e%pnode
               elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
               + e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*weigp(e%igaus)*pdsurf/e%hleng(1)*1e4_rp
            enddo
            call e%interpg(1_ip,eldisp,gpdis)     
            elrhs(1,inode) = elrhs(1,inode) + e%shape(inode,e%igaus)*weigp(e%igaus)*pdsurf*gpdis(1)/e%hleng(1)*1e4_rp
         enddo
      end subroutine

   
   end subroutine
   
   subroutine Finalizations
      implicit none
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      
      call a%Memor%dealloc(a%ndofn,e%mnode,eldisp,'eldisp','ale_elmope')
   end subroutine
   
end module
