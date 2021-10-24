module Mod_SmoothedFieldGradient

contains

   subroutine ComputeSmoothedFieldGradient(Mesh,Memor,ndofn,Displacement,SmoothedG)
      use typre
      use Mod_Mesh
      use Mod_Element
      use Mod_Memor
      implicit none
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: ndofn 
      real(rp) :: Displacement(ndofn,*)
      real(rp), target :: SmoothedG(*)
      
      real(rp), pointer :: SmoothedGradient(:,:,:) => NULL()
      
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable :: gradDisp(:,:),eldisp(:,:),elrhsSG(:,:,:)
      
      real(rp) :: dvol
      
      integer(ip) :: idofn,inode,ipoin,npoin,idime,ielem,nelem,igaus
      
      call Mesh%GetNpoin(npoin)
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','ComputeSmoothedGradient')
      call Memor%Alloc(ndofn,e%mnode,eldisp,'eldisp','plcd_EnditeElmope')
      call Memor%Alloc(ndofn,e%ndime,gradDisp,'gradDisp','plcd_EnditeElmope')
      
      call Memor%Alloc(ndofn,e%ndime,e%mnode,elrhsSG,'elrhsSG','plcd_EnditeElmope')

      SmoothedGradient(1:ndofn,1:e%ndime,1:npoin) => SmoothedG(1:ndofn*e%ndime*npoin)
      SmoothedGradient(:,:,:) = 0.0_rp
      
      
      call Mesh%GetNelem(nelem)
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)  
         
         call e%gather(ndofn,eldisp,Displacement)
         
         elrhsSG = 0.0_rp
         
         !Compute linear derivatives
         call e%elmdel
         
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
         
         
            call e%elmder
            dvol = e%detjm*e%weigp(e%igaus)

            call e%gradient(ndofn,eldisp,gradDisp)
            do inode = 1,e%pnode
               do idofn = 1,ndofn
                  do idime = 1,e%ndime
                     elrhsSG(idofn,idime,inode) = elrhsSG(idofn,idime,inode) + e%shape(inode,e%igaus)*dvol*gradDisp(idofn,idime)
                  enddo
               enddo
            enddo
            
         enddo gauss_points
         
         call Mesh%AssemblyToArray(e,ndofn*e%ndime,elrhsSG,SmoothedGradient)
      enddo elements
      
      
      call Memor%dealloc(ndofn,e%ndime,e%mnode,elrhsSG,'elrhsSG','SmoothGradientInitializations')
      
      call Memor%deAlloc(ndofn,e%mnode,eldisp,'eldisp','plcd_EnditeElmope')
      call Memor%deAlloc(ndofn,e%ndime,gradDisp,'gradDisp','plcd_EnditeElmope')
      
      call Mesh%Smooth(ndofn*e%ndime,SmoothedGradient)
      
      !DeallocateElement
      call Mesh%ElementDeAlloc(e,Memor,'ForceClosedRule','plcd_Elmope')

      
      
      
   end subroutine
   
   subroutine PostprocessGaussPointGradient(Mesh,Memor,ndofn,Displacement,fieldname,FilePostpr,istep,ctime)
      use typre
      use Mod_Memor
      use Mod_Element
      use Mod_Mesh
      use Mod_Postpr
      use Mod_RPMeshAllocator
      implicit none
      class(FemMesh) :: Mesh
      type(MemoryMan) :: Memor
      integer(ip) :: ndofn 
      real(rp) :: Displacement(ndofn,*)
      character(*) :: fieldname
      class(PostprFile) :: FilePostpr
      integer(ip) :: istep
      real(rp) :: ctime
      
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable :: gradDisp(:,:),eldisp(:,:)
      type(r3p), allocatable :: GaussPointField(:)
      
      
      real(rp) :: dvol
      
      integer(ip) :: idofn,inode,ipoin,npoin,idime,ielem,nelem,igaus,ndime
      
      call Mesh%GetNpoin(npoin)
      call Mesh%ElementAlloc(e,Memor,'ForceClosedRule','ComputeSmoothedGradient')
      call Memor%Alloc(ndofn,e%mnode,eldisp,'eldisp','plcd_EnditeElmope')
      call Memor%Alloc(ndofn,e%ndime,gradDisp,'gradDisp','plcd_EnditeElmope')


      call Mesh%GetNdime(ndime)
      call AllocR3P(Memor,Mesh,ndofn,ndime,GaussPointField)
      
      
      call Mesh%GetNelem(nelem)
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)  
         
         call e%gather(ndofn,eldisp,Displacement)

         
         !Compute linear derivatives
         call e%elmdel
         
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
         
         
            call e%elmder
            dvol = e%detjm*e%weigp(e%igaus)

            call e%gradient(ndofn,eldisp,gradDisp)
            
            GaussPointField(ielem)%a(:,:,e%igaus) = gradDisp
            
         enddo gauss_points

      enddo elements
      
      !Postprocess
      call FilePostpr%postgp(GaussPointField,fieldname,istep,ctime,Mesh)

      
      call Memor%deAlloc(ndofn,e%mnode,eldisp,'eldisp','plcd_EnditeElmope')
      call Memor%deAlloc(ndofn,e%ndime,gradDisp,'gradDisp','plcd_EnditeElmope')
      
      !DeallocateElement
      call Mesh%ElementDeAlloc(e,Memor,'ForceClosedRule','plcd_Elmope')


      
      call Mesh%GetNdime(ndime)
      call DeAllocR3P(Memor,Mesh,ndime,ndofn,GaussPointField)
   
   end subroutine
end module
