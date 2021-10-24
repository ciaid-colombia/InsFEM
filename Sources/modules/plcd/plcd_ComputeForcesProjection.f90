subroutine plcd_ComputeForcesProjection(a,ForcesProjection)
   use typre
   use Mod_PLCD
   use Mod_plcd_ExternalForces
   use Mod_Element
   implicit none
   class(PLCDProblem) :: a
   real(rp) :: ForcesProjection(:,:)
   
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: nelem,ielem,igaus
   real(rp) :: dvol,elext(4)
   
   real(rp), allocatable :: GaussElExternalForces(:,:)
   integer(ip) :: inode
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','L2Projector')
   
   call a%Memor%alloc(a%ndofn,e%mnode,GaussElExternalForces,'GaussElExternalForces','plcd_ComputeForcesProjection')
   
   ForcesProjection = 0.0_rp
   call a%Mesh%GetNelem(nelem)
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      call e%elmdel
      
      GaussElExternalForces = 0.0_rp
      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
      
         call GetExternalForces(e,a,elext)
         
         do inode = 1,e%pnode
            GaussElExternalForces(1:a%ndofn,inode) = GaussElExternalForces(1:a%ndofn,inode) + elext(1:a%ndofn)*e%shape(inode,e%igaus)*dvol
         enddo
      enddo
      
      call a%Mesh%AssemblyToArray(e,a%ndofn,GaussElExternalForces,ForcesProjection)
      
   enddo  
   
   call a%Mesh%Project(a%ndofn,ForcesProjection) 
   
   
   
   call a%Memor%dealloc(a%ndofn,e%mnode,GaussElExternalForces,'GaussElExternalForces','plcd_ComputeForcesProjection')
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','L2Projector')
   
   
   
   
end subroutine

