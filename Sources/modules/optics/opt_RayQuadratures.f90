subroutine opt_RayQuadratures(a,cn2,cn2u53,cn2_quadrature,cn2u53_quadrature,length_quadrature,r0,vwind,fg)
   use typre
   use Mod_Memor
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_Optics
   use Mod_Element
   use Mod_RayElementIntersection
   use MPI
   implicit none
   class(OpticsProblem) :: a
   real(rp) :: cn2(*), cn2u53(*)
   real(rp) :: cn2_quadrature(a%nbeams), length_quadrature(a%nbeams),cn2u53_quadrature(a%nbeams)
   real(rp) :: r0(a%nbeams), fg(a%nbeams), vwind(a%nbeams)
   
   integer(ip) :: ierr
   
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ibeam,ielem,nelem
   real(rp)    :: intersection(3,2)
   logical     :: kfl_intersect
   
   integer(ip) :: iinter
   
   real(rp), allocatable :: elcn2(:),elcn2u53(:)
   real(rp) :: gpcn2(1),gpcn2u53(1)
   
   integer(ip) :: npoinLocal,inode,istat
   real(rp)    :: weighting
   
   real(rp) :: length
   real(rp) :: length_vec(3)
   real(rp) :: xloc(3)
   
   real(rp) :: cn2_quadratureb(a%nbeams), length_quadratureb(a%nbeams),cn2u53_quadratureb(a%nbeams)
   
   
   character(150) :: fil_outpu
   integer(ip)    :: lun_outpu
   
   
   
   !Global values
   call a%Mesh%GetNpoinLocal(npoinLocal)
   
   !Allocations
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','opt_elmope')
   
   call a%Memor%alloc(e%mnode,elcn2,'elcn2','opt_TurnofComputations')
   call a%Memor%alloc(e%mnode,elcn2u53,'elcn2u53','opt_TurnofComputations')
   
   !Initializations
   cn2_quadratureb = 0.0_rp
   cn2u53_quadratureb = 0.0_rp
   length_quadratureb = 0.0_rp
   
   
   !Now loop through elements, see if each of the beams intersects the elements
   !and compute the contribution to the integrals in case of intersected elements
   call a%Mesh%GetNelem(nelem)
   do ielem = 1,nelem

      !Load Element
      call a%Mesh%ElementLoad(ielem,e)
      
      !Gather
      call e%gather(1,elcn2,cn2)
      call e%gather(1,elcn2u53,cn2u53)
      
      !Compute weighting factor for the quadrature
      !This takes into account the possibility of a ray traversing an element with nodes in multiple processors
      weighting = 0.0_rp
      do inode = 1,e%pnode
         if (e%lnods(inode) <= npoinLocal) then
            weighting = weighting+1.0_rp
         endif
      enddo
      weighting = weighting/real(e%pnode)
      
      !For each beam
      do ibeam = 1,a%nbeams
         call RayElementIntersection(e,a%beams(ibeam),kfl_intersect,intersection)
         if (kfl_intersect .eqv. .true.) then
            !Get integration length
            length_vec(1:e%ndime) = intersection(1:e%ndime,1)-intersection(1:e%ndime,2)
            call vecnor(length_vec,e%ndime,length,2)
         
            do iinter = 1,2
               !Get Gauss points coordinates
               !At this point Gauss points coordinates are the intersection points (linear elements)
               !Interpolate cn2 values, and velocity module, at gauss points coordinates
               !Add contribution to the beam cuadratures (cn2 and cn2*|u|**5/3)
               call e%isoparinv(intersection(:,iinter),xloc)
               call e%isopar(1,xloc,elcn2,gpcn2)
               call e%isopar(1,xloc,elcn2u53,gpcn2u53)
               
               cn2_quadratureb(ibeam) = cn2_quadratureb(ibeam) + gpcn2(1)*length/2*weighting;
               cn2u53_quadratureb(ibeam) = cn2u53_quadratureb(ibeam) + gpcn2u53(1)*length/2*weighting;
               length_quadratureb(ibeam) = length_quadratureb(ibeam) + length/2*weighting;
            enddo   
            
         endif   
      enddo
   enddo
   
   call  MPI_REDUCE(cn2_quadratureb, cn2_quadrature, a%nbeams, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE(cn2u53_quadratureb, cn2u53_quadrature, a%nbeams, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE(length_quadratureb, length_quadrature, a%nbeams, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   
   if (a%MPIrank == a%MPIroot) then
      do ibeam = 1,a%nbeams
         if (cn2_quadrature(ibeam) <= 0.0_rp .or. cn2u53_quadrature(ibeam) <= 0.0_rp) then
            write(*,*) 'Warning: bad values for integral of cn2 in beam ',ibeam
            r0(ibeam) = 0.0_rp
            vwind(ibeam) = 0.0_rp
            fg(ibeam) = 0.0_rp
         else
            r0(ibeam) = (16.6/(a%lambda**2)*cn2_quadrature(ibeam))**(-3.0_rp/5.0_rp)
            vwind(ibeam) = (cn2u53_quadrature(ibeam)/cn2_quadrature(ibeam) )**(3.0_rp/5.0_rp)
            fg(ibeam) = 0.43*vwind(ibeam)/r0(ibeam)
         endif   
     enddo
   endif
   
   !Memory deallocations
   call a%Memor%dealloc(e%mnode,elcn2,'elcn2','opt_TurnofComputations')
   call a%Memor%dealloc(e%mnode,elcn2u53,'elcn2u53','opt_TurnofComputations')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','opt_elmope')
   
end subroutine
