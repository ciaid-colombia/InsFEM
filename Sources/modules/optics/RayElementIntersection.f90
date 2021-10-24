module Mod_RayElementIntersection
contains
subroutine RayElementIntersection(e,beam,kfl_intersect,intersection)
   use typre
   use Mod_Element
   use Mod_Optics
   implicit none
   class(FiniteElement) :: e
   type(LightBeam)     :: beam
   logical             :: kfl_intersect
   real(rp)            :: intersection(3,2)
   
   
   real(rp)            :: aux_intersection(3,4)
   
   real(rp) :: coord(e%ndime,e%mnodb)
   integer(ip) :: nface,mnodb,cfael(e%mnodb,12)
   integer(ip) :: iface,intcount,iedge,nnodb
   logical     :: iflag
   real(rp) :: iIntersection(3)
   integer(ip) :: selected(2),iinter
   
   real(rp) :: ProjRay(4)
   
   kfl_intersect = .false.
   intcount = 0
   
   !Tetra elements
   if (e%ndime == 3 .and. e%pnode == 4) then
      !According to Platis and Theoharis
      !http://graphics.di.uoa.gr/Downloads/papers/journals/p19.pdf
      
      call dombou_nface(e%ndime,e%pnode,nface)           
      nnodb = 3
      call dombou_cfael(e%mnodb,nnodb,nface,e%ndime,e%pnode,cfael)
   
      do iface = 1,nface
         coord(1:e%ndime,1:nnodb) = e%elcod(1:e%ndime,cfael(1:nnodb,iface))
         
         call RayTriangleIntersection(beam,coord,iflag,iIntersection)
         if (iflag .eqv. .true.) then
            kfl_intersect = .true.
            intcount = intcount+1
            aux_intersection(1:e%ndime,intcount) = iIntersection(1:e%ndime)
         endif
      enddo
      if (intcount == 2) then
         intersection(1:e%ndime,1:2) = aux_intersection(1:e%ndime,1:2)
      elseif (intcount == 1) then
         call runend('RayElementIntersection: found a wrong number of intersections')
      elseif (intcount > 4) then
         call runend('RayElementIntersection: tetra, number of intersect cannot be larger than 4')
      elseif (intcount >2) then
         !Select the most separated ones, dismiss the one in the middle
         do iinter = 1,intcount
            ProjRay(iinter) = dot_product(aux_intersection(1:e%ndime,iinter)-beam%origin(1:e%ndime),beam%direct(1:e%ndime))
         enddo
         selected(1) = maxloc(ProjRay,1)
         selected(2) = minloc(ProjRay,1)
         intersection(1:e%ndime,1:2) = aux_intersection(1:e%ndime,selected(1:2))
         intcount = 2
      endif
      
      if (intcount == 2) then
         do iinter = 1,intcount
            ProjRay(iinter) = dot_product(intersection(1:e%ndime,iinter)-beam%origin(1:e%ndime),beam%direct(1:e%ndime))
         enddo
         
         
         if (minval(ProjRay(1:2)) >= 0.0_rp .and. maxval(ProjRay(1:2)) <= beam%length ) then
            kfl_intersect = .true.
         elseif (minval(ProjRay(1:2)) < 0.0_rp .and. maxval(ProjRay(1:2))>=0.0_rp) then
            intersection(1:e%ndime,minloc(ProjRay(1:2),1)) = beam%origin(1:e%ndime)
            kfl_intersect = .true.
         elseif ( minval(ProjRay(1:2)) >= 0.0_rp .and. maxval(ProjRay(1:2)) > beam%length) then
            intersection(1:e%ndime,maxloc(ProjRay(1:2),1)) = beam%origin(1:e%ndime) + beam%direct(1:e%ndime)*beam%length
            kfl_intersect = .true.   
         else
            kfl_intersect = .false.
         endif   
         
      endif
         
         
      
      
   else
      call runend('RayElementIntersection only ready for linear tetra')
   endif


end subroutine

subroutine RayTriangleIntersection(beam,coord,iflag,iIntersection)
   !According to Platis and Theoharis
   !http://graphics.di.uoa.gr/Downloads/papers/journals/p19.pdf
   use typre
   use Mod_Optics
   implicit none
   
   integer(ip), parameter :: ndime=3
   
   type(LightBeam) :: beam
   real(rp) :: coord(ndime,3)
   logical  :: iflag
   real(rp) :: iIntersection(ndime)
   
   integer(ip) :: iedge,ivertex
   
   integer(ip) :: edges(2,3)
   
   real(rp) :: pi_r(3,2),pi_e(3,2,3),epi_prods(3),upi_prods(3),sum_epi_prods
   
   edges(1,1) = 2
   edges(2,1) = 3
   
   edges(1,2) = 3
   edges(2,2) = 1
   
   edges(1,3) = 1
   edges(2,3) = 2
   
   iflag = .false.
   
   
   pi_r(:,1) = beam%direct
   call vecpro(beam%direct,beam%origin,pi_r(:,2),ndime)
   
   !pi for each edge
   do iedge = 1,3
      pi_e(:,1,iedge) = coord(1:ndime,edges(2,iedge))-coord(1:ndime,edges(1,iedge))
      call vecpro(pi_e(:,1,iedge),coord(1:ndime,edges(1,iedge)),pi_e(:,2,iedge),ndime)   
   enddo
   
   !pi_prod between pi ray and pi edges
   do iedge = 1,3
      call pi_prods(pi_r,pi_e(:,:,iedge),epi_prods(iedge))
   enddo
   
   !check wether there is aux_intersection or not
   if (minval(epi_prods) >= 0.0_rp .and. maxval (epi_prods) > 0.0_rp) iflag = .true.
   if (maxval(epi_prods) <= 0.0_rp .and. minval (epi_prods) < 0.0_rp) iflag = .true.
  
   if (minval(epi_prods) == 0.0_rp .and. maxval(epi_prods) == 0.0_rp) then
      write(*,*) 'warning, there where coplanar faces, should retest'
   endif
   
   !Compute aux_intersection point
   if (iflag .eqv. .true.) then
      sum_epi_prods = sum(epi_prods)
      upi_prods = epi_prods/sum_epi_prods  !equation 4 in platis theoharis
      iIntersection = 0.0_rp
      do ivertex = 1,3
         iIntersection = iIntersection + upi_prods(ivertex)*coord(:,ivertex)
      enddo
   endif





end subroutine


subroutine pi_prods(pi1,pi2,prod)
   use typre
   implicit none
   real(rp) :: pi1(3,2),pi2(3,2),prod
   
   prod = dot_product(pi1(:,1),pi2(:,2)) + dot_product(pi1(:,2),pi2(:,1))
end subroutine   
   
end module
   
