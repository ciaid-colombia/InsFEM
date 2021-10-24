subroutine ProjectionTest
   use typre
   use Mod_Element
   use Mod_Memor
   use def_domain
   use Mod_Postpr
   use def_master
   implicit none
   
   integer(ip) :: ndime,npoin,nelem,ipoin,idime,igaus,ielem
   real(rp), allocatable :: array(:,:)
   
   real(rp) :: gpcod(3), gpval,dvol
   
   type(MemoryMan) :: Memor
   class(FiniteElement), pointer :: e => NULL()
   
   call Mesh%GetNdime(ndime)
   call Mesh%GetNpoin(npoin)
   call Mesh%GetNelem(nelem)
   
   call Mesh%ElementAlloc(e,Memor,'DefaultRule','nsm_EnditeElmope')
   
   allocate(array(3,npoin))
   
   do iiter = 1,3
   
      array = 0.0_rp
      
      elements: do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)    
         
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         
         !Element length at center of gravity
         call e%elmlen

         !Gauss Point Loop
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            
            dvol = e%weigp(e%igaus)*e%detjm
            
            call e%interpg(e%ndime,e%elcod,gpcod)
            
            gpval = 12*sin(gpcod(1))*cos(gpcod(2))
            
            array(1,e%lnods(1:e%pnode)) = array(1,e%lnods(1:e%pnode)) + gpval*dvol*e%shape(1:e%pnode,e%igaus)
            array(2,e%lnods(1:e%pnode)) = array(2,e%lnods(1:e%pnode)) + 5*gpval*dvol*e%shape(1:e%pnode,e%igaus)
            array(3,e%lnods(1:e%pnode)) = array(3,e%lnods(1:e%pnode)) + 12*gpval*dvol*e%shape(1:e%pnode,e%igaus)
            
            
         enddo gauss_points
         
         
         
         
         
         
      enddo elements
      !enddo
      
      call Mesh%Project(3,array)

      !call Mesh%Smooth(1,array)
      
   enddo
   


   call FilePostpr%postpr(array,'array',1_ip,2.0_rp,Mesh)
   
   stop
   




end subroutine

   