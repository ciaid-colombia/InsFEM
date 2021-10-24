subroutine lev_surface(a,e,ielem,ndime,nipoin,dsurf,points,weipg)
!DESCRIPTION
!   This get the intersection surface of the level set and the coordinates
!   of the cut points
!-----------------------------------------------------------------------
   use typre
   use Mod_Mesh
   use Mod_Memor   
   use Mod_LevelSet
   use Mod_Element   
   use def_parame
   implicit none
   class(LevelSetProblem) :: a   
   class(FiniteElement) :: e
   integer(ip), intent(in)      :: ielem,ndime
   real(rp), intent(out)        :: dsurf
   integer(ip), intent(out)     :: nipoin
   
   integer(ip)                  :: inode,jnode,lpoin,hpoin,idime,ipoin
   real(rp)                     :: iratio(4),points(ndime,4),weipg(4)
   real(rp)                     :: t1(ndime),t2(ndime),cprod(ndime),rdsurf
   real(rp), pointer            :: coord1(:) => NULL(),coord2(:) => NULL()
   
   
   !Initialization
   weipg=0
   !To obtain the intersection points
   nipoin = 0 !number of intersection points
   
   if(ndime==2)then
   
      if(e%pnode==3)then   
         do inode = 1,e%pnode
            do jnode = 1+inode,e%pnode
            
            lpoin=e%lnods(inode)
            hpoin=e%lnods(jnode)
            
               if((((a%level(lpoin,3) <= 0.0_rp) .and. (a%level(hpoin,3) > 0.0_rp))).or. &
               (((a%level(lpoin,3) > 0.0_rp) .and. (a%level(hpoin,3) <= 0.0_rp))))then
                           
                  nipoin = nipoin + 1
                        
                  !ratio to define the intersection points
                  iratio(nipoin) = a%level(lpoin,3)/(a%level(lpoin,3) - a%level(hpoin,3))
                        
                  call a%Mesh%GetPointCoord(lpoin,coord1)
                  call a%Mesh%GetPointCoord(hpoin,coord2)                  
                  
                  points(:,nipoin) = coord1(:) + (coord2(:) - coord1(:))*iratio(nipoin)
                        
               end if
                  
            end do
         end do
         
      elseif(e%pnode==4)then
      
         if((((a%level(e%lnods(1),3) <= 0.0_rp) .and. (a%level(e%lnods(2),3) > 0.0_rp))).or. &
           (((a%level(e%lnods(1),3) > 0.0_rp) .and. (a%level(e%lnods(2),3) <= 0.0_rp))))then
           
            lpoin=e%lnods(1)
            hpoin=e%lnods(2)
                           
            nipoin = nipoin + 1
                        
            !ratio to define the intersection points
            iratio(nipoin) = a%level(lpoin,3)/(a%level(lpoin,3) - a%level(hpoin,3))
                        
            call a%Mesh%GetPointCoord(lpoin,coord1)
            call a%Mesh%GetPointCoord(hpoin,coord2)                  
                  
            points(:,nipoin) = coord1(:) + (coord2(:) - coord1(:))*iratio(nipoin)
                        
         end if  
         
         if((((a%level(e%lnods(1),3) <= 0.0_rp) .and. (a%level(e%lnods(4),3) > 0.0_rp))).or. &
           (((a%level(e%lnods(1),3) > 0.0_rp) .and. (a%level(e%lnods(4),3) <= 0.0_rp))))then
           
            lpoin=e%lnods(1)
            hpoin=e%lnods(4)
                           
            nipoin = nipoin + 1
                        
            !ratio to define the intersection points
            iratio(nipoin) = a%level(lpoin,3)/(a%level(lpoin,3) - a%level(hpoin,3))
                        
            call a%Mesh%GetPointCoord(lpoin,coord1)
            call a%Mesh%GetPointCoord(hpoin,coord2)                  
                  
            points(:,nipoin) = coord1(:) + (coord2(:) - coord1(:))*iratio(nipoin)
                        
         end if    
         
         if((((a%level(e%lnods(2),3) <= 0.0_rp) .and. (a%level(e%lnods(3),3) > 0.0_rp))).or. &
           (((a%level(e%lnods(2),3) > 0.0_rp) .and. (a%level(e%lnods(3),3) <= 0.0_rp))))then
           
            lpoin=e%lnods(2)
            hpoin=e%lnods(3)
                           
            nipoin = nipoin + 1
                        
            !ratio to define the intersection points
            iratio(nipoin) = a%level(lpoin,3)/(a%level(lpoin,3) - a%level(hpoin,3))
                        
            call a%Mesh%GetPointCoord(lpoin,coord1)
            call a%Mesh%GetPointCoord(hpoin,coord2)                  
                  
            points(:,nipoin) = coord1(:) + (coord2(:) - coord1(:))*iratio(nipoin)
                        
         end if 
         
         if((((a%level(e%lnods(3),3) <= 0.0_rp) .and. (a%level(e%lnods(4),3) > 0.0_rp))).or. &
           (((a%level(e%lnods(3),3) > 0.0_rp) .and. (a%level(e%lnods(4),3) <= 0.0_rp))))then
           
            lpoin=e%lnods(3)
            hpoin=e%lnods(4)
                           
            nipoin = nipoin + 1
                        
            !ratio to define the intersection points
            iratio(nipoin) = a%level(lpoin,3)/(a%level(lpoin,3) - a%level(hpoin,3))
                        
            call a%Mesh%GetPointCoord(lpoin,coord1)
            call a%Mesh%GetPointCoord(hpoin,coord2)                  
                  
            points(:,nipoin) = coord1(:) + (coord2(:) - coord1(:))*iratio(nipoin)
                        
         end if        
      
      else
         call runend('lev_Surface : level set only for linear elements')
      end if
      
      if(nipoin == 2)then !linear triangle or Cuadrilateral
      
         dsurf=0.0_rp
                  
         do idime = 1,ndime
            dsurf = dsurf + (points(idime,1) - points(idime,2))**2.0_rp    
         end do
               
         dsurf=sqrt(dsurf)
         
                  
         weipg(1) = 0.5_rp
         weipg(2) = 0.5_rp               
                              
      else
         write(*,*) 'e%pnode nipoin ielem  '
         write(*,*) e%pnode, nipoin, ielem            
         call runend('lev_Surface: wrong number of intersection points')
      endif
      
      
   elseif(ndime==3)then

      if(e%pnode==4)then   
         do inode = 1,e%pnode
            do jnode = 1+inode,e%pnode
            
            lpoin=e%lnods(inode)
            hpoin=e%lnods(jnode)
            
               if((((a%level(lpoin,3) <= 0.0_rp) .and. (a%level(hpoin,3) > 0.0_rp))).or. &
               (((a%level(lpoin,3) > 0.0_rp) .and. (a%level(hpoin,3) <= 0.0_rp))))then
                           
                  nipoin = nipoin + 1
                        
                  !ratio to define the intersection points
                  iratio(nipoin) = a%level(lpoin,3)/(a%level(lpoin,3) - a%level(hpoin,3))
                        
                  call a%Mesh%GetPointCoord(lpoin,coord1)
                  call a%Mesh%GetPointCoord(hpoin,coord2)                  
                  
                  points(:,nipoin) = coord1(:) + (coord2(:) - coord1(:))*iratio(nipoin)
                        
               end if
                  
            end do
         end do  
      else
         call runend('lev_Surface: level set only for linear tetrahedra in 3d ')
      end if
      
      
      if(nipoin == 3)then !linear tetrahedra
         
         dsurf=0.0_rp
         
         t1 = points(:,2) - points(:,1)
         t2 = points(:,3) - points(:,1)
         call vecpro(t1,t2,cprod,3)
         call vecnor(cprod,ndime,dsurf,2) 
         dsurf = dsurf*0.5_rp
         
               
         weipg(1) = 1.0_rp/3.0_rp
         weipg(2) = 1.0_rp/3.0_rp
         weipg(3) = 1.0_rp/3.0_rp
         
      elseif(nipoin==4)then

         rdsurf=0.0_rp
         dsurf=0.0_rp

         t1=points(:,1)-points(:,2)
         t2=points(:,1)-points(:,3)
         call vecpro(t1,t2,cprod,3)
         call vecnor(cprod,3,rdsurf,2)         
         
         dsurf= dsurf + abs(rdsurf)*0.5_rp         
         
         rdsurf=0.0_rp
         t1=points(:,1)-points(:,4)
         t2=points(:,1)-points(:,3)
         call vecpro(t1,t2,cprod,3)
         call vecnor(cprod,3,rdsurf,2)
                  
         dsurf= dsurf + abs(rdsurf)*0.5_rp         
         
         weipg(1) = 1.0_rp/4.0_rp
         weipg(2) = 1.0_rp/4.0_rp
         weipg(3) = 1.0_rp/4.0_rp
         weipg(4) = 1.0_rp/4.0_rp      

      else
         write(*,*) nipoin,ielem
         call runend('lev_Surface : wrong number of intersection')
         
      end if         
    
   endif


end subroutine
