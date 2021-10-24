subroutine tem_GetRefinementCriteria(a,markel)
   use typre
   use Mod_Temperature
   use Mod_Element
   use Mod_ZZErrorEstimator
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: markel(*)
   
   real(rp) :: miny,maxy,minx,maxx
   real(rp) :: TotalEstimatedError
   integer(ip) :: nelem,ielem
   class(FiniteElement), pointer :: e => NULL()
   
   real(rp), allocatable :: error(:)
   
   call a%Timer%Refine%Tic
   
   if (a%RefinerErrorEstimator /= 'TEST ') then
      
      call a%Mesh%GetNelem(nelem)
      call a%Memor%alloc(nelem,error,'error','tem_GetRefinementCriteria')
      if (a%RefinerErrorEstimator == 'ZZ   ') then 
         call ZZErrorEstimator(a%Mesh,a%Memor,1_ip,a%tempe,error)
      elseif (a%RefinerErrorEstimator == 'GRADI') then
         call GradientErrorEstimator(a%Mesh,a%Memor,1_ip,a%tempe,error)
      elseif (a%RefinerErrorEstimator == 'SUBSC') then
         call a%SpecificSubscalesRefCriteria(error,TotalEstimatedError)
      endif   
      
      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)
      call a%Memor%dealloc(nelem,error,'error','tem_GetRefinementCriteria')
      
   elseif (a%RefinerErrorEstimator == 'TEST ') then
      select case (a%istep)
         case(0)
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
         
         case(1)
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
         
         case(2)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(3)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(4)   
            call SelectRefine(0.6_rp,1.2_rp,0.6_rp,1.2_rp)
            
         case(5)   
            call SelectRefine(0.6_rp,1.2_rp,0.6_rp,1.2_rp)
            
         case(6)
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(7)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(8)   
            call SelectRefine(0.6_rp,1.2_rp,0.75_rp,0.9_rp)
            
         case(9)   
            call SelectRefine(0.6_rp,1.2_rp,0.75_rp,0.9_rp)
            
         case(10)   
            call SelectRefine(0.6_rp,1.2_rp,0.75_rp,0.9_rp)
            
         case(11)   
            call SelectRefine(0.3_rp,0.7_rp,0.3_rp,0.7_rp)
            
         case(12)   
            call SelectUnRefine(0.3_rp,1.7_rp,0.3_rp,1.7_rp)
            
         case(13)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(14)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(15)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(16)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(17)   
            call SelectRefine(0.5_rp,1.2_rp,0.5_rp,1.9_rp)
            
         case(18)   
            call SelectRefine(0.5_rp,1.2_rp,0.5_rp,1.9_rp)
            
         case(19)   
            call SelectRefine(0.5_rp,1.2_rp,0.5_rp,1.9_rp)
            
         case(20)   
            call SelectRefine(0.5_rp,1.2_rp,0.5_rp,1.9_rp)
            
         case(21)
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(22)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(23)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(24)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case default
      end select
   elseif (a%RefinerErrorEstimator == 'TEST2') then
      select case (a%istep)
         case(0)
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
         
         case(1)
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
         
         case(2)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(3)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(4)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(5)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(6)
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(7)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(8)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(9)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(10)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(11)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(12)   
            !call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(13)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(14)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(15)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(16)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(17)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(18)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(19)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(20)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(21)
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(22)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(23)   
            call SelectUnRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case(24)   
            call SelectRefine(0.0_rp,1.0_rp,0.0_rp,1.0_rp)
            
         case default
      end select

   
   elseif (a%RefinerErrorEstimator == 'TEST3') then
      
      call SelectRefine(1.0_rp-1/(2.0_rp**a%istep),1.0_rp,1.0_rp-1/(2.0_rp**a%istep),1.0_rp)
   endif   
   
   call a%Timer%Refine%Toc
contains

   subroutine SelectRefine(ominx,omaxx,ominy,omaxy)
      implicit none
      real(rp) :: ominx,omaxx,ominy,omaxy
   
      call a%Mesh%GetNelem(nelem)
      markel(1:nelem) = 0

      
      call a%Mesh%ElementAlloc(e,a%Memor,'ForceClosedRule','vmass')
      
      elements : do ielem = 1,nelem
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)    
         
         miny = minval(e%elcod(2,1:e%pnode))
         maxy = maxval(e%elcod(2,1:e%pnode))
         minx = minval(e%elcod(1,1:e%pnode))
         maxx = maxval(e%elcod(1,1:e%pnode))
         
         if (maxx > ominx .and. minx < omaxx .and. maxy > ominy .and. miny < omaxy) then
            markel(ielem) = 1
         endif
      enddo elements
      call a%Mesh%ElementDeAlloc(e,a%Memor,'ForceClosedRule','vmass')
      
   end subroutine
   
   subroutine SelectUnRefine(ominx,omaxx,ominy,omaxy)
      implicit none
      real(rp) :: ominx,omaxx,ominy,omaxy
   
      call a%Mesh%GetNelem(nelem)
      markel(1:nelem) = 0
      
      call a%Mesh%ElementAlloc(e,a%Memor,'ForceClosedRule','vmass')
      
      elements : do ielem = 1,nelem
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)    
         
         miny = minval(e%elcod(2,1:e%pnode))
         maxy = maxval(e%elcod(2,1:e%pnode))
         minx = minval(e%elcod(1,1:e%pnode))
         maxx = maxval(e%elcod(1,1:e%pnode))
         
         if (maxx > ominx .and. minx < omaxx .and. maxy > ominy .and. miny < omaxy) then
            markel(ielem) = -1
         endif
      enddo elements
      call a%Mesh%ElementDeAlloc(e,a%Memor,'ForceClosedRule','vmass')
      
   end subroutine
   

end subroutine
