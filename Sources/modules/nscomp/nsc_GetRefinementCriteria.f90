subroutine nsc_GetRefinementCriteria(a,markel)
   use typre
   use Mod_NSCompressible
   use Mod_Element

   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip) :: markel(*)
   
   real(rp) :: minz,maxz,miny,maxy,minx,maxx
   integer(ip) :: nelem,ielem
   class(FiniteElement), pointer :: e => NULL()
   
   call a%Timer%Refine%Tic

   call a%Mesh%GetNelem(nelem)
   
   if (a%RefinerErrorEstimator /= 'TEST ') then
        
      call a%SpecificNSCompRefinementCriteria(markel)
         
      if (a%RefinerErrorEstimator == 'INITI') then
         if (a%istep<=1) then
            markel(1:nelem) = 1
         else
            markel(1:nelem) = 0
         endif  
      endif 

   endif   
   
   call a%Timer%Refine%Toc
contains

   subroutine SelectRefine(ominx,omaxx,ominy,omaxy,ominz,omaxz)
      implicit none
      real(rp) :: ominx,omaxx,ominy,omaxy,ominz,omaxz
   
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
         minz = minval(e%elcod(3,1:e%pnode))
         maxz = maxval(e%elcod(3,1:e%pnode))
         
         if (maxx > ominx .and. minx < omaxx .and. maxy > ominy .and. miny < omaxy .and. maxz > ominz .and. minz < omaxz) then
            markel(ielem) = 1
         endif
      enddo elements
      call a%Mesh%ElementDeAlloc(e,a%Memor,'ForceClosedRule','vmass')
      
   end subroutine
   
   subroutine SelectUnRefine(ominx,omaxx,ominy,omaxy,ominz,omaxz)
      implicit none
      real(rp) :: ominx,omaxx,ominy,omaxy,ominz,omaxz
   
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
         minz = minval(e%elcod(3,1:e%pnode))
         maxz = maxval(e%elcod(3,1:e%pnode))
         
         if (maxx > ominx .and. minx < omaxx .and. maxy > ominy .and. miny < omaxy .and. maxz > ominz .and. minz < omaxz) then
            markel(ielem) = -1
         endif
      enddo elements
      call a%Mesh%ElementDeAlloc(e,a%Memor,'ForceClosedRule','vmass')
      
   end subroutine   
   

end subroutine
