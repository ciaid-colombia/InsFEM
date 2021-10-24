subroutine nsi_GetRefinementCriteria(a,markel)
   use typre
   use Mod_Element
   use Mod_ZZErrorEstimator
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip) :: markel(*)
   
   real(rp) :: miny,maxy,minx,maxx,minz,maxz
   real(rp) :: TotalEstimatedError
   integer(ip) :: nelem,ielem,ndime
   class(FiniteElement), pointer :: e => NULL()
   
   real(rp), allocatable :: error(:)
  
   call a%Timer%Refine%Tic
   
   call a%Mesh%GetNelem(nelem)
   call a%Memor%alloc(nelem,error,'error','nsi_GetRefinementCriteria')
   call a%Mesh%GetNdime(ndime)
   
   if (a%RefinerErrorEstimator /= 'TEST ') then
   
      if (a%RefinerErrorEstimator == 'ZZ   ') then 
         call ZZErrorEstimator(a%Mesh,a%Memor,ndime,a%veloc(:,:,1),error)
      elseif (a%RefinerErrorEstimator == 'GRADI') then
         call GradientErrorEstimator(a%Mesh,a%Memor,ndime,a%veloc(:,:,1),error)
      elseif (a%RefinerErrorEstimator == 'SUBSC') then
         call a%SpecificSubscalesRefCriteria(error,TotalEstimatedError)
      endif   
      
      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)   
      
      if (a%RefinerErrorEstimator == 'INITI') then
         if (a%istep==1) then
            call SelectRefine(-0.0055_rp,0.0055_rp,0.0075_rp,0.025_rp,-0.02_rp,-0.003_rp)
         elseif (a%istep==2 .or. a%istep==3) then
            call SelectRefine(-0.0055_rp,0.0055_rp,0.0177_rp,0.025_rp,-0.02_rp,-0.011_rp)         
         else
            markel(1:nelem) = 0
         endif  
      endif 
      
   else
      select case (mod(a%istep,4))
         case(0)
            call SelectRefine(0.0_rp,8.0_rp,0.0_rp,8.0_rp,0.0_rp,0.0_rp)
            
         case(1)
            call SelectRefine(0.0_rp,8.0_rp,0.0_rp,8.0_rp,0.0_rp,0.0_rp)   
         
         case(2)
            call SelectUnRefine(-1.0_rp,9.0_rp,-1.0_rp,9.0_rp,0.0_rp,0.0_rp)
            
         case(3)
            call SelectUnRefine(-1.0_rp,9.0_rp,-1.0_rp,9.0_rp,0.0_rp,0.0_rp)   
            
         case default
      end select
   endif 
   
   call a%Memor%dealloc(nelem,error,'error','nsi_GetRefinementCriteria')
   
   call a%Timer%Refine%Toc
   
contains

   subroutine SelectRefine(ominx,omaxx,ominy,omaxy,ominz,omaxz)
      implicit none
      real(rp) :: ominx,omaxx,ominy,omaxy,ominz,omaxz
      integer(ip) :: inode
   
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
