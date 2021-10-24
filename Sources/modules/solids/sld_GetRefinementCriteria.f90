subroutine sld_GetRefinementCriteria(a,markel)
   use typre
   use Mod_Element
   use Mod_ZZErrorEstimator
   use Mod_Solids

   implicit none
   class(SolidsProblem)          :: a
   class(FiniteElement), pointer :: e => NULL()
   integer(ip)                   :: markel(*)
   real(rp), allocatable         :: error(:)
   real(rp) :: miny,maxy,minx,maxx,minz,maxz
   integer(ip) :: nelem,ielem,ndime

   call a%Timer%Refine%Tic

   if (a%RefinerErrorEstimator /= 'TEST ') then

      call a%Mesh%GetNelem(nelem)
      call a%Memor%alloc(nelem,error,'error','sld_GetRefinementCriteria')
      call a%Mesh%GetNdime(ndime)
      if (a%RefinerErrorEstimator == 'ZZ   ') then 
          if (a%kfl_timei == 1) then
              call ZZErrorEstimator(a%Mesh,a%Memor,ndime,a%accel(:,:,1),error)
          else
              call ZZErrorEstimator(a%Mesh,a%Memor,ndime,a%disp(:,:,1),error)
          endif
      elseif (a%RefinerErrorEstimator == 'GRADI') then
          if (a%kfl_timei == 1) then
              call GradientErrorEstimator(a%Mesh,a%Memor,ndime,a%accel(:,:,1),error)
          else
              call GradientErrorEstimator(a%Mesh,a%Memor,ndime,a%disp(:,:,1),error)
          endif
      elseif (a%RefinerErrorEstimator == 'SUBSC') then
         !call a%SpecificSubscalesRefCriteria(error)
         call runend('SLD_REFINEMENT_CRITERIA: subscales not yet applied to sld module')
      endif   

      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)
      call a%Memor%dealloc(nelem,error,'error','sld_GetRefinementCriteria')

      if (a%RefinerErrorEstimator == 'INITI') then
         if (a%istep<=1) then
            markel(1:nelem) = 1
         else
            markel(1:nelem) = 0
         endif  
      endif 
   else
      select case (mod(a%istep,4))
         case(0)
            call SelectRefine(0.5_rp,1.0_rp,0.5_rp,1.0_rp,0.0_rp,0.0_rp)
            
         case(1)
            call SelectRefine(0.5_rp,1.0_rp,0.5_rp,1.0_rp,0.0_rp,0.0_rp)
         
         case(2)
            !call SelectUnRefine(-1.0_rp,9.0_rp,-1.0_rp,9.0_rp,0.0_rp,0.0_rp)
            
         case(3)
            !call SelectUnRefine(-1.0_rp,9.0_rp,-1.0_rp,9.0_rp,0.0_rp,0.0_rp)   
            
         case default
      end select
   endif 

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
         
         if (maxx > ominx .and. minx < omaxx .and. maxy > ominy .and. miny < omaxy )then
            markel(ielem) = 1
         endif

      enddo elements
      call a%Mesh%ElementDeAlloc(e,a%Memor,'ForceClosedRule','vmass')
            markel(1) = 0
            markel(2) = 1
      
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
