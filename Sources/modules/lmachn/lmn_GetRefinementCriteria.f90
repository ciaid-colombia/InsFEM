subroutine lmn_GetRefinementCriteria(a,markel)
   use typre
   use Mod_ZZErrorEstimator
   use Mod_Element
   use Mod_LowMach

   implicit none
   class(LowMachProblem)         :: a
   class(FiniteElement), pointer :: e => NULL()
   real(rp)                      :: TotalEstimatedError
   integer(ip)                   :: markel(*)
   integer(ip)                   :: ndime,nelem
   real(rp), allocatable         :: error(:)

   call a%Timer%Refine%Tic

   call a%Mesh%GetNelem(nelem)
   call a%Memor%alloc(nelem,error,'error','lmn_GetRefinementCriteria')
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
         if (a%istep<=1) then
            markel(1:nelem) = 1
         else
            markel(1:nelem) = 0
         endif  
      endif 
   endif   
      
   call a%Memor%dealloc(nelem,error,'error','lmn_GetRefinementCriteria')
   call a%Timer%Refine%Toc

end subroutine
