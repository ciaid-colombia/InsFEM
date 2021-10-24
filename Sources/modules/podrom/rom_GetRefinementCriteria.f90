subroutine rom_GetRefinementCriteria(a,markel)
   use typre
   use Mod_Element
   use Mod_ZZErrorEstimator
   use Mod_PodRom

   implicit none
   class(PodRomProblem)         :: a
   class(FiniteElement), pointer :: e => NULL()
   integer(ip)                   :: markel(*)
   integer(ip)                   :: ndime,nelem,idofr
   real(rp), allocatable         :: error(:),errorROM(:,:)

   call a%Problem%Timer%Refine%Tic
   
   call a%Mesh%GetNelem(nelem)
   call a%Memor%alloc(nelem,error,'error','rom_GetRefinementCriteria')
   call a%Memor%alloc(nelem,a%ndofr,errorROM,'errorROM','rom_GetRefinementCriteria')
   call a%Mesh%GetNdime(ndime)

   if (a%RefinerErrorEstimator == 'ZZ   ') then 
      do idofr = 1,a%ndofr
         call ZZErrorEstimator(a%Mesh,a%Memor,a%ndofn,a%Basis(:,:,idofr),errorROM(:,idofr))
      end do
      error = maxval(errorROM,2)
      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)
   elseif (a%RefinerErrorEstimator == 'GRADI') then
      do idofr = 1,a%ndofr
         call GradientErrorEstimator(a%Mesh,a%Memor,a%ndofn,a%Basis(:,:,1),errorROM(:,idofr))
      end do
      error = maxval(errorROM,2)
      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)
   elseif (a%RefinerErrorEstimator == 'PHYSI') then
      call a%Problem%GetRefinementCriteria(markel)
   endif   

   call a%Memor%dealloc(nelem,a%ndofr,errorROM,'errorROM','rom_GetRefinementCriteria')
   call a%Memor%dealloc(nelem,error,'error','rom_GetRefinementCriteria')

   call a%Problem%Timer%Refine%Toc

end subroutine
