subroutine rom_outerr(a)
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   
!   if (a%MPIsize > 1) then
!      call runend('PodRom not ready for parallel computations: at least, a parallel implementation of the Singular value decomposition is required')
!   endif 
   
   !if (a%kfl_itask == 1 .and. a%kfl_ROSNonLinear == 1) then
   !   call runend('kfl_ROSNonlinear not ready for run')
   !endif
   
end subroutine
