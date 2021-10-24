subroutine lmn_CrankNicolsonEndste(a)
   use typre
   use def_parame
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   
   a%veloc(:,:,1) = 2.0_rp*a%veloc(:,:,1)-a%veloc(:,:,3)
   a%tempe(:,1) = 2.0_rp*a%tempe(:,1)-a%tempe(:,3)
   a%press(:,1) = 2.0_rp*a%press(:,1)-a%press(:,3)
   a%pther(1) = 2.0_rp*a%pther(1)-a%pther(3)
   
end subroutine
