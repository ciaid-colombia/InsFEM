subroutine lmn_ComputeDensity(a)
   use typre
   use Mod_Element
   use Mod_LowMach
   implicit none
   class(LowMachProblem)  :: a
   real(rp)               :: aux1
   integer(ip)            :: npoin,ipoin

   !Initializations
   call a%Mesh%GetNpoin(npoin)

   select case (a%kfl_eqnst)
   case (1)     !Ideal gas
      aux1 = a%pther(1)/a%sgasc
      do ipoin=1,npoin
         a%densf(1:ipoin) = aux1/a%tempe(1:ipoin,1)
      end do
   case default !Constant density
   end select
end subroutine   

