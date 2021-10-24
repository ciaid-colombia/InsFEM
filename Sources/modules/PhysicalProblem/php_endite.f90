subroutine php_endite(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the internal loop iteration
   !-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: itask

   
   call a%Timer%Endite%Tic
   
   select case(itask)

   case(0)

      !Things to be done before the convergence check
      call a%SpecificEndite(10)   
         
      !Compute convergence residual of the internal iteration 
      call a%Cvgunk(one)
   
      call a%SpecificEndite(zero)   

   case(1)

      call a%SpecificEndite(one)

   case(2)

      !Compute convergence residual of the external iteration
      call a%Cvgunk(two)

      call a%SpecificEndite(two)

   case(4)

      !Compute convergence residual of case coupling iteration
      call a%Cvgunk(4)

      call a%SpecificEndite(4)

   end select
   
   call a%Timer%Endite%Toc

end subroutine php_endite
