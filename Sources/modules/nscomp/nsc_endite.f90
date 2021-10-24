subroutine nsc_endite(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   
   use def_parame
   use Mod_NSCompressible

   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip) :: itask
   

   
   select case(itask)

   case(0)

      call a%SpecificNSCompEndite(0)

   case(1)

      call a%SpecificNSCompEndite(1)

   case(2)

      call a%SpecificNSCompEndite(2)

   end select  

end subroutine nsc_endite
