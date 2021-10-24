subroutine nsc_RestarMaster(a,itask)
   use typre
   use def_parame
   use Mod_PhysicalProblem
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip), intent(in)    :: itask

   interface
      !subroutine nsc_nsirestar(a)
      !   use typre
      !   use Mod_NSCompressible
      !   implicit none
      !   class(NSCompressibleProblem) :: a
      !end subroutine

      subroutine php_restar(a,itask)
         use typre
         use Mod_PhysicalProblem
         implicit none
         class(PhysicalProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine
   end interface
   
   select case (itask)
   
   case (1)
   
      if (a%kfl_nsirstar == 1) then
         if (a%kfl_inter == 0) then
            call runend ('Restart from incompressible flow requires using an interpolated restart')
         else
            !call nsc_nsirestar(a)
         endif
      else
         call php_restar(a,itask)
      endif   
         
   case (2)

      call php_restar(a,itask)

   end select

end subroutine nsc_RestarMaster
