subroutine nsc_endste(a,itask)
   !This routine ends a time step of the compressible NS equations.
   use typre
   use def_parame
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip) :: itask

   if (itask == 1) then

      call a%SpecificNSCompEndste(1)
      
      call a%CalculateStatistics
         
   elseif (itask == 2) then

      call a%SpecificNSCompEndste(2)

      !Boundary operations at endste
      if((a%kfl_outfm == 1) .or. (a%kfl_nscnbc == 1))then
         call a%EndBouope
      end if

   endif

end subroutine nsc_endste
