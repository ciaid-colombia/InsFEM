subroutine ale_iniunk(a)
   !-----------------------------------------------------------------------
   !> This routine initializes the last component of the unknown at the boundaries
   !! at the beginning of the time step so the selected integration scheme can be 
   !! used later on. 
   !-----------------------------------------------------------------------
   use typre
   use Mod_Alemov
   implicit none
   class(AlemovProblem) :: a
   integer(ip) :: icomp,ipoin,npoin,idime,ndime
   
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   
   do ipoin=1,npoin
      do idime=1,ndime
         if((a%kfl_fixno(idime,ipoin)==1) .or. (a%kfl_fixno(idime,ipoin)==0)) then
            a%Displacement(idime,ipoin,a%ncomp) = a%bvess(idime,ipoin,1)
         end if
      end do
   end do
   
end subroutine
