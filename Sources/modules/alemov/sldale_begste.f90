subroutine sldale_begste(a)
   use typre
   use Mod_Timer
   use Mod_Memor  
   use Mod_Mesh
   use Mod_Element
   use Mod_sldAlemov
   use def_parame
   use MPI
   implicit none
   class(sldAlemovProblem) :: a


end subroutine

subroutine sldale_endste(a,itask)
   use typre
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem) :: a
   integer(ip) :: itask,i

   !Update displacement
   do i = a%ncomp,2,-1
      a%Displacement(:,:,i) = a%Displacement(:,:,i-1)   
   enddo

   !Update velocities
   a%Velocity(:,:,2) = a%Velocity(:,:,1)

end subroutine sldale_endste
