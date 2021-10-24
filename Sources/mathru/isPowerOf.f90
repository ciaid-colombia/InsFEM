subroutine isPowerOf(number,power,kfl_isPowerOf)
   use typre
   implicit none
   integer(ip) :: number,power
   logical :: kfl_isPowerOf
   
   integer(ip) :: auxnumber
   
   auxnumber = number
   do while (mod(auxnumber,power) == 0) 
      auxnumber = auxnumber/power
   enddo
   if (auxnumber > 1) then
      kfl_isPowerOf = .false.
   else
      kfl_isPowerOf = .true.
   endif
   
   
end subroutine