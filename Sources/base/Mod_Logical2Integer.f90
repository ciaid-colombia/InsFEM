subroutine Logical2Integer(ndime,logicalarray,Integerarray)
   use typre
   implicit none
   integer(ip) :: ndime
   logical logicalarray(ndime)
   
   integer :: IntegerArray(ndime)
   
   integer(ip) :: idime
   
   do idime = 1,ndime
      if (logicalarray(idime) .eqv. .true.) then
         IntegerArray(idime) = 1_ip
      else
         IntegerArray(idime) = 0.0_ip
      endif
   enddo
end subroutine

subroutine Integer2Logical(ndime,integerarray,LogicalArray)
   use typre
   implicit none
   integer(ip) :: ndime
   logical LogicalArray(ndime)
   integer :: integerarray(ndime)
   
   integer(ip) :: idime
   
   do idime = 1,ndime
      if (integerarray(idime) == 0.0_ip) then
         LogicalArray(idime) = .false.
      else
         LogicalArray(idime) = .true.
      endif
   enddo
end subroutine



