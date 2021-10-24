subroutine php_Outerr(a, itask)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: kfl_isManufactured
   logical :: isAdaptive
   integer(ip) :: ibopo, idofn,ipoin,npoin
   logical :: SomethingisWrong = .false.

   if (itask == 1) then
      call a%Mesh%IsManufactured(kfl_isManufactured)

      if (a%kfl_BoundaryConditionsReadStrategy /= 1 .and. kfl_isManufactured == 1) then
         call runend('Boundary Conditions for ', a%namod, ' need to be Manufactured for a Manufactured Mesh')
      endif
      
      !Specific Errors
      call a%SpecificOuterr
   
   elseif (itask == 2) then
      
      !If Adaptive we need to check that all the fixno
      !are on boundary nodes
      !Transmission of fixno info to child nodes does not work if 
      !fixno is on interior nodes
      call a%Mesh%IsAnAdaptiveMesh(isAdaptive)
      if (isAdaptive) then
         call a%Mesh%GetNpoin(npoin)
         do ipoin = 1,npoin
            call a%Mesh%GetIbopo(ipoin,ibopo)
            !Only for interior points
            if (ibopo == 0) then
               do idofn = 1,a%ndofbc
                  if (a%kfl_fixno(idofn,ipoin) > 0) then
                     SomethingisWrong = .true.
                     
                  endif
               enddo
            endif
         enddo
         if (SomethingisWrong) write(*,*) 'Interior point with a fixno, this is not compatible with adaptive'
      endif
   
   
   endif
   
end subroutine


