subroutine tem_ManufacturedBoundaryConditions(a)
   use typre
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
   
   integer(ip) :: npoin,ipoin,ibopo
   real(rp), pointer :: exnor(:,:) => NULL()
   
   if (a%ManufacturedBoundaryCondition == 1) then
      a%kfl_fixbo = -1
      a%kfl_fixno = -1
      a%bvess = 0
      
      call a%Mesh%GetNpoin(npoin)
      
      do ipoin = 1,npoin
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if (ibopo /= 0) then
            a%kfl_fixno(1,ipoin) = 1_ip
         endif
      enddo   

   endif


end subroutine
