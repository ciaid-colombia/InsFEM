subroutine tem_elmvel(e,kfl_advec_tem,gpvel)
   !This subroutine computes an analytical velocity field for the temperature problem
   use typre
   use Mod_Element
   use def_parame
   implicit none
   class(FiniteElement) :: e
   integer(ip)         :: kfl_advec_tem
   real(rp)            :: gpvel(e%ndime)
   
   select case(kfl_advec_tem)
   case(1)
      !This one is never used since it corresponds to external velocity (NSTINC or another)

   case(2)
      gpvel(1)=1.0_rp
      gpvel(2)=0.0_rp
   case(3)
      gpvel(1)=cos(pi/3_rp)
      gpvel(2)=sin(pi/3_rp)
   case(4)
      gpvel(1)=cos(pi/3_rp)*0.00001_rp
      gpvel(2)=sin(pi/3_rp)*0.00001_rp  
   case(6)
      gpvel(1)=5.0_rp
      gpvel(2)=-9.0_rp
   case(7)
      gpvel(1) = 0.0_rp
      gpvel(2) = 0.0_rp
      gpvel(3) = 1.0_rp
   case(8)
      gpvel(1) = 0.0_rp
      gpvel(2) = 0.0_rp
      gpvel(3) = -1.0_rp
   case(9)
      !This one corresponds to burgers equation, gpvel is set externally
   end select

end subroutine tem_elmvel
