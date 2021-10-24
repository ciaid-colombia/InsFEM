subroutine tem_endite(a,itask)
   use typre
   use Mod_Temperature
   use Mod_SmoothedFieldGradient
   use def_parame
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: ndime
   
   select case(itask)

   case(0)

      !Update unknowns
      call a%Mesh%GetNdime(ndime)
      a%tempe(:,1) = a%unkno(1,:)
      
   case(1)

      !Update the subscale and/or compute residual projection
      call a%EnditeElmope

   case(2)

      !Update u(n,i-1,*) <-- u(n,i,*)
      a%tempe(:,2) = a%tempe(:,1)
   
   end select  
   
end subroutine
