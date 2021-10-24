subroutine tem_endste(a,itask)
   use typre
   use Mod_Temperature
   use def_parame
   implicit none
   class(TemperatureProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: ielem,nelem,icomp
   
   interface
      subroutine tem_EndsteElmope(TempeProblem) 
         use Mod_Temperature
         implicit none
         class(TemperatureProblem), target :: TempeProblem
      end subroutine
   end interface
   
   !Pre-CrankNicolson
   if (itask == 1) then
   
      !Update the subscales
      !Compute dissipation 
      !and other endste elemental computations using u^(n+theta)
      call tem_EndsteElmope(a)
      
   !Post-CrankNicolson   
   elseif (itask == 2) then
   
      !Update a%velocity and a%pressure
      if (a%ncomp > 3) then
          do icomp = a%ncomp, 4,-1
            a%tempe(:,icomp) = a%tempe(:,icomp-1)
         enddo
      endif
      a%tempe(:,3) = a%tempe(:,1)  ! Vn = Vn+1
      

     !Boundary operation at the end
     if(a%kfl_outfm == 1)then
         call a%EndBouope
     end if
   
     !Update the subgrid velocity
     if(a%kfl_trasg/=0) then
        call a%Mesh%GetNelem(nelem)
       
        !Update the subscales
        !CrankNicolson
        if (a%kfl_tsche_1st_current == 'CN   ' .or. a%kfl_tsche_1st_current == 'CNOBS') then   
           !CrankNicolson
           do ielem=1,nelem
             a%tesgs(ielem)%a(2,:)=2.0_rp*a%tesgs(ielem)%a(1,:)-a%tesgs(ielem)%a(2,:)
           end do
        else
           do ielem=1,nelem
              a%tesgs(ielem)%a(2,:)=a%tesgs(ielem)%a(1,:)
           end do
        end if
     end if

  end if
   
end subroutine
