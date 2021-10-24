subroutine lmn_endste(a,itask)
   use typre
   use def_parame
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a

   integer(ip) :: itask
   integer(ip) :: ielem, nelem, icomp

   interface
      subroutine lmn_EndElmope(LMProblem,task)
         use Mod_LowMach
         implicit none
         class(LowMachProblem), target :: LMProblem
         character(6) :: task
      end subroutine
   end interface

   !Pre-CrankNicolson
   if (itask == 1) then
      !Compute dissipation and other endste elemental computations using u^(n+theta)
      !Also compute the dynamic subgrid scales
      call lmn_EndElmope(a,'Endste')

      !Post-CrankNicolson
   elseif (itask == 2) then

      ! Update a%velocity, a%pressure and a%temperature
      !Higher order components
      if (a%ncomp > 3) then
         do icomp = a%ncomp, 4,-1
            a%veloc(:,:,icomp) = a%veloc(:,:,icomp-1)  ! Vn-1 = Vn
            a%press(:,icomp)   = a%press(:,icomp-1)
            a%tempe(:,icomp)   = a%tempe(:,icomp-1)
            a%pther(icomp)     = a%pther(icomp-1)
         end do
      end if
      !Previous time step component
      a%veloc(:,:,3) = a%veloc(:,:,1)  ! Vn = Vn+1
      a%press(:,3)   = a%press(:,1)
      a%tempe(:,3)   = a%tempe(:,1)
      a%pther(3)     = a%pther(1)

      !Boundary operation at the end
      if(a%kfl_outfm == 1)then
         call a%EndBouope
      end if

      ! Update the subgrid scales
      call a%Mesh%GetNelem(nelem)
      if (a%kfl_tacsg>0) then
         !Crank Nicolson schemes
         if (a%kfl_tsche_1st_current == 'CN   ') then   
            do ielem=1,nelem
               a%vesgs(ielem)%a(:,2,:)=2.0_rp*a%vesgs(ielem)%a(:,1,:)-a%vesgs(ielem)%a(:,2,:)
               a%tesgs(ielem)%a(2,:)=2.0_rp*a%tesgs(ielem)%a(1,:)-a%tesgs(ielem)%a(2,:)
            end do
         else
            do icomp = a%ncomp-1, 2,-1
               do ielem=1,nelem
                  a%vesgs(ielem)%a(:,icomp,:)=a%vesgs(ielem)%a(:,icomp-1,:)
                  a%tesgs(ielem)%a(icomp,:)=a%tesgs(ielem)%a(icomp-1,:)
               end do
            end do
         end if
      end if
   endif

end subroutine lmn_endste
