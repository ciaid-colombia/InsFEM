subroutine lmn_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* lmachn/lmn_endite
   ! NAME 
   !    lmn_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   integer(ip) :: itask
   integer(ip) :: ndime
   real(rp)    :: rhsnorm
   
   interface
      subroutine lmn_EndElmope(LMProblem,task)
         use Mod_LowMach
         implicit none
         class(LowMachProblem), target :: LMProblem
         character(6) :: task
      end subroutine
   end interface
   
   select case(itask)

   case(0)

      call a%Mesh%GetNdime(ndime)
      a%veloc(:,:,1) = a%unkno(1:ndime,:)
      a%press(:,1)   = a%unkno(ndime+2,:)
      a%tempe(:,1)   = a%unkno(ndime+1,:)
   
   case(1)
      !Update the subscale and/or compute residual projection
      call lmn_EndElmope(a,'Endite')
      !Compute thermodynamic pressure and density field
      if (a%kfl_pther) call a%ComputeThermPressure

   case(2)

      !Update u(n,i-1,*) <-- u(n,i,*)
      a%veloc(:,:,2) = a%veloc(:,:,1)
      ! Update p(n,i-1,*) <-- p(n,i,*)
      a%press(:,2)   = a%press(:,1)
      ! Update T(n,i-1,*) <-- T(n,i,*)
      a%tempe(:,2)   = a%tempe(:,1)
      ! Update thermodynamic pressure
      a%pther(2)     = a%pther(1)


   end select  

end subroutine lmn_endite
