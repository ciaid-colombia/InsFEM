 module Mod_nsm_InterpolateGradients
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersInterpolateGradients
   
   type, extends(PointerSetter) :: SPInterpolateGradients
contains
      procedure :: SpecificSet => SpecificSetInterpolateGradients
   end type
   type(SPInterpolateGradients) :: SetPointersInterpolateGradients
   
contains
   
   subroutine SpecificSetInterpolateGradients(d)
      implicit none
      class(SPInterpolateGradients) :: d
         
      call ConcatenateProcedures(ProcHook_Interpolates,InterpolateGradients)
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine InterpolateGradients
      implicit none
      integer(ip) :: idime
      
      call e%gradient(e%ndime,elvel(:,:,1),grvel)    !Vel. gradient
      !Velocity divergence
      divvel = 0.0_rp
      do idime = 1,e%ndime
         divvel = divvel + grvel(idime,idime)
      enddo
      call e%gradient(1,elpre,grpre)            !Press. gradient
   end subroutine  
   
end module

