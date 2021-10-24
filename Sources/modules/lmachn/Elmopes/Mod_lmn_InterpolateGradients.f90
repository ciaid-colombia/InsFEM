module Mod_lmn_InterpolateGradients
   use typre
   use Mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersInterpolateGradients
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersInterpolateGradients(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGradients)
         endif  
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
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
      call e%gradient(1,eltem,grtem)            !Temp. gradient

   end subroutine  
   
end module

