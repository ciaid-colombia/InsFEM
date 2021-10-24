 module Mod_supm_InterpolateGradients
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersInterpolateGradients
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
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
            call ConcatenateProcedures(ProcHook_Interpolates,InterpolateGradients)
         endif  
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------
!SUBROUTINES

   subroutine InterpolateGradients
      implicit none
      call e%gradient(1   ,elpre(:,1)  ,grpre)       !Press. gradient
      call e%gradient(e%ndime,elvel(:,:,1),grvel)    !Vel. gradient
      call e%gradient(auxtens,elsig(:,:,1),grsig)    !Sig. gradient     
   end subroutine
   
end module
