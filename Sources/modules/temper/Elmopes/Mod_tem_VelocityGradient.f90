module Mod_tem_VelocityGradient
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_Advection
   implicit none
  private
  public SetPointersVelocityGradient
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersVelocityGradient(itask)
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
         
            !We need the advection velocity
            call SetPointersAdvectionVelocity(1)
            
            call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGrvel)
            call ConcatenateProcedures(ProcHook%Initializations,AllocGrvel)
            call ConcatenateProcedures(ProcHook%Finalizations,DeAllocGrvel)
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !--------------------------------------------------------------------
   !Actual computations
   
  
   
   !Interpolation of the velocity gradient
   subroutine AllocGrvel
      implicit none
      
      call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','tem_elmope')    !Vel. gradient
   end subroutine   
   
   subroutine DeAllocGrvel
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','tem_elmope')    !Vel. gradient
   end subroutine  
   
   subroutine InterpolateGrvel
      implicit none
      
      call e%gradient(e%ndime,elvel(:,:),grvel)    !Vel. gradient
   end subroutine   
end module