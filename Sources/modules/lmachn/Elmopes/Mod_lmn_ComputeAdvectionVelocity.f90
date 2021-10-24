module Mod_lmn_ComputeAdvectionVelocity
   use typre
   use Mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersAdvectionVelocity
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersAdvectionVelocity(itask)
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
         
            ProcPointer%ComputeAdvectionVelocity => ComputeAdvectionVelocity
            
            !Advection Velocity
            if (a%kfl_advec == 0) then
               ProcPointer%ComputeAdvectionVelocity => NULLSUB
            endif
            
            !NonLinear Subscales
            if (a%kfl_nolsg == 1) then
               call ConcatenateProcedures(ProcPointer%ComputeAdvectionVelocity,NonLinearSGSAdvectionVelocity)
            endif
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   !Actual Computations
   !AdvectionVelocity
   subroutine ComputeAdvectionVelocity
      implicit none
      
      gpadv = gpvel(:,1)
   end subroutine
   
   !Non-linear subscales
   subroutine NonLinearSGSAdvectionVelocity
      implicit none
      
      gpadv = gpadv + a%vesgs(ielem)%a(:,1,e%igaus)
   end subroutine
   
end module
