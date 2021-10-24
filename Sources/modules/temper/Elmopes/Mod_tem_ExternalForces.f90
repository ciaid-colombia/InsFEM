module Mod_tem_ExternalForces
   use typre
   use Mod_tem_BaseElmope
   implicit none
  private
   public SetPointersExternalForces
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersExternalForces(itask)
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
            
            ProcPointer%ExternalForces => ExternalForces
         
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !--------------------------------------------------------------------
   !Actual computations
   !---------------------------------------------------------------
   !For external Forces
   subroutine ExternalForces
      implicit none
      
      call tem_ComputeExternalForces(e,acsou,acsph,elext)
   end subroutine
   
   subroutine InterpolateCoord
      implicit none
      
      call e%interpg(e%ndime,e%elcod,gpcod)
   end subroutine
   
end module
