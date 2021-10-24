 module Mod_sld_ExternalForces
   use typre
   use Mod_sld_BaseElmope
   implicit none
   private
   public SetPointersExternalForces

   integer(ip), allocatable :: kfl_IsSet 
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
 
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
            
            ProcPointer%ExternalForces  => sldForces      
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine  sldForces
      implicit none    
      !Compute vector of external forces
      call sld_ComputeExternalForces(e,densi,a%grnor,a%gravi,a%traction,elext,dvol,a%force_factor)
   end subroutine   
        
end module
