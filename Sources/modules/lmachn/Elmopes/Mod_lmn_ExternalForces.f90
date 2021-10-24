module Mod_lmn_ExternalForces
   use typre
   use Mod_lmn_BaseElmope
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
            
            ProcPointer%ExternalForces  => lmnForces
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !Force term   
   !Non Exact Solution Case
   subroutine lmnForces
      implicit none    
      !Compute vector of external forces
      call lmn_ComputeExternalForces_mom(e,acden,a%grnor,a%gravi,elext_mom)  
      elext_pre = elext_mom 
      call lmn_ComputeExternalForces_ene(e,a%kfl_sourc,acsph,acsou,gpsou(1),elext_ene)  
   end subroutine  
end module

