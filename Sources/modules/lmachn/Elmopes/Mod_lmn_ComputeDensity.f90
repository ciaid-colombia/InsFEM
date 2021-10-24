module Mod_lmn_ComputeDensity
   use typre
   use Mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersDensity
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersDensity(itask)
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
         
            ProcPointer%ComputeDensity => ComputeDensity
            
            !NonLinear Subscales
            if (a%kfl_nolsg == 1) then
               ProcPointer%ComputeDensity => NonLinearSGSDensity
            endif
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   !Actual Computations
   !Density
   subroutine ComputeDensity
      implicit none
      integer(ip)    :: itime
      gpden(1) = a%pther(1)/a%sgasc/gptem(1)
      gpden(2) = a%pther(3)/a%sgasc/gptem(2)
      do itime = 3,nsteps !Time bdf2 and others
         gpden(itime) = a%pther(itime+1)/a%sgasc/gptem(itime)
      enddo 
   end subroutine
   
   !Non-linear subscales
   subroutine NonLinearSGSDensity
      implicit none 
      integer(ip)    :: itime
      gpden(1) = a%pther(1)/a%sgasc/(gptem(1) + a%tesgs(ielem)%a(1,e%igaus))
      gpden(2) = a%pther(3)/a%sgasc/(gptem(2) + a%tesgs(ielem)%a(2,e%igaus))
      do itime = 3,nsteps !Time bdf2 and others
         gpden(itime) = a%pther(itime+1)/a%sgasc/(gptem(itime) + a%tesgs(ielem)%a(itime,e%igaus))
      enddo 
   end subroutine
   
end module
