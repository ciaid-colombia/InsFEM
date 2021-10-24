module Mod_tem_TempeGradient
   use typre
   use Mod_tem_BaseElmope
   implicit none
   private
   public SetPointersComputeTempeGradient, grtem

   
   integer(ip), allocatable :: kfl_IsSet
   
   real(rp), allocatable :: grtem(:)  
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeTempeGradient(itask)
      implicit none
      integer(ip) :: itask

      integer(ip) :: kfl_nonlinear
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1) 
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            call ConcatenateProcedures(ProcHook%Initializations,Allocations)
            call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGrtem)
            call ConcatenateProcedures(ProcHook%Finalizations,Deallocations)
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
   
   
   !--------------------------------------------------------------------
   !Temperature Gradient interpolation
   subroutine Allocations
      implicit none
      
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,grtem,'grtem','tem_TempeGradient')
   end subroutine
   
   subroutine InterpolateGrtem
      implicit none

      call e%gradient(1,eltem,grtem)    !Temperature gradient
   end subroutine
   
   subroutine DeAllocations
      implicit none
      
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      call a%Memor%dealloc(ndime,grtem,'grtem','tem_TempeGradient')
   end subroutine
   
   
end module 
