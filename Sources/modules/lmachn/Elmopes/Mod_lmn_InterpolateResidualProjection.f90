module Mod_lmn_InterpolateResidualProjection
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_ComputeResidualProjection
   implicit none
   private
   public SetPointersInterpolateResidualProjection, gprep
   
   real(rp), allocatable :: elrep(:,:),gprep(:) !The residual projection

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersInterpolateResidualProjection(itask)
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
         
         if (a%kfl_repro /= 0) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocRep)
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocRep)
               call ConcatenateProcedures(ProcHook%Gathers,GatherRep)
               call ConcatenateProcedures(ProcHook%Interpolates,InterpolateRep)
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine AllocRep
      implicit none
      
      call a%Memor%alloc(a%ndofn,e%mnode,elrep,'elrep','lmn_EndElmope')
      call a%Memor%alloc(a%ndofn,gprep,'gprep','lmn_EndElmope')
   end subroutine
   
   subroutine DeallocRep
      implicit none
      
      call a%Memor%Dealloc(a%ndofn,e%mnode,elrep,'elrep','lmn_EndElmope')
      call a%Memor%Dealloc(a%ndofn,gprep,'gprep','lmn_EndElmope')
   end subroutine
   
   !Residual Interpolation Subroutines
   subroutine GatherRep
      implicit none
      
      call e%gather(a%ndofn,elrep,a%repro)
   end subroutine

   subroutine InterpolateRep
      implicit none
      
      !Interpolate
      call e%interpg(a%ndofn,elrep,gprep)
   end subroutine
   
end module
