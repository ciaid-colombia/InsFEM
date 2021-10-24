module Mod_tem_InterpolateResidualProjection
   use typre
   use mod_tem_BaseElmope
   implicit none
   private
   public SetPointersInterpolateResidualProjection,gprep

   integer(ip), allocatable :: kfl_IsSet
   
   !Residual Projections   
   real(rp), allocatable :: elrep(:)
   real(rp)              :: gprep(1)
   
contains
   
   !Set Pointers
   subroutine SetPointersInterpolateResidualProjection(itask)
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
         
            if (a%kfl_repro == 1) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocRep)
               call ConcatenateProcedures(ProcHook%Gathers,GatherRep)
               call ConcatenateProcedures(ProcHook%Interpolates,InterpolateRep)
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocRep)
            endif 
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
   !-------------------------------------------------------------------
   !FOR RESIDUAL PROJECTION
   subroutine AllocRep
      implicit none

      call a%Memor%alloc(e%mnode,elrep,'elrep','tem_elmope')
      gprep = 0.0_rp
   end subroutine

   subroutine DeallocRep
      implicit none

      call a%Memor%dealloc(e%mnode,elrep,'elrep','tem_elmope')
   end subroutine

   subroutine GatherRep
      implicit none
      
      call e%gather(1,elrep,a%repro)
   end subroutine

   subroutine InterpolateRep
      implicit none
      
      !Interpolate
      call e%interpg(1,elrep,gprep)
   end subroutine

 

end module