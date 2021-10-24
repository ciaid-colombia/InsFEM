module Mod_nsc_pr_ArbitraryLagrangianEulerian
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersArbitraryLagrangianEulerian, ALEGradV

   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet

   !ALE
   real(rp), pointer :: meshve(:,:,:) => NULL()
   real(rp), allocatable :: elmve(:,:),gpmve(:)
   real(rp), allocatable :: ALEGradV(:)
   
contains

   subroutine SetPointersArbitraryLagrangianEulerian(itask)
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
         
               call a%Mesh%GetMeshVeloc(meshve)
               call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateMeshVelocity)
               call ConcatenateProcedures(ProcHook_nsc_pr_Gathers,GatherMeshVelocity)
               call ConcatenateProcedures(ProcHook_nsc_pr_Interpolates,InterpolateMeshVelocity)
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeVariables,ComputeALEGradients)
               call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateMeshVelocity)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateMeshVelocity
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elmve,'elmve','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,gpmve,'gpmve','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,ALEGradV,'ALEGradV','nsc_pr_elmope')
   end subroutine
   
   
   subroutine GatherMeshVelocity
      implicit none
      
      call e%gather(e%ndime,elmve,meshve(:,:,1))
   end subroutine
   
   subroutine InterpolateMeshVelocity
      implicit none
      
      call e%interpg(e%ndime,elmve,gpmve)
   end subroutine
   
   subroutine ComputeALEGradients
      implicit none
      
      !Compute mom·grad(V)
      call ComputeAGradV(e,gpmve,ALEGradV)

   end subroutine

   subroutine DeallocateMeshVelocity
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,elmve,'elmve','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,gpmve,'gpmve','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,ALEGradV,'ALEGradV','nsc_pr_elmope')
   end subroutine

end module
