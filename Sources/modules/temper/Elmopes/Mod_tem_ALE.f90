module Mod_tem_ALE
   use typre
   use mod_tem_BaseElmope
   use Mod_tem_Advection
   implicit none
   private
   public SetPointersALE

   integer(ip), allocatable :: kfl_IsSet
   
   !ALE
   logical :: isALE
   real(rp), allocatable :: elmve(:,:),gpmve(:)
   real(rp), pointer :: meshve(:,:,:) => NULL()
   
contains
   
   !Set Pointers
   subroutine SetPointersALE(itask)
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
            !ALE
            call a%Mesh%GetALE(isALE)
            if (isALE) then
               !We need the advection velocity
               call SetPointersAdvectionVelocity(1)
               
               call a%Mesh%GetMeshVeloc(meshve)
               call ConcatenateProcedures(ProcHook%Initializations,AllocMeshVeloc)
               call ConcatenateProcedures(ProcHook%Gathers,GatherMeshVeloc)
               call ConcatenateProcedures(ProcHook%Interpolates,InterpolateMeshVeloc)
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocMeshVeloc)
            end if
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine
   
    
   !------------------------------------------------------
   !ALE
   subroutine AllocMeshVeloc
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elmve,'elmve','tem_elmope')
      call a%Memor%alloc(e%ndime,gpmve,'gpmve','tem_elmope')
   end subroutine
   
   subroutine GatherMeshVeloc
      implicit none
      
      call e%gather(e%ndime,elmve,meshve(:,:,1))
   end subroutine
   
   subroutine InterpolateMeshVeloc
      implicit none
      
      call e%interpg(e%ndime,elmve,gpmve)
   end subroutine
   
   subroutine DeallocMeshVeloc
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,elmve,'elmve','tem_elmope')
      call a%Memor%dealloc(e%ndime,gpmve,'gpmve','tem_elmope')
   end subroutine   
end module
