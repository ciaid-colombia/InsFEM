module Mod_nsm_ComputeAdvectionVelocity
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersAdvectionVelocity 

   type, extends(PointerSetter) :: SPAdvectionVelocity
contains
      procedure :: SpecificSet => SpecificSetAdvectionVelocity
   end type
   type(SPAdvectionVelocity) :: SetPointersAdvectionVelocity
 
   procedure() :: NULLSUB
   logical :: isALE
   real(rp), allocatable :: elmve(:,:),gpmve(:)
   real(rp), pointer     :: meshve(:,:,:) => NULL()
   
contains

   subroutine SpecificSetAdvectionVelocity(d)
      implicit none
      class(SPAdvectionVelocity) :: d
         
      ProcPointer_ComputeAdvectionVelocity => ComputeAdvectionVelocity
      
      !Advection Velocity
      if (a%kfl_advec == 0) then
         ProcPointer_ComputeAdvectionVelocity => NULLSUB
      endif
      
      if (a%kfl_advco == 1) ProcPointer_ComputeAdvectionVelocity => ConstantAdvectionVelocity
      
      !NonLinear Subscales
      if (a%kfl_nolsg == 1) then
         call ConcatenateProcedures(ProcPointer_ComputeAdvectionVelocity,NonLinearSGSAdvectionVelocity)
      endif
      
      !ALE
      call a%Mesh%GetALE(isALE)
      if (isALE) then
         call a%Mesh%GetMeshVeloc(meshve)
         call ConcatenateProcedures(ProcHook_Initializations,AllocMeshVeloc)
         call ConcatenateProcedures(ProcHook_Gathers,GatherMeshVeloc)
         call ConcatenateProcedures(ProcHook_Interpolates,InterpolateMeshVeloc)
         call ConcatenateProcedures(ProcPointer_ComputeAdvectionVelocity,ALEAdvectionVelocity)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocMeshVeloc)
      end if
      
   end subroutine
   
   !-----------------------------------------------------------------------
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
  
   !Oseen 
   subroutine ConstantAdvectionVelocity
      implicit none
      
      gpadv(1:e%ndime) = a%advco(1:e%ndime)
   end subroutine
   
   !ALE
   subroutine AllocMeshVeloc
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elmve,'elmve','nsm_AdvectionVelocity')
      call a%Memor%alloc(e%ndime,gpmve,'gpmve','nsm_AdvectionVelocity')
   end subroutine
   
   subroutine GatherMeshVeloc
      implicit none
      
      call e%gather(e%ndime,elmve,meshve(:,:,1))
   end subroutine
   
   subroutine InterpolateMeshVeloc
      implicit none
      
      call e%interpg(e%ndime,elmve,gpmve)
   end subroutine
   
   subroutine ALEAdvectionVelocity
      implicit none
      
      gpadv = gpadv-gpmve
   end subroutine
   
   subroutine DeallocMeshVeloc
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,elmve,'elmve','nsm_AdvectionVelocity')
      call a%Memor%dealloc(e%ndime,gpmve,'gpmve','nsm_AdvectionVelocity')
   end subroutine
   
end module
