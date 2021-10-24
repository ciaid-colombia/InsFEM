module Mod_supm_ComputeAdvectionVelocity
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersAdvectionVelocity
   !procedure() :: NULL()
   integer(ip), allocatable :: kfl_IsSet
   logical :: isALE
   real(rp), allocatable :: elmve(:,:),gpmve(:)
   real(rp), pointer     :: meshve(:,:,:) => NULL()

contains

   subroutine SetPointersAdvectionVelocity(itask)
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
             ProcPointer%ComputeAdvectionVelocity_sup => ComputeAdvectionVelocity
             ProcPointer%ComputeAdvectionVelocityNorm_sup => ComputeAdvectionVelocityNorm
             ProcPointer%ComputeAGradV_sup =>ComputeAGradVAdvec
         endif

            !ALE
            call a%Mesh%GetALE(isALE)
            if (isALE) then
               call a%Mesh%GetMeshVeloc(meshve)
               call ConcatenateProcedures(ProcHook_Initializations,AllocMeshVeloc)
               call ConcatenateProcedures(ProcHook_Gathers,GatherMeshVeloc)
               call ConcatenateProcedures(ProcHook_Interpolates,InterpolateMeshVeloc)
               call ConcatenateProcedures(ProcPointer%ComputeAdvectionVelocity_sup,ALEAdvectionVelocity)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocMeshVeloc)
            end if
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   !Actual Computations
   subroutine ComputeAdvectionVelocity
      implicit none
      
      gpadvec = 0.0_rp
      gpadvec_VE = gpvel(:,1) !advection used in constitutive equation
      

      if (a%kfl_advec == 0) then
         a%staco(2)=0.0_rp
      else
         gpadvec=gpadvec_VE !advection used in momentum equation
      endif
     
   end subroutine
   
   subroutine ComputeAdvectionVelocityNorm
      implicit none
      
      call vecnor(gpadvec_VE,e%ndime,gpvno_VE,2)
      if (a%kfl_advec==1) then
         gpvno=gpvno_VE
      endif
      
   end subroutine   
   
   subroutine ComputeAGradVAdvec
      implicit none
      
      
      call ComputeAGradV(e,gpadvec_VE,AGradV_VE)
      
      if (a%kfl_advec==1) then
         AGradV=AGradV_VE
      endif
      
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
      
      gpadvec = gpadvec-gpmve
   end subroutine
   
   subroutine DeallocMeshVeloc
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,elmve,'elmve','nsm_AdvectionVelocity')
      call a%Memor%dealloc(e%ndime,gpmve,'gpmve','nsm_AdvectionVelocity')
   end subroutine
   
end module
