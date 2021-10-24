module Mod_tem_ComputeResidualProjection
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_ComputeGPResidual
   implicit none
   private
   public SetPointersComputeResidualProjection
   
   integer(ip), allocatable :: kfl_IsSet
   
   real(rp), allocatable :: elres(:)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeResidualProjection(itask)
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
            
            !Residual Projection
            call ConcatenateProcedures(ProcHook%PreLoop,PreLoopRep)
            call ConcatenateProcedures(ProcHook%Initializations,InitializationsRep)
            call ConcatenateProcedures(ProcHook%PreGauss,ElmatsToZeroRep)
            call SetPointersComputeGpResidual(1)
            if (kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeRHSNonLinearRes)
            call ConcatenateProcedures(ProcHook%InGaussElmats,GpresAssemblyToelres)
            
            call ConcatenateProcedures(ProcHook%Assembly,AssemblyResidual)
            call ConcatenateProcedures(ProcHook%Finalizations,FinalizationsRep)
            call ConcatenateProcedures(ProcHook%PostLoop,ProjectResidual)
            
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
   !-------------------------------------------------------------------
   !FOR RESIDUAL PROJECTION
   subroutine PreLoopRep
      implicit none
      
      a%repro = 0.0_rp
   end subroutine
   
   subroutine InitializationsRep
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%mnode,elres,'elres','nsm_EnditeElmope')
   end subroutine

   subroutine FinalizationsRep
      implicit none
      
      !Matrices dealloc
      call a%Memor%dealloc(e%mnode,elres,'elres','nsm_EnditeElmope')
   end subroutine

   subroutine ElmatsToZeroRep
      implicit none
      
      elres = 0.0_rp
   end subroutine
   
  
   
   subroutine GpresAssemblyToelres
      implicit none
      call tem_elmrep(e,dvol,gpres,elres)
      
   end subroutine
   
   subroutine ComputeRHSNonLinearRes
      implicit none
   
      call tem_elmrfe_oto_nonlinear(e,acvis,eltem,gpres)
      
   end subroutine
   
   subroutine AssemblyResidual
      implicit none
      
      !a%repro(e%lnods(1:e%pnode)) = a%repro(e%lnods(1:e%pnode)) + elres(1:e%pnode)
      call a%Mesh%AssemblyToArray(e,1_ip,elres,a%repro)
   end subroutine

   subroutine ProjectResidual
      implicit none
      call a%Project(1,a%repro) 
      !call a%FilePostpr%postpr(a%repro,'repro',a%istep,a%ctime,a%Mesh)
   end subroutine
   
   
end module 
