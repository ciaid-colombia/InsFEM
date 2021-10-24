module Mod_tem_ShockCapturing
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_TempeGradient
   use Mod_tem_Advection
   implicit none
   private
   public SetPointersShockCapturing
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
   integer(ip) :: ndime
   real(rp), allocatable :: elgrtem(:,:), elvis(:)
   
   real(rp) :: Projectedgrtem(3), OrthogonalGrtem(3),NormGrTem,NormOGrTem

   
contains

   subroutine SetPointersShockCapturing(itask)
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
            
            if (a%kfl_shock == 0) then
            
            elseif (a%kfl_shock == 3) then
               call ConcatenateProcedures(ProcHook%Initializations,SCInitializations)
               call SetPointersComputeTempeGradient(1)
               call SetPointersAdvectionVelocity(1)
               call ConcatenateProcedures(ProcHook%Gathers,SCGathers)
               call ConcatenateProcedures(ProcHook%Interpolates,SCInterpolates)
               
               call ConcatenateProcedures(ProcHook%PhysicalProp,SCPhysicalProperties)
               
               call ConcatenateProcedures(ProcHook%Finalizations,SCFinalizations)
            else
               call runend('This type of shock capturing is not implemented') 
            endif
           
            
            
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)  
      
      end select
   end subroutine
   
   
   !--------------------------------------------------------------------
   !Actual computations
   !---------------------------------------------------------------
   subroutine SCInitializations
      implicit none
      
      call a%Mesh%GetNdime(ndime)
      
      call a%Memor%alloc(ndime,e%mnode,elgrtem,'elgrtem','ShockCapturing')
   end subroutine
   
   subroutine SCGathers
      implicit none
      
      call e%Gather(ndime,elgrtem,a%gradient)
   end subroutine
   
   subroutine SCInterpolates
      implicit none
      
      call e%interpg(ndime,elgrtem,Projectedgrtem)
      OrthogonalGrTem(1:ndime) = grtem(1:ndime) - Projectedgrtem(1:ndime)
   end subroutine   
   
   subroutine SCPhysicalProperties
      implicit none
      
      real(rp) :: physpar,alpha
      
      call elmchl(e,a%kfl_advec,elvel,chale)
      
      call vecnor(ProjectedGrTem,ndime,NormGrTem,2_ip)
      call vecnor(OrthogonalGrTem,ndime,NormOGrTem,2_ip)
      
      alpha = max(0.0_rp,a%shock - (2*acvis)/(gpvno*chale(1)))
      if (NormgrTem >  1e-12_rp) then
         physpar = alpha*0.5_rp*(acvis + gpvno*chale(1) + acrcp*e%hleng(2)**2)*NormOGrtem/NormGrTem
      else
         physpar = 0.0_rp
      endif
      
      acvis = acvis + physpar
      
      !For Postprocessing
      if (a%npp_stepi(8) /= 0) then
         a%ShockCapturingViscosity(ielem)%a(e%igaus) = physpar
         a%GradientGaussPoints(ielem)%a(:,e%igaus) = grtem(1:ndime)
      endif
      
   end subroutine
   
   
   subroutine SCFinalizations
      implicit none
      
      call a%Memor%dealloc(ndime,e%mnode,elgrtem,'elgrtem','ShockCapturing')
   end subroutine
  
end module