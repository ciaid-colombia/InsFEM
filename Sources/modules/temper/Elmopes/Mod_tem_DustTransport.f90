 module Mod_tem_DustTransport
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_Advection
   implicit none
   private
   public SetPointersDustTransport,SetPointersDustTransportBoundary,ComputeSettlingVelocity
   
   integer(ip), allocatable :: kfl_IsSet, kfl_IsSetOnBoundary

   
contains

   !Computing the Settling Velocity, done once at the beginning of the simulation
   subroutine ComputeSettlingVelocity(a)
      class(TemperatureProblem) :: a
      
      real(rp) :: nu, wt,g,rhoair,Cd,d,wt0,Ret,raux
      integer(ip) :: j,ndime
      real(rp) :: acden,acsph,actco,acrea,acsou
      
      call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
      
      nu = a%DT_MediumDynamicViscosity/a%DT_MediumDensity
      wt = 1e-3
      d = a%DT_ParticleSize
      rhoair = a%DT_MediumDensity
      call a%Mesh%GetNdime(ndime)
      call vecnor(a%DT_GravityForce,ndime,g,2)
      
      do j = 1,100 !converges really fast, no need to check convergence
         wt0 = wt
         
         !Reynolds number of the particle flow
         Ret = wt*d/nu
         !Drag coefficient
         Cd = (24/Ret)*(1+0.15*(Ret**0.687))
         if (Cd < 0.2) Cd = 0.2
         wt = sqrt((4*(acden-rhoair)*g*d)/(3*rhoair*Cd))
         
         raux = abs((wt-wt0)/wt)
         if (raux < 1.0e-6 ) exit
      enddo

      
      a%DT_DustDepositionVelocityNorm = wt
      a%DT_DustDepositionVelocityVector = a%DT_GravityForce/g*wt
      
!       !Threshold wind velocity
!       if (a%DT_SurfaceWetness < 0.5) then
!          rhoair = a%DT_MediumDensity
!          a%ThresholdWindVelocity = 6.5*sqrt((acden-rhoair)*g*d/rhoair)*(1.2+0.2*log10(a%DT_SurfaceWetness))
!       else
!          a%ThresholdWindVelocity = 1e12_rp
!       endif
      
      !ThresholdFrictionVelocity
      a%DT_ThresholdFrictionVelocity = sqrt(0.0123*(   (acden/a%DT_MediumDensity)*g*d+3e-4_rp/(acden*d)   )) 
      
      
      
      
      
      
   end subroutine




   
!----------------------------------------------------------------------------
! ON ELEMENTS
!----------------------------------------------------------------------------

   !Setting Pointers
   subroutine SetPointersDustTransport(itask)
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
            
            !DustTransport
            if (a%DT_kfl_DustTransport >= 1)  then
               if (a%kfl_advec == 0) call runend('Dust Transport but no advection')
               
               !Advection is required first
               call SetPointersAdvectionVelocity(1)
               call ConcatenateProcedures(ProcHook%Interpolates,DustInterpolates)
               
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !
   subroutine DustInterpolates
      implicit none
      
      gpadv = gpadv + a%DT_DustDepositionVelocityVector(1:size(gpadv))
   end subroutine
   
!----------------------------------------------------------------------------
! ON BOUNDARIES
!----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersDustTransportBoundary(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSetOnBoundary)
         call a%Memor%allocObj(0,'kfl_IsSetOnBoundary','InitProcedurePointer',1)
         kfl_IsSetOnBoundary = -1
      
      case(1)
         if (kfl_IsSetOnBoundary == -1) then
            kfl_IsSetOnBoundary = 1
            
            !DustTransport
            if (a%DT_kfl_DustTransport >= 1)  then
               if (a%kfl_advec == 0) call runend('Dust Transport but no advection')
               
               !Advection is required first
               call SetPointersAdvectionVelocity(1)
               call ConcatenateProcedures(ProcHook%Interpolates,DustBoundaryFluxes)
               
            endif
         endif  
      case(100)
         deallocate(kfl_IsSetOnBoundary)
         call a%Memor%deallocObj(0,'kfl_IsSetOnBoundary','InitProcedurePointer',1)
      
      end select
   end subroutine  
   
   subroutine DustBoundaryFluxes
      real(rp) :: raux, delta,tvel(e%ndime),tveno,vikin,velfr,Fi
      
      !Dust Transport
      if (a%kfl_fixbo(iboun) == 5) then
         
         delta = 1000 !Distance to the wall, we consider the asymptotic behavior of the wall law
         
         raux = dot_product(e%baloc(1:e%ndime,e%ndime),gpvel(1:e%ndime))
         
         !Tangent velocity norm 
         tvel = gpvel - e%baloc(1:e%ndime,e%ndime)*raux
         call vecnor(tvel,e%ndime,tveno,2)
         vikin = a%DT_MediumDynamicViscosity/a%DT_MediumDensity
         
         !Friction Velocity
         call frivel(delta,tveno,vikin,velfr)                  ! U*
      
         !Erosion dust flux
         if (velfr > a%DT_ThresholdFrictionVelocity) then
            Fi = 1.4e-6_rp*((velfr**4)*(1-a%DT_ThresholdFrictionVelocity/velfr))
         else
            Fi = 0.0_rp
         endif
         
         tract = tract + acden*Fi
      endif      
         
      
      
   
   end subroutine
   

   
end module

 
