module Mod_tem_Advection
   use typre
   use Mod_tem_BaseElmope
   implicit none
  private
   public SetPointersAdvectionVelocity
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
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
         
            !--------------------------------------------------------------------
            !Advection
            if (a%kfl_advec /= 0) then
               !External velocity
               if (a%kfl_advec == 1) then 
                  call ConcatenateProcedures(ProcHook%Gathers,GatherExternalVeloc)
                  call ConcatenateProcedures(ProcHook%Interpolates,VelocInterpolates)
                  
                  if (a%kfl_ExternalAdvectionSGS == 1) then
                     call ConcatenateProcedures(ProcHook%Interpolates,VelocSGSInterpolates)
                  endif
               
               !Burgers equation
               elseif (a%kfl_advec == 9) then
                  call ConcatenateProcedures(ProcHook%Interpolates,BurgersVeloc)
               
               !Analytical velocity
               else   
                  call ConcatenateProcedures(ProcHook%Interpolates,GaussPointAnalyticalVeloc)
               endif
               call ConcatenateProcedures(ProcPointer%VelocityAndAgradV,VelocityNormAndAGradV)
            endif
            
            
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !--------------------------------------------------------------------
   !Actual computations
   subroutine GatherExternalVeloc
      !Gathers
      call e%gather(e%ndime,elvel(:,:),a%veloc(:,:))
   end subroutine
   
   subroutine VelocInterpolates
      !Interpolate
      call e%interpg(e%ndime,elvel(:,:),gpvel(:))
      gpadv(1:e%ndime) = gpvel(1:e%ndime)
   end subroutine
   
   subroutine VelocSGSInterpolates
      implicit none
      gpadv(1:e%ndime) = gpadv(1:e%ndime) + a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
   end subroutine
   
   subroutine VelocityNormAndAGradV
      implicit none
      
      !Advection velocity norm
      call vecnor(gpadv,e%ndime,gpvno,2)
      
      !Compute a·grad(V)
      call ComputeAGradV(e,gpadv,AGradV)
   end subroutine

   subroutine GaussPointAnalyticalVeloc
      implicit none
      
      call tem_elmvel(e,a%kfl_advec,gpvel)
      gpadv = gpvel
   end subroutine
   
   !Interpolation of the velocity gradient
   subroutine AllocGrvel
      implicit none
      
      call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','tem_elmope')    !Vel. gradient
   end subroutine   
   
   subroutine DeAllocGrvel
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','tem_elmope')    !Vel. gradient
   end subroutine  
   
   subroutine InterpolateGrvel
      implicit none
      
      call e%gradient(e%ndime,elvel(:,:),grvel)    !Vel. gradient
   end subroutine   
   
   subroutine BurgersVeloc
      implicit none
      
      gpvel = 0.0_rp
      gpvel(1) = gptem(1)
      gpadv = gpvel
   end subroutine
   

end module