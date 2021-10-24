subroutine lmn_output(a,itask)
   use typre
   use Mod_Postpr
   use Mod_Stream
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   integer(ip)           :: ndime,itask
   integer(ip), save     :: dopost(12)

   interface
      subroutine nsi_outtpo(a)
         use typre
         import LowMachProblem
         implicit none         
         class(LowMachProblem) :: a
      end subroutine
   end interface
   
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)

   if (dopost(1) == 1) then
      call a%FilePostpr%postpr(a%veloc(:,:,1),'Velocity',a%istep,a%ctime,a%Mesh)
   end if

   if (dopost(2) == 1) then
      call a%FilePostpr%postpr(a%press(:,1),'Pressure',a%istep,a%ctime,a%Mesh)
   end if

   if (dopost(3) == 1) then
      call a%FilePostpr%postpr(a%tempe(:,1),'Temperature',a%istep,a%ctime,a%Mesh)
   end if

   if (dopost(4) == 1) then
      call a%ComputeDensity
      call a%FilePostpr%postpr(a%densf(:),'Density',a%istep,a%ctime,a%Mesh)
   end if

   if (dopost(5) == 1) then
      call a%FilePostpr%postpr(a%pther(1),'Therm_Pressure',a%istep,a%ctime,a%Mesh)
   end if

   if (dopost(7) == 1) then
      call a%Mesh%GetNdime(ndime)
      call stream(a%veloc,ndime,a%istep,a%ctime,a%Mesh,a%Memor,a%FilePostpr,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPISize)
   end if

   !Subscales
   if (dopost(8) == 1) then
      call a%FilePostpr%postgp(a%vesgs,'VelocitySGS',a%istep,a%ctime,a%Mesh)
   end if
   if (dopost(9) == 1) then
      call a%FilePostpr%postgp(a%tesgs,'TemperatureSGS',a%istep,a%ctime,a%Mesh,'SCALAR')
   end if
   if (dopost(10) == 1) then
      call a%FilePostpr%postgp(a%prsgs,'PressureSGS',a%istep,a%ctime,a%Mesh)
   end if

   !Residual
   if (dopost(11) == 1) then
      call a%FilePostpr%postgp(a%residualU,'ResidualU',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postgp(a%residualP,'ResidualP',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postgp(a%residualT,'ResidualT',a%istep,a%ctime,a%Mesh)
   end if    

   !Repro
   if (dopost(12) == 1) then
      call a%Mesh%GetNdime(ndime)
      call a%FilePostpr%postpr(a%repro(1:ndime,:),'ResidualProjU',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%repro(ndime+2,:),'ResidualProjP',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postpr(a%repro(ndime+1,:),'ResidualProjT',a%istep,a%ctime,a%Mesh)
   end if   

   select case(itask)
   case(0)
      a%pos_alrea=0
      !Tracking of points.
      if(a%nptra > 0) then
         call lmn_outtpo(a)
      end if
   end select

end subroutine lmn_output
