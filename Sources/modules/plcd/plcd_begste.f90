subroutine plcd_begste(a)
   ! DESCRIPTION
   !    This routine prepares for a new time step of the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_PLCD
   use Mod_Mesh
   use Mod_TimeIntegrator
   use Mod_plcd_Stages
   implicit none
   class(PLCDProblem), target :: a
   real(rp), pointer :: TimeStep => NULL()
   integer(ip) :: ndime

   interface
      subroutine plcd_ComputeForcesVectors(a)
         import
         implicit none
         class(PLCDProblem) :: a
      end subroutine
   end interface

   if (a%CurrentStage > a%NumberOfStages) call runend('PLCD current time outside of all stages')
   
   a%cs => a%Stages(a%CurrentStage)
   a%css => a%cs%Substages(a%cs%CurrentSubstage)

   a%css%PreviousLoadFactor = a%css%CurrentLoadFactor
   a%css%CurrentLoadFactor = a%css%InitialLoadFactor + (a%ctime-a%css%Initime)/a%css%TimeInterval*(a%css%FinalLoadFactor-a%css%InitialLoadFactor)
   a%css%LoadFactorIncrement = a%css%CurrentLoadFactor-a%css%PreviousLoadFactor

   !Updating variables each time step
   
   if (a%kfl_TransientProblem == 1) then
      TimeStep => a%css%TimeStep
      a%Displacement(:,:,3) = a%Displacement(:,:,1) + TimeStep*a%Velocity(:,:,1) + TimeStep*TimeStep*0.5_rp*(1.0_rp - 2.0_rp*a%Beta)*a%Acceleration(:,:,1)
      a%Displacement(:,:,1) = a%Displacement(:,:,3)
      a%Velocity(:,:,3) = a%Velocity(:,:,1) + (1.0_rp - a%Gamma)*TimeStep*a%Acceleration(:,:,1)
      
      a%Acceleration(:,:,1) = (1.0_rp/(a%Beta*TimeStep*TimeStep))*(a%Displacement(:,:,1)-a%Displacement(:,:,3))
      a%Velocity(:,:,1) = a%Velocity(:,:,3) + a%Gamma*TimeStep*a%Acceleration(:,:,1)
      
   endif
   
   !Compute the required Forces Vectors
   a%itera = 0_ip  !So that we know that we are not at endite, but at begste
   call plcd_ComputeForcesVectors(a)

   

end subroutine plcd_begste

