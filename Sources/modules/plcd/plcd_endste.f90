subroutine plcd_endste(a,itask)
   !This routine ends a time step of the PLCD equations.
   use typre
   use def_parame
   use Mod_PLCD
   use Mod_plcd_stages 
   use Mod_int2str
   use Mod_SmoothedFieldGradient
   implicit none
   class(PLCDProblem), target :: a

   integer(ip) :: itask
   
   type(Stage), pointer :: cs => NULL()
   type(Substage), pointer :: css => NULL()
   
   integer(ip) :: imaterial, ndime
   real(rp) :: CurrentLoadFactor
   
   if (itask == 2) then
      !Update databases and postprocess data
      call a%EndsteElmope
      
      !if (a%UseSmoothedDisplacementGradient) call a%FilePostpr%postpr(a%SmoothedDisplacementGradient,'SmoothGrad'//int2str(a%itera),a%istep,a%ctime,a%Mesh)
      !call a%Mesh%GetNdime(ndime)
      !call PostprocessGaussPointGradient(a%Mesh,a%Memor,ndime,a%Displacement,'Gradient',a%FilePostpr,a%istep,a%ctime)
      
      !Update stages
      if (a%ctime > a%css%EndTime) then
         if (a%cs%CurrentSubstage < a%cs%NumberOfSubstages) then
            a%cs%CurrentSubstage =  a%cs%CurrentSubstage +1
            CurrentLoadFactor = a%css%CurrentLoadFactor
            a%css => a%cs%Substages(a%cs%CurrentSubstage)
            a%css%CurrentLoadFactor = CurrentLoadFactor
            if (a%TDData%kfl_WhenPerformTopologyOptimization == 2) a%TDData%kfl_PerformTopologyIteration = 1
         else
            a%CurrentStage = a%CurrentStage+1
            if (a%TDData%kfl_WhenPerformTopologyOptimization == 3) a%TDData%kfl_PerformTopologyIteration = 1
         endif
      endif
      
      if (a%CurrentStage > a%NumberOfStages+1) then
         call runend('current Stage not defined for plcd')
         a%kfl_stead = 1
      endif
      
      if (a%TDData%kfl_WhenPerformTopologyOptimization == 1) a%TDData%kfl_PerformTopologyIteration = 1
      
   endif

 
end subroutine plcd_endste
