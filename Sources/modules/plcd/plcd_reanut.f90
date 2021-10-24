subroutine plcd_reanut(a,itask)
   use typre
   use Mod_Listen
   use Mod_PLCD
   implicit none

   class (PLCDProblem) :: a
   integer(ip) :: itask

   integer(ip) :: istab,ipart,rpalo

   !For Listener
   real(rp), pointer     :: param(:)
   character(5), pointer :: words(:)
   integer(ip), pointer  :: nnpar,nnwor

   integer(ip) :: ipara


   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !  Initializations (defaults)
      a%UseSmoothedDisplacementGradient = .false.
      a%staco(1)  = 4.0_rp
      a%ErrorEstimatorTypeOfSubscales = 0  !0: Orthogonal, 1: ASGS. Default is orthogonal subscales


   elseif (itask == 1) then     !Inside the Numerical Treatment Block
      if (words(1) == 'SMOOT' ) then
         if (a%Listener%exists('ON   ')) then
            a%UseSmoothedDisplacementGradient = .true.
         else
            a%UseSmoothedDisplacementGradient = .false.
         endif
      elseif (words(1) == 'UPFOR' ) then
         if (a%Listener%exists('ON   ')) then
            a%UseUPFormulation = .true.
         else
            a%UseUPFormulation = .false.
         endif
      elseif (words(1) == 'LINES') then
         a%kfl_LineSearch = 1
         a%LineSearchRadiusOfConvergence = param(1)
         a%LineSearchMaxIterations = param(2)
      !For Error estimator
      elseif(words(1)=='ERROR')then
         a%NLayersRefinement = a%Listener%getint('NLAYE',0_ip,'OutLayers')    !Preliminary frequency
         a%GeneralRefinementLevels = a%Listener%getint('GENER',0_ip,'GeneralRefinementLevels')    !Preliminary frequency
         a%InterfaceRefinementLevels = a%Listener%getint('INTER',0_ip,'InterfaceRefinementLevels')    !Preliminary frequency

      elseif (words(1) == 'ERRSU') then
         if (a%Listener%exists('ORTHO')) then
            a%ErrorEstimatorTypeOfSubscales = 0
         elseif (a%Listener%exists('ASGS ')) then
            a%ErrorEstimatorTypeOfSubscales = 1
         endif
      elseif(words(1)=='STACO') then
         do ipara = 1,size(a%staco)
            a%staco(ipara) = param(ipara)
         end do
      endif

   elseif (itask == 100) then   !Finalize reading operations

   endif

end subroutine plcd_reanut

