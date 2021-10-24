subroutine lev_reanut(a,itask)
   use typre
   use Mod_Listen
   use Mod_LevelSet
   implicit none
   
   class (LevelSetProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   integer(ip) :: istab
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   
   !Initializations
   if(itask == 0)then
   
      !to work whith the CDR options 
      a%staco(1)  = 4.0_rp                           ! Diffusive term
      a%staco(2)  = 2.0_rp                           ! Convective term
      a%staco(3)  = 1.0_rp                           ! Reactive term
      
      a%kfl_ReinitLevel=0
      a%NstepsReinitLevel = 50
      a%kfl_CutCriteria=0
      a%kfl_MassCorrection = 0
   
   
   elseif(itask == 1)then
   
      if(words(1)=='STABI')then
         do istab = 1,3
            a%staco(istab) = param(istab)
         end do
      
      !Initial redistancing
      elseif(words(1) == 'INITR') then
         if(a%Listener%exists('ON   ')) then 
            a%kfl_InitialRedistance = 1
         else
            a%kfl_InitialRedistance = 0
         endif
      
      elseif(words(1)=='REINI')then
         if(words(2)=='ZIGZA')then
            a%kfl_ReinitLevel=1
         elseif(words(2)=='DISTA')then
            a%kfl_ReinitLevel=2
         elseif(words(2)=='POISS') then
            a%kfl_ReinitLevel=3
         end if

      elseif(words(1)=='NSTEP')then
         a%NstepsReinitLevel = param(1)

      elseif(words(1)=='CUTCR')then
         if(words(2)=='BLEND') a%kfl_CutCriteria=1

      !For Error estimator      
      elseif(words(1)=='ERROR')then
         a%OutLayers = a%Listener%getrea('PARA1',0.0_rp,'OutLayers')    !Preliminary frequency
         a%GeneralRefinementLevels = a%Listener%getrea('PARA2',0.0_rp,'GeneralRefinementLevels')    !Preliminary frequency
         a%InterfaceRefinementLevels = a%Listener%getrea('PARA3',0.0_rp,'InterfaceRefinementLevels')    !Preliminary frequency
      
      !For FixedMeshALE options
      elseif (words(1) == 'NOALE') then
         if(a%Listener%exists('FORCE')) a%kfl_ForceEulerianAdvection = 1
         
      !For Mass Correction
      elseif (words(1) == 'MASSC') then
         if(a%Listener%exists('ON   ')) a%kfl_MassCorrection = 1   
         
      
      end if
      
   elseif(itask == 100)then  

   
   endif

   
end subroutine   
