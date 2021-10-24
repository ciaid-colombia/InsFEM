!New Driver Keys should be added here
subroutine CreateDriverFromKey(key,DriverCreator)
   use Mod_DriverCreator

#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_NSTINC) || !defined MODULE_SELECTION_ON
   use Mod_NstincDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_PODROM) || !defined MODULE_SELECTION_ON   
   use Mod_PodRomDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_TEMPE) || !defined MODULE_SELECTION_ON
   use Mod_TempeDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_OPTICS) || !defined MODULE_SELECTION_ON   
   use Mod_OpticsDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_ALE) || !defined MODULE_SELECTION_ON   
   use Mod_ALEDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_LOWMACH) || !defined MODULE_SELECTION_ON   
   use Mod_LowMachDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_NSCOMP) || !defined MODULE_SELECTION_ON   
   use Mod_NSCompDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_SOLIDS) || !defined MODULE_SELECTION_ON   
   use Mod_SolidsDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_LEVEL) || !defined MODULE_SELECTION_ON   
   use Mod_LevelDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_PLCD) || !defined MODULE_SELECTION_ON   
   use Mod_PLCDDriver
#endif   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_RVEPROP) || !defined MODULE_SELECTION_ON   
   use Mod_RVEpropDriver
#endif   
   
   
   implicit none
   character(5) :: key
   type(DriverCreatorType) :: DriverCreator
   
   !Name of the Driver
   DriverCreator%key = key   
   DriverCreator%Driver => NULL()
   
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_NSTINC) || !defined MODULE_SELECTION_ON   
   if (key == 'NSTIN') then
      DriverCreator%Driver => NstincDriver_Const()
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_TEMPE) || !defined MODULE_SELECTION_ON   
   if (key == 'TEMPE') then
      DriverCreator%Driver => TempeDriver_Const()   
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_OPTICS) || !defined MODULE_SELECTION_ON   
   if (key == 'OPTIC') then
      DriverCreator%Driver => OpticsDriver_Const()   
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_LOWMACH) || !defined MODULE_SELECTION_ON   
   if (key == 'LMACH') then
      DriverCreator%Driver => LowMachDriver_Const()  
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_NSCOMP) || !defined MODULE_SELECTION_ON   
   if (key == 'NSCOM') then
      DriverCreator%Driver => NSCompDriver_Const()   
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_ALE) || !defined MODULE_SELECTION_ON   
   if (key == 'ALEPR') then
      DriverCreator%Driver => ALEDriver_Const()    
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_LEVEL) || !defined MODULE_SELECTION_ON   
   if (key == 'LEVSE') then
      DriverCreator%Driver => LevelDriver_Const()   
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_SOLIDS) || !defined MODULE_SELECTION_ON   
   if (key == 'SOLID') then
      DriverCreator%Driver => SolidsDriver_Const()   
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_PODROM) || !defined MODULE_SELECTION_ON   
   if (key == 'PODRO') then
      DriverCreator%Driver => PodRomDriver_Const()
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_PLCD) || !defined MODULE_SELECTION_ON   
   if (key == 'PLCDP') then
      DriverCreator%Driver => PLCDDriver_Const()
   endif
#endif    
     
#if (defined MODULE_SELECTION_ON && defined SELECT_MODULE_RVEPROP) || !defined MODULE_SELECTION_ON   
   if (key == 'RVEPR') then
      DriverCreator%Driver => RVEpropDriver_Const()
   endif
#endif    
   
   
   if (.not. associated(DriverCreator%Driver)) then
      DriverCreator%key = 'NODRI'
   else
      !Also compose the words which signal the end of the reading for the driver
      DriverCreator%endword(1:3) = 'END'
      DriverCreator%endword(4:5) = key(1:2)
   endif
         
end subroutine      



