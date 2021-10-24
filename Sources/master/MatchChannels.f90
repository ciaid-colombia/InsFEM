module Mod_MatchChannels
   use def_master
   use Mod_DCHashCharSize
   use Mod_DistributedContainer
   use Mod_GeneralCase
   use Mod_DriverInterface
   use Mod_InChannel
   use Mod_DC_GeneralCase
   use Mod_DC_Driver
   use Mod_DC_InChannel
   use Mod_CaseInterpolator
   use Mod_DC_CaseInterpolator
   use Mod_MeshInterpolator
   implicit none
   private
   public LoopInCases
   
   character(DCHashCharSize), save :: CurrentInCaseKey, CurrentInDriverKey, CurrentInChannelKey
   type(DistributedContainerList), pointer :: CurrentInCase_CaseInterpolatorList => NULL()
   
contains

   subroutine LoopInCases(myDC)
      class(DistributedContainer), pointer :: myDC
      type(GeneralCase), pointer :: myCase => NULL()
      
      
      !Store the current InCase key just in case we need it
      call myDC%Getkey(CurrentInCaseKey)

      !Extract case
      call ExtractGeneralCase(myDC,myCase)
      
      !Loop through all the case interpolators
      call myCase%CaseInterpolatorList%LoopList(LoopInCaseInterpolators)
      CurrentInCase_CaseInterpolatorList => myCase%CaseInterpolatorList
      
      !Loop through all my drivers
      call myCase%DriverList%LoopList(LoopInDrivers)
   end subroutine
   
   subroutine LoopInCaseInterpolators(myDC)
      class(DistributedContainer), pointer :: myDC
      type(CaseInterpolator), pointer :: myCaseInterpolator => NULL()
      
      class(DistributedContainer), pointer :: NeighborCaseDC => NULL()
      character(DCHashCharSize) :: NeighborCaseKey
      type(GeneralCase), pointer :: NeighborCase => NULL()
      
      call ExtractCaseInterpolator(myDC,myCaseInterpolator)
         
      !See who is the other case
      call myDC%GetKey(NeighborCaseKey)
      call CaseList%GetFromKey(NeighborCaseKey,NeighborCaseDC)
      if (.not. associated(NeighborCaseDC)) call runend('NeighborCase not found when looking for interpolators')
      call ExtractGeneralCase(NeighborCaseDC,NeighborCase)
      
      call myCaseInterpolator%SetNeighborCase(NeighborCase)
   end subroutine   
   
   subroutine LoopInDrivers(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
   
      !Extract key just in case we need it
      call myDC%Getkey(CurrentInDriverKey)
      !Extract driver
      call ExtractDriver(myDC,myDriver)
      
      !Match the Channels
      call myDriver%InChannelList%LoopList(LoopInChannels)
         
   end subroutine
   
   subroutine LoopInChannels(myDC)
      class(DistributedContainer), pointer :: myDC
      type(InChannel), pointer :: myInChannel => NULL()
      
      class(DistributedContainer), pointer :: myDC2 => NULL()
      character(DCHashCharSize) :: OutCaseKey,OutDriverKey,OutChannelKey
      
      type(GeneralCase), pointer          :: myOutCase => NULL()
      class(DriverInterface), pointer      :: myOutDriver => NULL()
      class(DistributedContainer), pointer :: myOutChannel => NULL()
      type(CaseInterpolator), pointer :: myCaseInterpolator => NULL()
      type(Interpolator), pointer :: InterpolatorPointer => NULL()
      integer(ip), pointer        :: npoinPointer => NULL()
      
      !call myDC%GetNick(CurrentInChannelKey)
      call myDC%GetKey(CurrentInChannelKey)
      call ExtractInChannel(myDC,myInChannel)
      
      !Here I am, looping for each InChannel
      call myInChannel%GetCaseKey(OutCaseKey)
      !If the OutCaseKey was not specified, then it is an intracase channel
      if (OutCaseKey == "") then
         OutCaseKey = CurrentInCaseKey
      endif
      
      !Find the specified case
      call CaseList%GetFromKey(OutCaseKey,myDC2)
      if (.not. associated(myDC2)) call SendRunendMessage(CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey)
      call ExtractGeneralCase(myDC2,myOutCase)
      
      !Find the specified driver
      call myInChannel%GetDriverNick(OutDriverKey)
      !call myInChannel%GetDriverKey(OutDriverKey)
      call myOutCase%DriverList%GetFromKey(OutDriverKey,myDC2)
      if (.not. associated(myDC2)) call SendRunendMessage(CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey,OutDriverkey)
      call ExtractDriver(myDC2,myOutDriver)
      
      !Find the specified channel
      call myInChannel%GetChannelKey(OutChannelKey)
      call myOutDriver%OutChannelList%GetFromKey(OutChannelKey,myDC2)
      if (.not. associated(myDC2) .and. myInChannel%Priority == 'OPTIONAL') then 
          call SendOptionalNPMessage(CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey,OutDriverKey,OutChannelKey)
      elseif (.not. associated(myDC2) .and. myInChannel%Priority == 'MAIN') then 
          call SendRunendMessage(CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey,OutDriverKey,OutChannelKey)
      endif
      
      !Point the channel
      myInChannel%DC => myDC2

      !We also need to Initialize the interpolator info for the channel if the channel is interpolated
      !Only for intercase channels
      if (OutCaseKey /= CurrentInCasekey) then
         call CurrentInCase_CaseInterpolatorList%GetFromKey(OutCaseKey,myDC2)
         if (.not. associated(myDC2)) call runend('Case for a channel is not defined in interpolation')
         call ExtractCaseInterpolator(myDC2,myCaseInterpolator)
         call myCaseInterpolator%GetInterpolatorAndNpoin(InterpolatorPointer,npoinPointer)
         call myInChannel%InitializeInterpolation(InterpolatorPointer,npoinPointer)
      endif
   end subroutine
   
   subroutine SendRunendMessage(CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey,OutDriverKey,OutChannelKey)
      character(DCHashCharSize) :: CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey
      character(DCHashCharSize), optional :: OutDriverKey,OutChannelKey
      
      character(1) :: mychar
      
      character(150) :: errormessage
      
      errormessage = 'Matchchannels, failed to match '//new_line(mychar)//' From: case '//trim(adjustl(CurrentInCaseKey)) 
      errormessage =  trim(adjustl(errormessage))//', driver '//trim(adjustl(CurrentInDriverKey))
      errormessage =  trim(adjustl(errormessage))//', channel '//trim(adjustl(CurrentInChannelKey))
      errormessage = trim(adjustl(errormessage))//''//new_line(mychar)//' To: case  '//trim(adjustl(OutCaseKey))
      if (present(OutDriverkey))  errormessage = trim(adjustl(errormessage))//', driver  '//trim(adjustl(OutDriverKey))
      if (present(OutChannelkey)) errormessage = trim(adjustl(errormessage))//', channel  '//trim(adjustl(OutChannelKey))
      
      call runend(errormessage)
   end subroutine

   subroutine SendOptionalNPMessage(CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey,OutDriverKey,OutChannelKey)
      character(DCHashCharSize) :: CurrentInCasekey, CurrentInDriverKey, CurrentInChannelKey,OutCaseKey
      character(DCHashCharSize), optional :: OutDriverKey,OutChannelKey
      
      character(1) :: mychar
      
      character(250) :: errormessage
      
      errormessage = 'Matchchannels, failed to match '//new_line(mychar)//' From: case '//trim(adjustl(CurrentInCaseKey)) 
      errormessage =  trim(adjustl(errormessage))//', driver '//trim(adjustl(CurrentInDriverKey))
      errormessage =  trim(adjustl(errormessage))//', channel '//trim(adjustl(CurrentInChannelKey))
      errormessage = trim(adjustl(errormessage))//' '//new_line(mychar)//' To: case  '//trim(adjustl(OutCaseKey))
      if (present(OutDriverkey))  errormessage = trim(errormessage)//', driver  '//trim(adjustl(OutDriverKey))
      if (present(OutChannelkey)) errormessage = trim(adjustl(errormessage))//', channel  '//trim(adjustl(OutChannelKey))
      errormessage = trim(errormessage)//', '//new_line(mychar)//' Ignoring OPTIONAL channel  '//trim(adjustl(OutChannelKey))
      if(MPIrank==MPIroot) then
          write(*,901) 'WARNING: '//adjustl(trim(errormessage))
      endif
      901 format(/,5x,a,/)
   end subroutine


end module

subroutine MatchChannels
   use def_master
   use Mod_MatchChannels
   implicit none
   
   
   !This subroutine matches all the channels, even between multiple cases
   
   call CaseList%LoopList(LoopInCases)
   if(MPIrank==MPIroot) write(*,*) '----: Channels matched'
   
end subroutine
