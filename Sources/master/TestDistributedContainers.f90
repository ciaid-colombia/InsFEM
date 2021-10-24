subroutine TestDistributedContainers
   use typre
   use Mod_DCHashCharSize
   use Mod_DistributedContainer
   use Mod_DistributedContainerList
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   use Mod_DC_i1
   use Mod_DC_Physical
   use Mod_InChannel
   use Mod_DriverInterface
   use Mod_PhysicalProblemDriver
   use Mod_NstincDriver
   use Mod_GeneralCase
   
   implicit none
   
   type (DistributedContainerList), target :: DCList
   type (DistributedContainerList) :: InChannelList
   type (DistributedContainerList) :: InStringList
   type (DistributedContainerList), pointer :: OutChannelList => NULL()
   

   
   class(DistributedContainer), pointer :: myDC => NULL(), theirsDC, myInDC,myPointedDC
   
   type(DistributedContainerList), pointer :: DCListPointer

   
   
   
   
   integer(ip) :: iparray(5), iparray2(3)
   integer(ip), pointer :: ippointer(:)
   
   class(PhysicalProblem), pointer :: PhysicalPointer
   type(NavierStokesProblem) :: NavierStokes
   character(DCHashCharSize) :: namod

   call DCList%Initialize
   
   iparray = 0
   iparray(1) = 23
   iparray(2) = 35
   
   iparray2 = 0
   iparray2(1) = 423
   iparray2(3) = 234
   
   myDC => DC_i1(iparray)
   call myDC%SetKey('dc1   ')
   call DCList%Add(myDC)
   
   call NavierStokes%SetExmod
   myDC => DC_Physical(NavierStokes)
   call myDC%SetKey('dc2   ')
   call DCList%Add(myDC)
   
   myDC => DC_i1(iparray2)
   call myDC%SetKey('dc3   ')
   call DCList%Add(myDC)
   
   call DCList%GoToFirst
   call DCList%GetNext(myDC)
   do while (associated(myDC))
   
      select type (myDC)
      
      type is (DC_i1)
         call myDC%GetValue(ippointer)
         
         write(*,*) ippointer
         
      type is (DC_Physical)
         call myDC%GetValue(PhysicalPointer)
         
         !should not access this since it should be private
         !just for testing
         write(*,*) PhysicalPointer%exmod
      end select
   
      call DCList%GetNext(myDC)
   enddo
   
   call DCList%GetFromKey('dc2   ',myDC)
   write(*,*) 'dc2   '
   select type (myDC)
   type is (DC_i1)
      write(*,*) myDC%a
   type is (DC_Physical)
      write(*,*) myDC%a%exmod
   end select
   
   call DCList%GetFromKey('dc1   ',myDC)
   write(*,*) 'dc1   '
   select type (myDC)
   type is (DC_i1)
      write(*,*) myDC%a
   type is (DC_Physical)
      write(*,*) myDC%a%exmod
   end select 
   
   call DCList%GetFromKey('dc3   ',myDC)
   write(*,*) 'dc3   '
   select type (myDC)
   type is (DC_i1)
      write(*,*) myDC%a
   type is (DC_Physical)
      write(*,*) myDC%a%exmod
   end select  
   
   
   
   !----------------------------------------------------------------------
   !Testing the channels
   
   call InStringList%Initialize

   
   myDC => DistributedContainer()
   call myDC%SetKey('dc3   ')
   
   call InStringList%Add(myDC)
   
   myDC => DistributedContainer()
   call myDC%SetKey('dc1   ')
   call InStringList%Add(myDC)

   myDC => DistributedContainer()
   call myDC%SetKey('dc2   ')
   call InStringList%Add(myDC)

   !We will store the InContainers here
   call InChannelList%Initialize
   
   !We use as the OutContainers the ones previously defined
   OutChannelList => DCList
   
   !We loop through the strings
   call InStringList%GoToFirst
   call InStringList%GetNext(MyDC)
   do while (associated(myDC))
      call OutChannelList%GetFromKey(myDC%key,theirsDC)
      if (associated(theirsDC)) then
         call InChannelList%Add(theirsDC)
      endif
  
      call InStringList%GetNext(myDC)
   enddo
   
   !The string list is no longer needed
   call InStringList%Finalize
   
   !Now we loop through the incomming channel list
   !We loop through the strings
   call InChannelList%GoToFirst
   call InChannelList%GetNext(myDC)
   do while (associated(myDC))
      
      select type (myDC)
         type is (DC_i1)
      
         type is (DC_Physical)
      
      end select
      
      call InChannelList%GetNext(myDC)
   enddo
   
   call InChannelList%Finalize(.true.)
   call DCList%Finalize
   
   
   
!   call mycommunicator1%Initialize
!   call mycommunicator2%Initialize
!   
!   call mycommunicator1%AddInChannelName('dc1')
!   call mycommunicator1%AddInChannelName('dc3')
!   
!   myDC => DC_i1(iparray)
!   call myDC%SetKey('dc1   ')
!   call mycommunicator2%AddOutChannel(myDC)
!   
!   call NavierStokes%SetExmod
!   myDC => DC_Physical(NavierStokes)
!   call myDC%SetKey('dc2   ')
!   call mycommunicator2%AddOutChannel(myDC)
!   
!   myDC => DC_i1(iparray2)
!   call myDC%SetKey('dc3   ')
!   call mycommunicator2%AddOutChannel(myDC)
!   
!   
!   call mycommunicator1%LinkToOutChannel(mycommunicator2%OutChannelList)
!   
!   call mycommunicator1%GetInChannelList(DCListPointer)
!   
!   call DCListPointer%GoToFirst
!   myDC => DCListPointer%GetNext()
!   do while (associated(myDC))
!      write(*,*) 'Testing an InChannel'
!      write(*,*) 'My key: ', myDC%key 
!      
!      select type (myDC)
!         type is (DC_i1)
!      
!            write(*,*) 'Integer value: ', myDC%a
!      
!         type is (DC_Physical)
!      
!            write(*,*) 'Physical value: ', myDC%a%exmod
!      
!      end select
!      write(*,*)
!         
!      
!      myDC => DCListPointer%GetNext()
!   enddo
!   
   
    
  
end subroutine
