module Mod_InChannel
    use typre
    use Mod_DChash
    use Mod_DistributedContainerList
    use Mod_DistributedContainer
    use Mod_MeshInterpolator
    use Mod_DC_rp
    use Mod_Memor
    use Mod_Listen
    use Mod_BroadCastBuffer

    implicit none

    private

    public ::  InChannel, InChannelGroupInfo, InChannelGroupInfoAddToBuffer, InChannelGroupInfoFromListener, InChannel_Const
    
    !InChannelGroupInfo
    type :: InChannelGroupInfo
      logical :: IsOn = .false.
      character(DCHashCharSize) :: ChannelKey="",DriverKey="",DriverNick="",CaseKey=""
    end type
    
    !InChannel
    type :: InChannel
      character(DCHashCharSize)   :: ChannelKey
      character(DCHashCharSize)   :: DriverKey 
      character(DCHashCharSize)   :: DriverNick
      character(DCHashCharSize)   :: CaseKey
      character(DCHashCharSize)   :: Priority
      class(DistributedContainer), pointer :: DC => NULL()
      
      integer(ip) :: kfl_InterpolateType
      
      logical :: isInterpolated = .false.
      type(Interpolator), pointer :: Interpolator => NULL()
      integer(ip), pointer :: npoinInterpolator => NULL()
      
      type(DC_rp) :: InterpolatedDC
      
      contains
      procedure :: InitializeInterpolation
      procedure :: SetupInterpolation
      procedure :: UpdateInterpolation
      procedure :: FinalizeInterpolation
      procedure :: GetDC
      procedure :: GetDriverKey
      procedure :: GetDriverNick
      procedure :: GetCaseKey
      procedure :: GetChannelKey
      procedure :: GetChannelType
      procedure :: SetDefaultInterpolatedArrayA1
      
      generic :: SetDefaultInterpolatedArray => SetDefaultInterpolatedArrayA1
      
    end type
    
    interface InChannel_Const
      procedure InChannelConstructor
    end interface

contains

   !-----------------------------------------------------------------------
   !InChannelGroupInfo
   subroutine InChannelGroupInfoFromListener(Listener,InInfo)
      type(InChannelGroupInfo) :: InInfo
      type(ListenFile) :: Listener
      character(DCHashCharSize) :: Casekey
      
      if (.not. Listener%exists('OFF  ')) then
          InInfo%IsOn = .true.
          InInfo%DriverKey  = Listener%words(2)
          InInfo%Casekey    = Listener%words(4)
          InInfo%DriverNick = Listener%words(6)

          if (Listener%words(7) == "DRIVE") then
             InInfo%DriverKey  = Listener%words(8)
          endif
          if (InInfo%DriverNick == "") then
              InInfo%DriverNick = InInfo%DriverKey
          endif
      endif
   end subroutine
   
   subroutine InChannelGroupInfoAddToBuffer(InInfo,bb)
      type(InChannelGroupInfo) :: InInfo
      type(BroadCastBuffer) :: bb
      
      call bb%Add(InInfo%IsOn)
      call bb%Add(DCHashCharSize,InInfo%DriverKey)
      call bb%Add(DCHashCharSize,InInfo%DriverNick)
      call bb%Add(DCHashCharSize,InInfo%CaseKey)
   end subroutine
   
   !-----------------------------------------------------------------------------------------  
   !InChannel
   function InChannelConstructor(InInfo,ChannelKey,InterpolationKey,ChannelPriority)
      type(InChannel), pointer     :: InChannelConstructor
      type(InChannelGroupInfo) :: InInfo
      character(*) :: InterpolationKey, ChannelKey
      character(*),optional :: ChannelPriority
      
      allocate(InChannelConstructor)
      
      InChannelConstructor%ChannelKey = ChannelKey
      InChannelConstructor%DriverKey  = InInfo%DriverKey
      InChannelConstructor%DriverNick = InInfo%DriverNick
      InChannelConstructor%CaseKey    = InInfo%CaseKey

      if(present(ChannelPriority)) then
          InChannelConstructor%Priority   = ChannelPriority
      else
          InChannelConstructor%Priority   = "MAIN"
      endif


      
      if (InterpolationKey == 'INTERPOLABLE') then
         InChannelConstructor%kfl_InterpolateType = 1
      elseif (InterpolationKey == 'DO_NOT_INTERPOLATE') then
         InChannelConstructor%kfl_InterpolateType = 2
      elseif (InterpolationKey == 'SKIP_IF_INTERPOLATE') then
         InChannelConstructor%kfl_InterpolateType = 3
      elseif (InterpolationKey == 'RUNEND_IF_INTERPOLATE') then
         InChannelConstructor%kfl_InterpolateType = 4
      else
         call runend ('Mod_InChannel::InChannelConstructor: wrong interpolation key')
      endif

   end function
   
   subroutine GetCaseKey(a,CaseKey)
      class(InChannel) :: a
      character(DCHashCharSize) :: Casekey
      
      CaseKey = a%CaseKey
   end subroutine

   subroutine GetDriverNick(a,DriverNick)
      class(InChannel) :: a
      character(DCHashCharSize) :: DriverNick
      
      DriverNick = a%DriverNick
   end subroutine
   
   subroutine GetDriverKey(a,DriverKey)
      class(InChannel) :: a
      character(DCHashCharSize) :: Driverkey
      
      DriverKey = a%DriverKey
   end subroutine
   
   subroutine GetChannelKey(a,ChannelKey)
      class(InChannel) :: a
      character(DCHashCharSize) :: ChannelKey
      
      ChannelKey = a%ChannelKey
   end subroutine
   
      
   subroutine InitializeInterpolation(a,ChannelInterpolator,npoinInterpolator)
      class(InChannel) :: a
      type(Interpolator), target :: ChannelInterpolator
      integer(ip), pointer :: npoinInterpolator
   
      if (a%kfl_interpolateType == 1) then
         a%Interpolator => ChannelInterpolator
         a%npoinInterpolator => npoinInterpolator
         a%isInterpolated = .true.
      elseif (a%kfl_interpolateType == 2) then
         a%isInterpolated = .false.
      elseif (a%kfl_interpolateType == 3) then
         a%isInterpolated = .false.
         !call runend('Dont know how to manage skip_if_interpolable yet')   
      
      elseif (a%kfl_interpolateType == 4) then
         call runend('Mod_InChannel::InitializeInterpolation: Trying to interpolate a non-interpolable type')
      endif
   end subroutine
      
   subroutine SetupInterpolation(a,Memor)
      class(InChannel) :: a
      type(MemoryMan) :: Memor
      
      class(DistributedContainer), pointer :: DC => NULL()
      
      if (a%isInterpolated .eqv. .false.) then
         !Nothing to do if it is not inerpolated
      else
         DC => a%DC
         select type (DC)
         type is (DC_rp)
            if (associated(DC%a1)) then
               if (associated(a%InterpolatedDC%a1)) then
                  call Memor%pdealloc(size(a%InterpolatedDC%a1,1),a%InterpolatedDC%a1,'InterpolatedDC','InChannel')
               endif
               call Memor%palloc(a%npoinInterpolator,a%InterpolatedDC%a1,'InterpolatedChannel','InChannel')
            elseif (associated(DC%a2)) then
               if (associated(a%InterpolatedDC%a2)) then
                  call Memor%pdealloc(size(a%InterpolatedDC%a2,1),size(a%InterpolatedDC%a2,2),a%InterpolatedDC%a2,'InterpolatedChannel','InChannel')
               endif
               call Memor%palloc(size(DC%a2,1),a%npoinInterpolator,a%InterpolatedDC%a2,'InterpolatedChannel','InChannel')
            else
               call runend('Mod_InChannel::SetupInterpolation: Error connecting type of channel')   
            endif
         class default
            call runend('Mod_InChannel::SetupInterpolation: Wrong DC type for interpolated inChannel')
         end select
      endif
   end subroutine   
   
   subroutine UpdateInterpolation(a)
      class(InChannel) :: a
      
      class(DistributedContainer), pointer :: DC => NULL()
      integer(ip) :: icomp
      
      if (a%isInterpolated .eqv. .false.) then
         !Nothing to do if it is not interpolated
      else 
         
         DC => a%DC
         
         select type (DC)
         type is (DC_rp) 
            if (associated(DC%a1)) then
               call a%Interpolator%Interpolate(1,DC%a1,a%InterpolatedDC%a1)
            elseif (associated(DC%a2)) then
               call a%Interpolator%Interpolate(size(DC%a2,1),DC%a2,a%InterpolatedDC%a2)
            else
               call runend('Mod_InChannel::UpdateInterpolation: Error updating channel, Type not found')   
            endif
         class default
            call runend('Mod_InChannel::UpdateInterpolation: Wrong DC type')   
         end select
      endif
   end subroutine
   
   subroutine GetDC(a,myDC)
      class(InChannel), target :: a
      class(DistributedContainer), pointer :: myDC
      
      if (a%isInterpolated .eqv. .false.) then
         myDC => a%DC
      else
         call a%UpdateInterpolation !maybe needs to be thought
         myDC => a%InterpolatedDC
      endif
   end subroutine
   
   subroutine SetDefaultInterpolatedArrayA1(a,array)
      class(InChannel) :: a
      real(rp) :: array(:)
      
      class(DistributedContainer), pointer :: DC => NULL()
      
      if (a%isInterpolated .eqv. .false.) then
         !Nothing to do if it is not interpolated
      else 
         DC => a%DC
         select type (DC)
         type is (DC_rp) 
            if (associated(DC%a1)) then
               a%InterpolatedDC%a1 = array
            elseif (associated(DC%a2)) then
               call runend('Mod_InChannel::SetDefaultInterpolatedArrayA1: Wrong type of array passed or whatever')
            else
               call runend('Mod_InChannel::SetDefaultInterpolatedArrayA1: Error Setting Default Interpolated Array, Type not found')   
            endif
         class default
            call runend('Mod_InChannel::SetDefaultInterpolatedArrayA1: Wrong DC type')     
         end select
      endif
   end subroutine
   
   subroutine GetChannelType(a,kfl_interpolateType)
      class(InChannel) :: a
      integer(ip) :: kfl_interpolateType
      
      kfl_interpolateType = a%kfl_interpolateType
   end subroutine
   
   subroutine FinalizeInterpolation(a,Memor)
      class(InChannel) :: a
      type(MemoryMan) :: Memor
      
      if (a%isInterpolated .eqv. .false.) then
         !Nothing to do if it is not inerpolated
      else
         if (associated(a%InterpolatedDC%a1)) then
            call Memor%pdealloc(size(a%InterpolatedDC%a1,1),a%InterpolatedDC%a1,'InterpolatedDC','InChannel')
         elseif (associated(a%InterpolatedDC%a2)) then
            call Memor%pdealloc(size(a%InterpolatedDC%a2,1),size(a%InterpolatedDC%a2,2),a%InterpolatedDC%a2,'InterpolatedDC','InChannel')
         else
            call runend('Mod_InChannel::FinalizeInterpolation: Error Setting Finalizing Interpolation, Type not found')   
         endif
      endif      
   end subroutine  
      
end module 
