module Mod_DriverInterface
    use typre
    use def_parame
    use Mod_DChash
    use Mod_BroadCastBuffer
    use Mod_MPIObject
    use Mod_CaseVariables
    use Mod_InChannel
    use Mod_DC_InChannel
    use Mod_PhysicalProblem
    implicit none
    private
    
    !This is an abstract Driver.
    !It will be extended either to a physical problem driver
    !Or to a ROM Driver

    public :: DriverInterface
    
    type OutDataType
      integer(8) :: Current,MaxMemo,TotalMemo,TotalMax
      real(rp) :: cputim
    end type


    type, extends(MPIObject), abstract :: DriverInterface

        ! General
        integer(ip) ::    &
            kfl_delay,    &       ! Delay NSI problem
            kfl_conve,    &       ! NSI Convergence required
            kfl_exist,    &       ! Existence of the problem
            kfl_preli,    &       ! Module preliminary run
            kfl_rstar,    &       ! Module restart
            kfl_inter,    &       ! Module interpolated restart
            ndela ,       &       ! Steps to delay the problem
            nprit,           &    ! Preliminary output time frequency
            consec_delay=0          ! Used to delay modules consecutively

        real(rp) ::   cpres= 1.0_rp
        
        !Driver communicator 
        type(DistributedContainerList), pointer :: OutChannelList => null()
        type(DistributedContainerList), pointer :: InChannelList => null()
            
        !Outdata for Storing info to deliver
        type(OutDataType) :: OutData
        integer(ip) :: kfl_algor = 1  !Default algorithm

        
    contains

        procedure                                  :: Initialize
        procedure                                  :: SetDriverMPI
        procedure                                  :: InitializeChannelLists
        procedure                                  :: FinalizeChannelLists  
        procedure                                  :: GetChannelLists
        procedure                                  :: Reapro      !This one is entered by root only, it reads
        procedure                                  :: ReaproMPI   !This one is entered by everyone, it broadcasts the read info
        procedure(Lev1Reapro),            deferred :: Lev1Reapro
        procedure(Lev1ReaproMPI),         deferred :: Lev1ReaproMPI
        procedure(BaseInterface),         deferred :: SetInChannels
        procedure(SetOutChannels),        deferred :: SetOutChannels
        procedure(BaseInterface),         deferred :: UpdateChannels
        
        !These are usually a call to a physical problem
        procedure(UsualInterface),        deferred :: TurnOn
        procedure(UsualInterface),        deferred :: GetRefinementCriteria
        procedure                                  :: PreRefine => NULLSUB
        procedure(UsualInterface),        deferred :: Refine
        procedure(UsualInterface),        deferred :: Rebalance
        procedure(MeshProjections),       deferred :: MeshProjections
        procedure(MeshAdvections),        deferred :: MeshAdvections
        procedure(UsualInterface),        deferred :: SetTimeStep
        procedure(UsualInterface),        deferred :: GetTimeStep
        procedure(UsualInterface),        deferred :: Begste
        procedure(UsualInterface),        deferred :: DoIter
        procedure(Convergence),           deferred :: Convergence
        procedure(UsualInterface),        deferred :: Endste
        procedure(UsualInterface),        deferred :: Output
        procedure(UsualInterface),        deferred :: Turnof 
        
        procedure                                  :: GetRemeshingCriteria
        procedure                                  :: CoupledConvergence
        procedure                                  :: SetCoupledResidual
        procedure                                  :: GetCoupledResidual
       
        procedure                                  :: GetTime
        procedure                                  :: GetRestartFlags
        procedure(BaseInterface),         deferred :: WriteTimes
        procedure                                  :: GetFromInChannelList
        procedure                                  :: AddToInChannelList
        procedure                                  :: GetPhysicalProblem
        procedure                                  :: GetPhysicalAlgorithm

      procedure :: SetDefaultInterpolatedArrayA1

      generic :: SetDefaultInterpolatedArray => SetDefaultInterpolatedArrayA1

    end type DriverInterface

    abstract interface

      subroutine BaseInterface(a)
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
      end subroutine
       
      subroutine UsualInterface(a,c)
         use Mod_caseVariables
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         type(caseVariables) :: c
      end subroutine

      subroutine Lev1Reapro(a,c)
         use Mod_BroadCastBuffer
         use Mod_caseVariables
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         type(caseVariables) :: c
      end subroutine
       
      subroutine Lev1ReaproMPI(a,c,bb)
         use Mod_BroadCastBuffer
         use Mod_caseVariables
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         type(caseVariables) :: c
         type(BroadCastBuffer) :: bb
      end subroutine
       
      subroutine MeshProjections(a,c,itask)
         use typre
         use Mod_CaseVariables
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         type(caseVariables) :: c
         integer(ip) :: itask
      end subroutine
       
       
      subroutine MeshAdvections(a,c,itask)
         use typre
         use Mod_CaseVariables
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         type(caseVariables) :: c
         integer(ip) :: itask
      end subroutine 
       
      subroutine SetOutChannels(a,task)
         use typre
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         character(*) :: task
      end subroutine
       
      subroutine Convergence(a,c,glres)
         use typre
         use Mod_caseVariables
         import DriverInterface
         implicit none
         class(DriverInterface) :: a
         type(caseVariables) :: c
         real(rp) :: glres
      end subroutine
       
     
   end interface

contains
  
      subroutine Initialize(a,c)
         class(DriverInterface) :: a
         type(caseVariables), intent(in) :: c
         
         call a%SetLOutputFile(c%masterVars%lun_memor,c%masterVars%lun_outpu)
         call a%SetInputUnit(c%masterVars%lun_pdata)
         call a%SetMPI(c%mastervars%MPIcomm,c%masterVars%MPIsize,c%masterVars%MPIroot,c%masterVars%MPIrank)
      end subroutine

      subroutine SetDriverMPI(a,c)
         class(DriverInterface) :: a
         type(caseVariables), intent(in) :: c
         
         call a%SetMPI(c%mastervars%MPIcomm,c%masterVars%MPIsize,c%masterVars%MPIroot,c%masterVars%MPIrank)
      end subroutine
   
      subroutine Reapro(a,c,endword)
         class(DriverInterface) :: a
         type(caseVariables), intent(inout) :: c
         character(5) :: endword

         real(rp), pointer     :: param(:) => NULL()
         character(5), pointer :: words(:)
         integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
         integer(ip) :: MyIostat1
         character(150) :: title

         a%kfl_exist= 1                                      ! Problem is: on-off
         a%kfl_conve= 0                                      ! Should prob converge
         a%kfl_rstar= 0                                      ! Restart off
         a%kfl_preli= 0                                      ! Preliminary run off
         a%nprit = 20
         a%kfl_inter= 0                                      ! Preliminary run off
         a%kfl_delay= 0                                      ! Do not delay problem
         a%ndela= 0                                          ! Steps of delay 

         !Let us read
         call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
         call a%Listener%getarrs(words,param,nnpar,nnwor)  
         call a%Listener%listen('RPRODA')
         do while (words(1)/=endword)
            if(words(1)=='DELAY') then
               a%kfl_delay  = 1
               a%ndela      = int(param(1))
            elseif (words(1)=='CONVE') then
               if(a%Listener%exists('YES  ')) a%kfl_conve = 1
            elseif (words(1)=='PRELI') then
               if(a%Listener%exists('ON   ')) a%kfl_preli = 1  
               a%nprit = a%Listener%getint('FREQU',1,'#Preliminary frequency')
            elseif (words(1)=='RESTA') then
               if(a%Listener%exists('INITI')) then
                  a%kfl_rstar = 1
               elseif(a%Listener%exists('INTER')) then
                  a%kfl_inter = 1
                  a%kfl_rstar = 1
                  read(a%Listener%nunit,*,IOStat=myIoStat1) title, c%masterVars%oldnamda
               endif      
            elseif (words(1)=='CONSE') then
               a%consec_delay = int(param(1))
            else
               call a%Lev1Reapro(c)
            endif
            call a%Listener%listen('RPRODA')  
         enddo
       end subroutine
     
       subroutine ReaproMPI(a,c)
         class(DriverInterface) :: a
         type(CaseVariables) :: c
         type(BroadCastBuffer) :: bb
       
         !MPI Communications     
         call bb%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
         call bb%SetLOutputFile(a%lun_memo,a%lun_outpu)
         call bb%Initialize(100,100)

         call bb%Add(a%kfl_exist)
         call bb%Add(a%kfl_conve)
         call bb%Add(a%kfl_delay)
         call bb%Add(a%kfl_preli)
         call bb%Add(a%kfl_rstar)     
         call bb%Add(a%kfl_inter) 
         call bb%Add(a%nprit) 
         call bb%Add(150,c%masterVars%oldnamda)                             
         call bb%Add(a%ndela)
         call bb%Add(a%consec_delay)

         !Particular read
         call a%Lev1ReaproMPI(c,bb)

         call bb%BroadCast
         call bb%Dealloc
      end subroutine

       subroutine InitializeChannelLists(a)
         class(DriverInterface) :: a
         
         allocate(a%OutChannelList,a%InChannelList)
         call a%Memor%allocObj(0,'ChannelLists','InitializeChannelLists',2_ip)
         
         call a%OutChannelList%Initialize
         call a%InChannelList%Initialize
      end subroutine
      
      subroutine FinalizeChannelLists(a)
         class(DriverInterface) :: a

         call a%OutChannelList%Finalize
         call a%InChannelList%Finalize(.true.)
         
         deallocate(a%OutChannelList,a%InChannelList)
         call a%Memor%deallocObj(0,'ChannelLists','InitializeChannelLists',2_ip)
      end subroutine
      
      subroutine GetChannelLists(a,OutChannelList,InChannelList)
         class(DriverInterface) :: a
         type(DistributedContainerList), pointer :: OutChannelList
         type(DistributedContainerList), pointer :: InChannelList
         
         OutChannelList => a%OutChannelList
         InChannelList => a%InChannelList
      end subroutine
      
      subroutine GetTime(a,cputim)
         class(DriverInterface) :: a
         real(rp) :: cputim
         
         cputim = a%OutData%cputim
         
      end subroutine
      
      subroutine GetRestartFlags(a,kfl_rstar,kfl_inter,kfl_preli)
         class(DriverInterface) :: a
         integer(ip) :: kfl_rstar,kfl_inter,kfl_preli
        
         kfl_rstar = a%kfl_rstar
         kfl_inter = a%kfl_inter
         kfl_preli = a%kfl_preli
         
      end subroutine
      
      subroutine GetMemo(a,Current,MaxMemo,TotalMemo,TotalMax)
         implicit none
         class(DriverInterface) :: a
         integer(8) :: Current,MaxMemo,TotalMemo,TotalMax
         
         Current = a%OutData%Current
         MaxMemo = a%OutData%MaxMemo
         TotalMemo = a%OutData%TotalMemo
         TotalMax = a%OutData%TotalMax
      end subroutine
      
      subroutine SetDefaultInterpolatedArrayA1(a,akey,array)
         class(DriverInterface) :: a
         character(*) ::akey
         real(rp) :: array(:)
         
         class(DistributedContainer), pointer :: auxDC => NULL()
         type(InChannel), pointer :: myInChannel => NULL()
         character(DCHashCharSize) :: key
         
         key = akey
         
         call a%InChannelList%GetFromKey(key,auxDC)
         call ExtractInChannel(auxDC,myInChannel)
         call myInChannel%SetDefaultInterpolatedArray(array)
      end subroutine
      
      subroutine GetFromInChannelList(a,akey,myDC)
         class(DriverInterface) :: a
         character(*) ::akey
         class(DistributedContainer), pointer :: myDC
         
         class(DistributedContainer), pointer :: auxDC => NULL()
         type(InChannel), pointer :: myInChannel => NULL()
         character(DCHashCharSize) :: key
         
         key = akey
         
         call a%InChannelList%GetFromKey(key,auxDC)
         call ExtractInChannel(auxDC,myInChannel)
         call myInChannel%GetDC(myDC)
      end subroutine
      
      subroutine AddToInChannelList(a,key,ChannelKey,InterpolateType,InGroupInfo,ChannelPriority)
         class(DriverInterface)   :: a
         character(*)             :: key,ChannelKey,InterpolateType
         character(*),optional    :: ChannelPriority
         type(InChannelGroupInfo) :: InGroupInfo
         
         type(InChannel), pointer :: myInChannel => NULL()
         class(DistributedContainer), pointer :: myDC => NULL()
         
         myInChannel => InChannel_Const(InGroupInfo,adjustl(trim(ChannelKey)),adjustl(trim(InterpolateType)),adjustl(trim(ChannelPriority)))
         myDC => DC_InChannel_Const(myInChannel)
         call myDC%SetKey(key)
         call a%InChannelList%Add(myDC)
      end subroutine

      subroutine CoupledConvergence(a,c,kfl_flag,cpres)
          implicit none
          class(DriverInterface) :: a
          type(caseVariables) :: c
          logical             :: kfl_flag
          real(rp)            :: cpres

          kfl_flag = .true.
          cpres    = 0.0_rp 
      end subroutine

      subroutine SetCoupledResidual(a,cpres)
          implicit none
          class(DriverInterface) :: a
          real(rp)               :: cpres

          a%cpres    = cpres
      end subroutine

      subroutine GetCoupledResidual(a,cpres)
          implicit none
          class(DriverInterface) :: a
          real(rp)               :: cpres

          cpres    = a%cpres
      end subroutine
      
      subroutine GetRemeshingCriteria(a,DoRemesh)
         class(DriverInterface) :: a
         integer(ip) :: DoRemesh
         
         call runend('This driver does not have a criteria for remeshing')
      end subroutine
      
       subroutine GetPhysicalProblem(a,PhysicalPro)
         implicit none
         class(DriverInterface), target :: a
         class(PhysicalProblem), pointer :: PhysicalPro
         
         call runend('DriverInterface, GetPhysicalProblem not implemented')
      end subroutine

       subroutine GetPhysicalAlgorithm(a,PhysicalAlgor)
         implicit none
         class(DriverInterface), target :: a
         integer(ip) :: PhysicalAlgor

         PhysicalAlgor = a%kfl_algor

     end subroutine
      
       subroutine NULLSUB(a,c)
         implicit none
         class(DriverInterface):: a
         type(CaseVariables) :: c
         
      end subroutine

  end module Mod_DriverInterface
