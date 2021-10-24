!> This module contains the MPIObject type. It contains the MPI information, plus additional common info
module Mod_MpiObject
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_int2str
   use MPI
   implicit none
   private
   public MpiObject, MPIObjectLight
   
   type MPIObjectLight
      integer(ip)    :: MPIrank
      integer(ip)    :: MPIroot
      integer(ip)    :: MPIsize
      integer(ip)    :: MPIcomm
      
contains      
      procedure :: SetMPI
   end type
   
   
   !> This is the MPI object
   type, extends(MPIObjectLight) ::  MpiObject
      !MPI


      integer(ip)    :: lun_pdata         !Logical unit for reading data
      integer(ip)    :: lun_outpu
      integer(ip)    :: lun_memo
      character(150) :: namda             !FileName for reading domain data
      character(150) :: InputFolder = 'data', OutputFolder = 'results'
      integer(ip)    :: kfl_flush = 0
      integer(ip)    :: kfl_iofor = 0
      integer(ip)    :: kfl_ReadType = -1      !0: SERIAL 1: PARTITIONED
      
      !Auxiliary
      type(MemoryMan)  :: Memor
      type(ListenFile) :: Listener
      
      
contains

      procedure :: SetReadType
      procedure :: SetInputFolder
      procedure :: SetInputFile
      procedure :: SetInputUnit
      
      procedure :: SetOutputFolder
      procedure :: SetOutputFile
      procedure :: SetLOutputFile
      generic   :: SetOutputFiles => SetOutputFile, SetLOutputFile
      procedure :: SetFlush
      procedure :: SetIOFormat
      
      procedure :: GetMemo

   end type
   
contains
   
   subroutine SetMPI(a,comm,size,root,irank)
      implicit none
      class(MpiObjectLight) :: a
      integer(ip) :: comm,irank,size,root
      
      a%MPIcomm = comm
      a%MPIrank = irank
      a%MPIsize = size
      a%MPIroot = root
      a%MPIcomm = comm
   
   
   end subroutine
   
   !> Sets the input file
   !> @param namda The input file name
   
   subroutine SetReadType(a,ReadTypeString)
      class(MPIObject) :: a
      character(*) :: ReadTypeString
      
      if (trim(ReadTypeString) == 'SERIAL') then
         a%kfl_ReadType = 0
      elseif (trim(ReadTypeString) == 'PARTITIONED') then
         a%kfl_ReadType = 1
      else
         call runend('MPIObject: Read Type unknown')
      endif
   end subroutine
   
   subroutine SetInputFolder(a,InputFolder) 
      use Mod_Listen
      implicit none
      class(MPIObject) :: a
      character(150) :: InputFolder
      if (a%kfl_ReadType == -1) then
         call runend('MPIObject: setting input folder before setting read type')
      endif
      
      if (a%kfl_ReadType == 0) then
         a%InputFolder = InputFolder
      elseif (a%kfl_ReadType == 1) then
         a%InputFolder = trim(InputFolder)//'/data'//int2str(a%MPIrank)
      else 
      endif
      
      call a%Listener%SetReadingFolder(a%InputFolder)
   end subroutine
   
   subroutine SetInputFile(a,namda)
      implicit none
      class(MpiObject) :: a
      character(150) :: namda
      
      a%namda = namda
   end subroutine
   
   subroutine SetInputUnit(a,lun_pdata)
      implicit none
      class(MPIObject) :: a
      integer(ip) :: lun_pdata
      
      a%lun_pdata = lun_pdata
   end subroutine
   
   subroutine SetLOutputFile(a,lun_memo,lun_outpu)
      implicit none
      class(MpiObject) :: a
      integer(ip)    :: lun_memo, lun_outpu
      
      a%lun_outpu = lun_outpu
      a%lun_memo = lun_memo
      call a%Memor%init(lun_memo,lun_outpu)
   
   end subroutine
   
   subroutine SetFlush(a,kfl_flush)
      implicit none
      class(MPIObject) :: a
      integer(ip) :: kfl_flush
      
      a%kfl_flush = kfl_flush
   end subroutine
  
   subroutine SetIOFormat(a,kfl_iofor)
      implicit none
      class(MPIObject) :: a
      integer(ip) :: kfl_iofor
      
      a%kfl_iofor = kfl_iofor
   end subroutine
  
   subroutine SetOutputFolder(a,OutputFolder)
      implicit none
      class(MPIObject) :: a
      character(150) :: OutputFolder
      a%OutputFolder = OutputFolder
   end subroutine
   
   subroutine SetOutputFile(a,namdamemo,namdaoutpu)
      implicit none
      class(MpiObject) :: a
      character(150) :: namdamemo,namdaoutpu
      
      call runend('Initb for type Mesh to be defined')
   end subroutine
   
   subroutine GetMemo(a,Current,MaxMemo,TotalMemo,TotalMax)
      implicit none
      class(MpiObject) :: a
      integer(8) :: Current,MaxMemo,TotalMemo,TotalMax
      
      call a%Memor%GetValue(Current,MaxMemo,TotalMemo,TotalMax)
   
   
   end subroutine
   
   
   
   
end module
