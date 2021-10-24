module Mod_LocalDimensionsAndIRenumbering
   use typre
   use Mod_Mesh
   implicit none

contains
   
   subroutine LocalDimensionsAndIRenumberingSerial(a)
      use MPI
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh), intent(inout) :: a

      integer(ip), allocatable :: aux_ijasizes(:,:)
      
      !MPI
      integer, parameter :: mtag1 = 1, mtag2 = 2, mtag3 = 3
      integer status(MPI_STATUS_SIZE),irequest1(a%MPIsize)
      integer :: ierr
      integer(ip) :: ipoin,irank,ijasizes(2)
      
      !Nodes
      call a%Memor%alloc(a%gnpoin,a%igPointGlobNumber,'igPointGlobNumber','ScatterData')
      if (a%MPIrank == a%MPIroot) then
!         call a%Memor%alloc(a%gnpoin,a%igPointGlobNumber,'igPointGlobNumber','ScatterData')
         call a%Memor%alloc(2,a%MPIsize,aux_ijasizes,'aux_ijasizes','ScatterData')
         !aux_ijasizes contains:  1 npoinLocal
         !                        2 global numbering of the first point
         
         !Inverse numbering 
         do ipoin = 1,a%gnpoin
            a%igPointGlobNumber(a%gPointGlobNumber(ipoin)) = ipoin
            !Local dimensions
            aux_ijasizes(1,a%gPointProcNumber(ipoin)+1) = aux_ijasizes(1,a%gPointProcNumber(ipoin)+1)+1
         enddo
         
         !Send the local dimensions
         aux_ijasizes(2,1) = 0
         do irank = 0,a%MPIsize-1
            call MPI_ISEND(aux_ijasizes(:,irank+1),2,MPI_INTEGER4,irank,mtag1,a%MPIcomm,irequest1(irank+1),ierr)
            if (irank < a%MPIsize-1) aux_ijasizes(2,irank+2) = aux_ijasizes(2,irank+1) + aux_ijasizes(1,irank+1) !Initial point, next proc
         enddo
         
         !Statistics
         a%MaxNpoinPerProc = maxval(aux_ijasizes(1,:))
      endif

      CALL MPI_BCAST(a%igPointGlobNumber,a%gnpoin,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      
      !Receive the Local dimensions 
      call MPI_RECV(ijasizes,2,MPI_INTEGER4,a%MPIroot,mtag1,a%MPIcomm,status,ierr)
      a%npoinLocal = ijasizes(1)
      a%poin0 = ijasizes(2) 
      if (a%MPIrank == a%MPIroot) then
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest1(irank+1), status, ierr)
         enddo   
         call a%Memor%dealloc(2,a%MPIsize,aux_ijasizes,'aux_ijasizes','ScatterData')
      endif
   end subroutine   
   
   subroutine LocalDimensionsAndIRenumberingPartitioned(a)
      use MPI
      use typre
      use Mod_Mesh
      implicit none
      class(FemMesh), intent(inout) :: a
      integer(ip), allocatable :: npoinLocals(:),poin0s(:)
      integer :: ierr
      integer(ip) :: irank
      
      !We need only:
      !poin0
      !npoinLocal
      a%npoinLocal = sum(a%RP_NpoinINeed)
      
      !To be used only by root
      call a%Memor%alloc(a%MPIsize,npoinLocals,'npoinLocals','LocalDimensionsAndIRenumberingPartitioned')
      call a%Memor%alloc(a%MPIsize,poin0s,'poin0s','LocalDimensionsAndIRenumberingPartitioned')

      call MPI_Gather(a%npoinLocal,1,MPI_INTEGER4,npoinLocals,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr); 

      if (a%MPIrank == a%MPIroot) then
         poin0s(1) = 0
         do irank = 0,a%MPIsize-2
            poin0s(irank+2) = poin0s(irank+1)+npoinLocals(irank+1)
         enddo
         a%MaxNpoinPerProc = maxval(npoinLocals)
      endif
      call MPI_Scatter(poin0s,1_ip,MPI_INTEGER4,a%poin0,1_ip,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      
      !Deallocates
      call a%Memor%dealloc(a%MPIsize,npoinLocals,'npoinLocals','LocalDimensionsAndIRenumberingPartitioned')
      call a%Memor%dealloc(a%MPIsize,poin0s,'poin0s','LocalDimensionsAndIRenumberingPartitioned')
      
   end subroutine   
   
end module

subroutine LocalDimensionsAndIRenumbering(a)
   use MPI
   use typre
   use Mod_Mesh
   use Mod_LocalDimensionsAndIRenumbering
   implicit none
   class(FemMesh), intent(inout) :: a
   
   call a%Timer%LocalDimensionsAndIRenumbering%Tic

   !Serial
   if (a%kfl_ReadType == 0) then
      call LocalDimensionsAndIRenumberingSerial(a)
   
   !Partitioned
   elseif (a%kfl_ReadType == 1) then
      call LocalDimensionsAndIRenumberingPartitioned(a)
   endif
   
   call a%Timer%LocalDimensionsAndIRenumbering%Toc
end subroutine
