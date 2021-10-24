module Mod_PorousMedia
   use typre
   use Mod_MPIObject
   use Mod_Listen

   implicit none
   type, extends(MPIObject) :: PorousMedia
      real(rp) :: perm1, perm2, minx, maxX, minY, maxY, Dis 
contains
      procedure ReadData
      procedure ReadDataMPI
   end type

contains
   subroutine ReadData(a,Listener)
      class(PorousMedia) :: a
      type(ListenFile) :: Listener
 
      call Listener%Listen('Porous')
      do while (Listener%words(1) /= 'ENDPR')
         if (Listener%words(1) == 'MINX ') then
            a%minX = Listener%param(1)
         elseif (Listener%words(1) == 'MAXX ') then
            a%maxX = Listener%param(1)
         elseif (Listener%words(1) == 'MINY ') then
            a%minY = Listener%param(1)
         elseif (Listener%words(1) == 'MAXY ') then
            a%maxY = Listener%param(1)
         elseif (Listener%words(1) == 'DIST ') then
            a%Dis = Listener%param(1)
         elseif (Listener%words(1) == 'PERM1') then
            a%perm1 = Listener%param(1)   
         elseif (Listener%words(1) == 'PERM2') then
            a%perm2 = Listener%param(1)         
         endif
         call Listener%listen('Porous')   
      enddo
   end subroutine
   
   subroutine ReadDataMPI(a)
      use MPI
      class(PorousMedia) :: a
      integer(ip) :: ierr
      
      CALL MPI_BCAST(a%minX,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%maxX,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%minY,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%maxY,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%Dis,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%perm1,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%perm2,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
   end subroutine
   
end module
