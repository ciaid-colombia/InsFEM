module Mod_ExitBodyForces
   use typre
   use Mod_MPIObject
   use Mod_Listen
   implicit none
   
   type, extends(MPIObject) :: ExitBodyForce
      real(rp) :: tx0, txsc,tampt, ObjectiveVelocity(3),density

      
contains
      procedure ReadData
      procedure ReadDataMPI
      procedure GetForce
      procedure SetDensity

   end type

contains
   subroutine ReadData(a,Listener)
      class(ExitBodyForce) :: a
      type(ListenFile) :: Listener
   
      call Listener%Listen('Turbul')
      do while (Listener%words(1) /= 'ENDEX')
         if (Listener%words(1) == 'TX0  ') then
            a%tx0 = Listener%param(1)
         elseif (Listener%words(1) == 'TXSC ') then
            a%txsc = Listener%param(1)
         elseif (Listener%words(1) == 'TAMPT') then
            a%tampt = Listener%param(1)
         elseif (Listener%words(1) == 'OBJEC') then
            a%ObjectiveVelocity(1) = Listener%param(1)   
            a%ObjectiveVelocity(2) = Listener%param(2)   
            a%ObjectiveVelocity(3) = Listener%param(3)   
         endif
         call Listener%listen('Turbul')   
      enddo
   end subroutine
   
   subroutine ReadDataMPI(a)
      use MPI
      class(ExitBodyForce) :: a
      
      integer(ip) :: ierr
      
      CALL MPI_BCAST(a%tx0,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%txsc,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%tampt,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%ObjectiveVelocity,3,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
      
      !We do the initialization here
      
      
   end subroutine
   
   subroutine SetDensity(a,density)
      class(ExitBodyForce) :: a
      real(rp) :: density
      a%density = density
   end subroutine
   
   subroutine GetForce(a,ctime,coord,gpvel,force)
      class(ExitBodyForce) :: a
      real(rp) :: ctime, coord(3), force(3),gpvel(3)
      
      integer(ip) :: i
      real(rp) :: p, b,r, f2, fzt, x, y, z
      
      force = 0.0_rp
      
      x = coord(1)
      y = coord(2)
      z = coord(3)
      
      force(:) = force(:) + a%density*a%tampt*exp(-((x-a%tx0)/a%txsc)**2)*(a%ObjectiveVelocity-gpvel)
   end subroutine   
   
   


end module