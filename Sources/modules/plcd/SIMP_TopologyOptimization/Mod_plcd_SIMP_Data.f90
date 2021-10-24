module Mod_plcd_SIMPData
   use typre
   use Mod_MPIObject
   use Mod_Listen
   use Mod_Element
   implicit none

   type, extends(MPIObject) :: SIMPDataType
      real(rp) :: p = 4.0_rp
      real(rp) :: eta = 0.5_rp   !Default value 0.5  -> larger -> more dynamic, less safe
      real(rp) :: xi = 0.2_rp    !Default value 0.2  -> larger -> more dynamic, less safe
      real(rp) :: VolumeFraction = 0.2_rp
      real(rp) :: chimin = 1e-5_rp, chimax = 1.0_rp
      integer(ip) :: ChiSpaceChoice = 1 !0: ElementWise, 1: PointWise


      !Auxiliary to be used during the iterative process
      !secant method for the volume restriction condition
      real(rp) :: HistoryMeanChi(3), HistoryObjectiveVolumeFraction(3)

      real(rp), pointer :: NodalChi(:) => NULL()

contains
      procedure :: ReadData
      procedure :: ScatterData
    end type

contains

   subroutine ReadData(a,Listener)
   implicit none
      class(SIMPDataType) :: a
      type(ListenFile) :: Listener

      do while(Listener%words(1)/='ENDTO')
         call Listener%listen('plcd_reaphy')
         if(Listener%words(1) == 'P    ') then
            a%p = Listener%param(1)
         elseif(Listener%words(1) == 'ETA  ')   then
            a%eta = Listener%param(1)
         elseif(Listener%words(1) == 'XI   ')   then
            a%xi = Listener%param(1)
         elseif(Listener%words(1) == 'VOLUM')   then
            a%VolumeFraction = Listener%param(1)
         elseif(Listener%words(1) == 'CHIMI')   then
            a%chimin = Listener%param(1)
         elseif(Listener%words(1) == 'CHIMA')   then
            a%chimax = Listener%param(1)
         elseif(Listener%words(1) == 'CHISP')   then
            a%ChiSpaceChoice = Listener%param(1)
         endif
       enddo
       !This is just so that in the external loop there are no problems
       Listener%words(1) = '     '
   end subroutine

   subroutine ScatterData(a)
      use MPI
      implicit none
      class(SIMPDataType) :: a

      integer(ip) :: ierr

      call MPI_BCAST(a%p,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%eta,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%xi,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%VolumeFraction,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%chimin,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%chimax,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%ChiSpaceChoice,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
   end subroutine


end module
