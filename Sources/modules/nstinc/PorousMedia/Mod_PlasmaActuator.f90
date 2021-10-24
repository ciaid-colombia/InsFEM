module Mod_PlasmaActuator 
   use typre
   use Mod_MPIObject
   use Mod_Listen
   implicit none

   type, extends(MPIObject) :: PlasmaAct
      
      integer(ip) :: nboxes = 1
      real(rp) :: boxcoord(3,10), boxdirection(3,10) = 1.0_rp
      real(rp) :: XforceConstant = 2.88e-3, yForceConstant = 1.44e-3_rp
      real(rp) :: VoltageInKV = 4.0_rp, DDistanceInm = 0.00025, BDistanceInm = 0.003,CDIstanceInm=0.0015
      real(rp) :: FrequencyInKHZ = 3.0_rp
      real(rp) :: Eb_KV_m =3000
      real(rp) :: k1,k2,E0
      
 contains
       procedure ReadData
       procedure ReadDataMPI
       Procedure PlasmaForce
 
    end type
 
  contains
    subroutine ReadData(a,Listener)
       class(PlasmaAct) :: a
       type(ListenFile) :: Listener
       integer(ip)      :: ibox

 
      call Listener%Listen('Plasma')
      do while (Listener%words(1) /= 'ENDPL')
         if (Listener%words(1) == 'NBOX ') then
            a%nboxes = Listener%param(1)
            do ibox = 1,a%nboxes
               a%boxcoord(:,ibox) = Listener%param(2+(ibox-1)*6:1+ibox*6-1)
               a%boxdirection(:,ibox) = Listener%param(2+(ibox-1)*6+3:1+ibox*6-1+3)
            enddo
         elseif (Listener%words(1) == 'DDIST') then
            a%DDistanceInm = Listener%param(1)
         elseif (Listener%words(1) == 'BDIST') then
            a%BDistanceInm = Listener%param(1)
         elseif (Listener%words(1) == 'CDIST') then
            a%CDistanceInm = Listener%param(1)
         elseif (Listener%words(1) == 'VOLTA') then
            a%VoltageInKV = Listener%param(1)
         elseif (Listener%words(1) == 'FREQU') then   
            a%FrequencyInKHZ = Listener%param(1)
         endif
         call Listener%listen('Plasma')   
      enddo
   end subroutine
   
   subroutine ReadDataMPI(a)
      use MPI
      class(PlasmaAct) :: a
      
      integer(ip) :: ierr
      
      CALL MPI_BCAST(a%nboxes,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%boxcoord,30,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%boxdirection,30,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%VoltageInKv,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%FrequencyInKHZ,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%DDistanceInm,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%BDistanceInm,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%CDistanceInm,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
      a%E0 = a%VolTageInKV/a%DDistanceInm
      a%k1 = (a%E0-a%Eb_KV_m)/a%BDistanceInm
      a%k2 = (a%E0-a%Eb_KV_m)/a%CDistanceInm
   end subroutine

   subroutine PlasmaForce(a,gpcod,Pforce)
      class(PlasmaAct) :: a
      real(rp)         :: Pforce(3), gpcod(3)
      real(rp)         :: E, freqcorrect
      integer(ip)      :: ibox
      
      do ibox = 1,a%nboxes
         E = a%E0 - sqrt((a%k1*abs(gpcod(1)-a%boxcoord(1,ibox)))**2 + (a%k2*abs(gpcod(2)-a%boxcoord(2,ibox)))**2)
         if (E > a%Eb_KV_m) then
            !Frequency correction factor
            freqcorrect = a%FrequencyInKHZ/3.0_rp
            Pforce(1) = Pforce(1) + E*a%XforceConstant*freqcorrect*a%boxdirection(1,ibox)
            Pforce(2) = PForce(2) + E*a%XForceConstant*freqcorrect*a%boxdirection(2,ibox)
            Pforce(3) = PForce(3) + E*a%XForceConstant*freqcorrect*a%boxdirection(3,ibox)
         endif
      enddo
      
   end subroutine 
end module
