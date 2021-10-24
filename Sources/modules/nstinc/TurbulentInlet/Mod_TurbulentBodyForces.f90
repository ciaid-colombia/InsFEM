module Mod_TurbulentBodyForces
   use typre
   use Mod_MPIObject
   use Mod_Listen
   implicit none
   
   integer(ip), parameter :: nfunctions = 2
   
   real(rp) :: ftmean = 0.0_rp, ftmax=-1e6, ftmin=1e6
      integer(ip) :: nftmean = 0
   
   type, extends(MPIObject) :: TurbulentInletForce
      real(rp) :: tx0, ty0, txsc, tysc,tampt,density,tdt = 1e-6_rp,zdt = 1e-6_rp,xwidth=1e6_rp,ywidth = 1e6_rp
      
      real(rp) :: p
      integer(ip) :: i
      
      real(rp) :: FourierSeriesFrequencies(3,nfunctions), FourierSeriesPhase(nfunctions)
      real(rp) :: minmax(2) = 0.0_rp
      
      integer(ip) :: kfl_isRandom = 0
      
      
contains
      procedure ReadData
      procedure ReadDataMPI
      procedure GetForce
      procedure SetDensity
      procedure ResetRandomSeed

   end type

contains
   subroutine ReadData(a,Listener)
      class(TurbulentInletForce) :: a
      type(ListenFile) :: Listener
   
      call Listener%Listen('Turbul')
      do while (Listener%words(1) /= 'ENDTU')
         if (Listener%words(1) == 'TX0  ') then
            a%tx0 = Listener%param(1)
         elseif (Listener%words(1) == 'TY0  ') then
            a%ty0 = Listener%param(1)   
         elseif (Listener%words(1) == 'TXSC ') then
            a%txsc = Listener%param(1)
         elseif (Listener%words(1) == 'TYSC ') then
            a%tysc = Listener%param(1)         
         elseif (Listener%words(1) == 'TAMPT') then
            a%tampt = Listener%param(1)
         elseif (Listener%words(1) == 'TDT  ') then
            a%tdt = Listener%param(1)   
         elseif (Listener%words(1) == 'ZDT  ') then
            a%zdt = Listener%param(1)
         elseif (Listener%words(1) == 'XWIDT') then
            a%xwidth = Listener%param(1)
         elseif (Listener%words(1) == 'YWIDT') then
            a%ywidth = Listener%param(1)   
         elseif (Listener%words(1) == 'RANDO') then
            if (Listener%exists('ON   ')) a%kfl_isRandom = 1
         endif
         call Listener%listen('Turbul')   
      enddo
   end subroutine
   
   subroutine ReadDataMPI(a)
      use MPI
      class(TurbulentInletForce) :: a
      
      integer(ip) :: ierr,i,j
      
      CALL MPI_BCAST(a%tx0,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%ty0,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%txsc,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%tysc,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%tampt,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%tdt,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%zdt,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%xwidth,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%ywidth,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%kfl_isRandom,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      
      
      !We do the initialization here
      do i = 1,nfunctions
         do j = 1,3
            call random_number(a%FourierSeriesFrequencies(j,i))
            a%FourierSeriesFrequencies(j,i) = (a%FourierSeriesFrequencies(j,i) + 0.5)*2.0_rp*3.14_rp/real(i)
         enddo
         call random_number(a%FourierSeriesPhase(i))
      enddo
      a%FourierSeriesPhase=a%FourierSeriesPhase*2.0_rp*3.14_rp
      
      write(*,*) 'FourierSeriesFrequencies: ', a%FourierSeriesFrequencies
      write(*,*) 'FourierSeriesPhase: ', a%FourierSeriesPhase
      write(*,*) 'a%xwidth: ', a%xwidth, a%ywidth
      
      
   end subroutine
   
   subroutine SetDensity(a,density)
      class(TurbulentInletForce) :: a
      real(rp) :: density
      a%density = density
   end subroutine
   
   subroutine ResetRandomSeed(a,istep)
      class(TurbulentInletForce) :: a
      integer(ip) :: istep
      integer(ip) :: getsize,i
      integer(ip),allocatable :: seed(:)
      
      call RANDOM_SEED(size=getsize)
      
      allocate(seed(getsize))
      do i = 1,getsize
         seed(i) = istep+i
      enddo
      call RANDOM_SEED(put=seed)
      deallocate(seed)
      
      if (nftmean /= 0) ftmean = ftmean / nftmean
      write(*,*) 'ftmean = ', ftmean
      write(*,*) 'ftmax = ', ftmax
      write(*,*) 'ftmin = ', ftmin
      ftmean = 0.0_rp
      ftmax = -1e6_rp
      ftmin = 1e6_rp
      nftmean = 0
   end subroutine
   
   subroutine GetForce(a,ctime,coord,force)
      class(TurbulentInletForce) :: a
      real(rp) :: ctime, coord(3), force(3)
      
      integer(ip) :: i,getsize,j
      real(rp) :: p, b, f2, fcoordt(3),ft, x, y, z,timeaux,spaceaux,myaux

      x = coord(1)
      y = coord(2)
      z = coord(3)
      
      if (abs(x-a%tx0) < a%xwidth .and. (abs(y-a%ty0) < a%ywidth)) then
         if (a%kfl_isRandom == 1) then
         
            call RANDOM_NUMBER(ft)
            ft = ft-0.5_rp
         else
            fcoordt = 0.0_rp
            do i = 1,nfunctions
               do j = 1,3
                  timeaux = ctime/a%tdt
                  spaceaux = coord(j)/a%zdt
                  myaux = a%FourierSeriesFrequencies(j,i)*(timeaux+spaceaux)
                  
                  fcoordt(j) = fcoordt(j) + (myaux+a%FourierSeriesPhase(i)  )
               enddo
            enddo
            ft = sum(fcoordt)
            ft = cos(ft)
         endif
      
         force = 0.0_rp
         f2 = a%density*a%tampt*ft*exp(-((x-a%tx0)/a%txsc)**2-((y-a%ty0)/a%tysc)**2)
         
         force(2) = force(2) + f2
         ftmean = ftmean + f2
         ftmax = max(ftmax,f2)
         ftmin = min(ftmin,f2)
         nftmean = nftmean+1
         
      endif
      
   end subroutine   
   
end module
