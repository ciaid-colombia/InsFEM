module Mod_tem_TemperatureMaterial
   use typre
   use Mod_MPIObject
   use Mod_Listen
   implicit none
   
   
   type, extends(MPIObjectLight) :: TemperatureMaterial
      !Physical properties
      real(rp) :: &
            densi = 1.0_rp,&                  ! Density (rho)
            sphea = 1.0_rp,&                  ! Specific heat (Cp)
            tcond = 1.0_rp                    ! Thermal conduction (k)
contains
      procedure :: ReadData
      procedure :: ScatterData
      
      procedure :: GetDensity
      procedure :: GetSpecificHeat
      procedure :: GetThermalConductivity
      procedure :: SetDensity
      procedure :: SetSpecificHeat
      procedure :: SetThermalConductivity
      
   end type
contains

   subroutine ReadData(a,Listener)
      implicit none
      class(TemperatureMaterial), target :: a
      type(ListenFile) :: Listener
      
      do while(Listener%words(1)/='ENDMA')
         call Listener%listen('tempe_reaphy')
         if(Listener%words(1) == 'DENSI') then
            a%densi = Listener%param(1)
         elseif(Listener%words(1) == 'SPECI') then
            a%sphea = Listener%param(1)
         elseif(Listener%words(1) == 'THERM') then
            a%tcond = Listener%param(1)
         endif
      enddo
      !This is just so that in the external loop there are no problems
      Listener%words(1) = '     '
   end subroutine
   
   subroutine ScatterData(a)
      use MPI
      implicit none
      class(TemperatureMaterial), target :: a

      integer(ip) :: ierr
      call MPI_BCAST(a%densi,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%sphea,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      call MPI_BCAST(a%tcond,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
   end subroutine
   
   subroutine GetDensity(a,density)
      implicit none
      class(TemperatureMaterial), target :: a
      real(rp) :: density
      
      density = a%densi
   end subroutine
   
   subroutine GetSpecificHeat(a,sphea)
      implicit none
      class(TemperatureMaterial), target :: a
      real(rp) :: sphea
      
      sphea = a%sphea
   end subroutine
   
   subroutine GetThermalConductivity(a,tcond)
      implicit none
      class(TemperatureMaterial), target :: a
      real(rp) :: tcond
      
      tcond = a%tcond
   end subroutine
   
   subroutine SetDensity(a,density)
      implicit none
      class(TemperatureMaterial), target :: a
      real(rp) :: density
      
      a%densi = density
   end subroutine
   
   subroutine SetSpecificHeat(a,sphea)
      implicit none
      class(TemperatureMaterial), target :: a
      real(rp) :: sphea
      
      a%sphea = sphea
   end subroutine
   
   subroutine SetThermalConductivity(a,tcond)
      implicit none
      class(TemperatureMaterial), target :: a
      real(rp) :: tcond
      
      a%tcond = tcond
   end subroutine
   
   



end module