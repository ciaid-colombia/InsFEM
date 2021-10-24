module Mod_TurbulentInletBoundaryConditions
   use typre
   use Mod_Listen
   use Mod_MPIObject 
   use Mod_Octree
   use Mod_Memor
   use Mod_LinkedList
   implicit none
   
   real(rp), parameter :: VonKarmanK = 0.41   !Von Karman Constant
   real(rp), parameter :: Cplus = 5.0 !Rugosity constant
   
   !This implements the methodology described in:
   !Jarrin, Nicolas, et al. "A synthetic-eddy-method for generating inflow conditions for large-eddy simulations." 
   !International Journal of Heat and Fluid Flow 27.4 (2006): 585-593.

   !Variable names are chosen following the paper naming
   
   type, extends(MPIObject) ::  TurbulentInletBoundaryConditions 
      !Data
      real(rp) :: EddyLengthScales(3) = 0.0_rp     !Eddies Length scale
      real(rp) :: R(3,3) = 0.0_rp        !Reynolds stress tensor for generating data
      real(rp) :: U0             !Eddies transport velocity
      real(rp) :: nu = 1.0_rp, delta = 1.0_rp     !For wall law
      integer(ip) :: kfl_InletType = 1   !Default is wall law
      real(rp) :: TurbulenceIntensity = -1.0_rp
      
      
      type(MemoryMan), pointer :: MemorPointer => NULL()
      
      !Variables
      
      
      real(rp) :: lPlaneCoordinates(2,3),gPlaneCoordinates(2,3)     !xmin ymin zmin
      real(rp) :: gGenerationVolumeDimensions(3), gGenerationVolume                                             
      real(rp) :: lGenerationVolumeDimensions(3), lGenerationVolume                                             
                                            
      integer(ip) :: LocalN, GlobalN, sqrtGlobalN           !Number of eddies generated
      real(rp), allocatable :: EddiesCenter(:,:), EddiesSign(:,:)
      real(rp) :: EddyRadius
      
      type(Octree) :: EddiesOctree
      
      real(rp) :: SignalNormalizationCoefficient
      
      real(rp) :: aij_matrix(3,3)
      
      logical :: IsPeriodicY = .false., IsPeriodicZ = .false.
      
      real(rp) :: UTau
      
      !For Prescribed Ys
      integer(ip) :: nPrescribedY
      real(rp)    :: PrescribedY(10)
      
      
contains
      
      procedure :: ReadData
      procedure :: ReadDataMPI
      procedure :: SetPlaneCoordinates
      procedure :: SetGlobalPlaneCoordinates
      procedure :: Initialize
      procedure :: GetVelocity
      procedure :: SetArePeriodicBoundariesYZ
      procedure :: SetViscosity
      procedure :: SetMemor
      procedure :: Finalize
      

   end type
   
   type, extends(RangeGiver) :: TIBCRangeGiver
      
      type(TurbulentInletBoundaryConditions), pointer :: TIBC => NULL()
      
contains
      procedure :: GiveMeRange => GiveRangeTIBC
      procedure :: InitializeTIBC
   
   end type


contains

   subroutine ReadData(a,Listener)
      class(TurbulentInletBoundaryConditions) :: a
      type(ListenFile) :: Listener
      
      integer(ip) :: iPresc
   
      call Listener%Listen('Turbul')
      do while (Listener%words(1) /= 'ENDIN')
         if (Listener%words(1) == 'LENGT') then
            a%EddyLengthScales(1:3) = Listener%param(1:3)
         elseif (Listener%words(1) == 'INTEN') then
            a%TurbulenceIntensity = Listener%param(1)
         elseif (Listener%words(1) == 'REYNO') then   
            a%R(1,1:3) = Listener%param(1:3)
            a%R(2,1:3) = Listener%param(4:6)
            a%R(3,1:3) = Listener%param(7:9)
         elseif (Listener%words(1) == 'U0   ') then
            a%U0 = Listener%param(1)
         elseif (Listener%words(1) == 'TYPE ') then 
            if (Listener%words(2) == 'WALLL') then
               a%kfl_InletType = 1
            elseif (Listener%words(2) == 'WALLC') then
               a%kfl_InletType = 3   
            elseif (Listener%words(2) == 'UNIFO') then
               a%kfl_InletType = 2
            elseif (Listener%words(2) == 'PRESC') then
               a%kfl_InletType = 4
               a%nPrescribedY =  Listener%GetInt('PRESC',0,'NUMBER OF PRESCRIBED Y')
               do iPresc = 1,a%nPrescribedY
                  call Listener%Listen('Turbul')
                  a%PrescribedY(iPresc) = Listener%param(1)
               enddo
            endif
         elseif (Listener%words(1) == 'DELTA') then
            a%Delta = Listener%param(1)
            
            
         endif
   
         call Listener%listen('Turbul')
      enddo
   end subroutine
   
   subroutine ReadDataMPI(a)
      use MPI
      class(TurbulentInletBoundaryConditions) :: a
      
      integer(ip) :: ierr
      
      CALL MPI_BCAST(a%EddyLengthScales,3,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%R,9,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%U0,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%Delta,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%kfl_InletType,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%TurbulenceIntensity,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%nPrescribedY,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
      CALL MPI_BCAST(a%PrescribedY,a%nPrescribedY,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)
      
   end subroutine
   
   subroutine SetMemor(a,Memor)
      class(TurbulentInletBoundaryConditions) :: a
      type(MemoryMan), target :: Memor
      
      a%MemorPointer => Memor
   end subroutine

   subroutine SetPlaneCoordinates(a,PlaneCoordinates)
      class(TurbulentInletBoundaryConditions) :: a
      real(rp) :: PlaneCoordinates(2,3)     !xmin ymin zmin
                                            !xmax ymax zmax
                                            
      a%lPlaneCoordinates = PlaneCoordinates
   end subroutine
   
   subroutine SetGlobalPlaneCoordinates(a,gPlaneCoordinates)
      class(TurbulentInletBoundaryConditions) :: a
      real(rp) :: gPlaneCoordinates(2,3)     !xmin ymin zmin
                                            !xmax ymax zmax
                                            
      a%gPlaneCoordinates = gPlaneCoordinates
   end subroutine
   
   subroutine SetArePeriodicBoundariesYZ(a,IsPeriodicY,IsPeriodicZ)
      class(TurbulentInletBoundaryConditions) :: a
      logical :: IsPeriodicY,IsPeriodicZ
      
      a%IsPeriodicY = IsPeriodicY
      a%IsPeriodicZ = IsPeriodicZ
      
   end subroutine   
   
   subroutine SetViscosity(a,nu)
      class(TurbulentInletBoundaryConditions) :: a
      real(rp) :: nu
      
      a%nu = nu
   end subroutine
      
   
   subroutine Initialize(a)
      class(TurbulentInletBoundaryConditions) :: a
      
      real(rp) :: lGenerationVolumeCoordinates(2,3), gGenerationVolumeCoordinates(2,3), localInletSurface, globalInletSurface, EddySurface
      integer(ip) :: ieddy,idime
      real(rp) :: rand_num
      
      type(TIBCRangeGiver) :: TIBCRG
      
      if (a%lPlaneCoordinates(1,1) > a%lPlaneCoordinates(2,1)) return   !nothing to do for this processor
      
      !Check that x coordinate is the same everywhere
      if (a%gPlaneCoordinates(1,1)-a%gPlaneCoordinates(2,1) > 1e-8) call runend('TurbulentInletBoundaryConditions need to be applied in a x=cte plane')
      
      !Generation Volume dimensions
      lGenerationVolumeCoordinates(1,:) = a%lPlaneCoordinates(1,:) - a%EddyLengthScales(:)
      lGenerationVolumeCoordinates(2,:) = a%lPlaneCoordinates(2,:) + a%EddyLengthScales(:)
      a%lGenerationVolumeDimensions(:) = lGenerationVolumeCoordinates(2,:) - lGenerationVolumeCoordinates(1,:)
      a%lGenerationVolume = a%lGenerationVolumeDimensions(1)*a%lGenerationVolumeDimensions(2)*a%lGenerationVolumeDimensions(3)
      
      !GLOBAL Generation Volume dimensions
      gGenerationVolumeCoordinates(1,:) = a%gPlaneCoordinates(1,:) - a%EddyLengthScales(:)
      gGenerationVolumeCoordinates(2,:) = a%gPlaneCoordinates(2,:) + a%EddyLengthScales(:)
      a%gGenerationVolumeDimensions(:) = gGenerationVolumeCoordinates(2,:) - gGenerationVolumeCoordinates(1,:)
      a%gGenerationVolume = a%gGenerationVolumeDimensions(1)*a%gGenerationVolumeDimensions(2)*a%gGenerationVolumeDimensions(3)
      
      
      
      !Compute the number of eddies
      globalInletSurface = abs(a%gGenerationVolumeDimensions(2)*a%gGenerationVolumeDimensions(3))
      EddySurface = 4*a%EddyLengthScales(2)*a%EddyLengthScales(3)
      a%globalN = ceiling(globalInletSurface/EddySurface)
      a%sqrtGlobalN = sqrt(real(a%GlobalN));
      a%EddyRadius = (a%EddyLengthScales(1)*a%EddyLengthScales(2)*a%EddyLengthScales(3))**(1.0_rp/3.0_rp)
      
      LocalInletSurface = abs(a%lGenerationVolumeDimensions(2)*a%lGenerationVolumeDimensions(3))
      a%localN = ceiling(LocalInletSurface/EddySurface)
      
      !Compute the coordinates for the center of each eddy
      allocate(a%EddiesCenter(3,a%localN))
      allocate(a%EddiesSign(3,a%localN))
      
      !Random generation of eddies
      !Does not need to be communicated because srand(0) generates the same sequence of random numbers everywhere
      call srand(0)
      do ieddy = 1,a%localN
         do idime = 1,3
            call RANDOM_NUMBER(rand_num) 
            a%EddiesCenter(idime,ieddy) = rand_num*a%lGenerationVolumeDimensions(idime)+lGenerationVolumeCoordinates(1,idime)
            
            call RANDOM_NUMBER(rand_num)
            if (rand_num > 0.5_rp) then
               a%EddiesSign(idime,ieddy) = 1.0_rp
            else
               a%EddiesSign(idime,ieddy) = -1.0_rp
            endif
         enddo
      enddo
      
      !Prepare an octree for looking for the eddies efficiently
      !Prepare an Octree of Subdomains
      call TIBCRG%InitializeTIBC(a)
      call a%EddiesOctree%InitializeElementOctree(2,a%localN,TIBCRG,5,a%MemorPointer)
   
      !If we only have turbulence intensity, we build the Reynolds stress tensor from it
      if (a%TurbulenceIntensity > 0.0_rp) then
         a%R = 0.0_rp
         a%R(1,1) = (a%TurbulenceIntensity*a%U0)**2
         a%R(2,2) = (a%TurbulenceIntensity*a%U0)**2
         a%R(3,3) = (a%TurbulenceIntensity*a%U0)**2
      endif
      
      if (maxval(abs(a%R)) == 0.0_rp) then
         a%aij_matrix = 0.0_rp
      else
         !aij_matrix
         a%aij_matrix = 0.0_rp
         a%aij_matrix(1,1) = sqrt(a%R(1,1))
         a%aij_matrix(2,1) = a%R(2,1)/a%aij_matrix(1,1)
         a%aij_matrix(2,2) = sqrt(a%R(2,2)-a%aij_matrix(2,1)**2.0_rp)
         a%aij_matrix(3,1) = a%R(3,1)/a%aij_matrix(1,1)
         a%aij_matrix(3,2) = (a%R(3,2)-a%aij_matrix(2,1)*a%aij_matrix(3,1))/a%aij_matrix(2,2)
         a%aij_matrix(3,3) = sqrt(a%R(3,3)-a%aij_matrix(3,1)**2-a%aij_matrix(3,2)**2)
      endif
      
      
       !Eddies are assumed to be Gaussian distributions,
      !The support goes from p = 0.001, to p = 0.999
      
      !We work with mean = 0, deviation = 1
      !Then we scale it so that the integral on the volume is 1
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !MATLAB CODE FOR THE NORMALIZATION OF THE FUNCTION
!      
!       %Integral on a Sphere of radius 3 (p = 0.999 for the normal distribution
!       %with 0 mean, standard deviation one, of normpdf(r)^2.
! 
!       %Integrating over the sphere means:
! 
!       %Int_sphere(f^2) = 4*pi*Int_r(f^2 r^2))
! 
!       fun = @(r) (r.*normpdf(r)).^2;
!       q = integral(fun,0,3.09);
!       result = q*4*pi
!       !The result for the previous integral is 0.8860   
!
!
!       R=200;
!       Fun = @(r) (normpdf(r/R*3.09)*sqrt((3.09/R)^3/0.8860));
!       Fun2 = @(r) (Fun(r).*Fun(r).*r.*r);
!       q=integral(Fun2,0,R);
!       result2=q*4*pi
!       !The result for the previous integral is 1.0000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !For a given EddyLengthScales, we need to addimensionalize it by multiplying the function
      ! by sqrt((3.09/R)^3/0.8860), and then entering with r/R*3.09 instead of r
      
      !This ensures that the integral is one
      !However, we require that the integral DIVIDED BY THE TOTAL VOLUME is 1, so we multiply it by the total volume
      !Also we will have N eddies, so we divide it by sqrt of N
      
      a%SignalNormalizationCoefficient = (a%gGenerationVolume)*sqrt(((3.09)**3)/(0.8860*(a%EddyRadius**3)))/a%sqrtGlobalN 
      
      
      !For the wall law profile
      call FindUTauFromWallLawParameters(a)
      
   end subroutine
   
   subroutine Finalize(a)
      class(TurbulentInletBoundaryConditions) :: a
      
      if (a%lPlaneCoordinates(1,1) > a%lPlaneCoordinates(2,1)) return   !nothing to do for this processor
      
      deallocate(a%EddiesCenter)
      deallocate(a%EddiesSign)
      
      call a%EddiesOctree%Dealloc(a%Memor)
   end subroutine
   
   subroutine GetVelocity(a,Coord,ctime,Velocity)
      use Mod_Prob
      implicit none
      class(TurbulentInletBoundaryConditions) :: a
      real(rp), intent(in) :: Coord(3),ctime
      real(rp) :: Velocity(3),y
      
      real(rp) :: TurbulentVelocity(3)=0.0_rp, FlowVelocity(3)=0.0_rp
      integer(ip) :: iPresc
      
      call GetTurbulentVelocity(a,Coord,ctime,TurbulentVelocity)
      
      !Wall law
      if (a%kfl_InletType == 1) then
         y = abs(Coord(2)-a%gPlaneCoordinates(1,2))
         call GetWallLawVelocity(a,y,ctime,FlowVelocity)
      !Channel wall law
      elseif (a%kfl_InletType == 3) then
         y = min(abs(Coord(2)-a%gPlaneCoordinates(1,2)),abs(Coord(2)-a%gPlaneCoordinates(2,2)))
         call GetWallLawVelocity(a,y,ctime,FlowVelocity)
      
      !Uniform
      elseif (a%kfl_InletType == 2) then
         FlowVelocity = 0.0_rp
         FlowVelocity(1) = a%U0
      
      !PrescribedY
      elseif (a%kfl_InletType == 4) then   
         y = 1e12
         do iPresc = 1,a%nPrescribedY
            y = min(y,abs(Coord(2) - a%PrescribedY(iPresc)))
         enddo
         call GetWallLawVelocity(a,y,ctime,FlowVelocity)
      endif
      
      
      Velocity = FlowVelocity + TurbulentVelocity
      
   end subroutine    
   

   subroutine GetTurbulentVelocity(a,Coord,ctime,Velocity)
      use Mod_Prob
      implicit none
      class(TurbulentInletBoundaryConditions) :: a
      real(rp), intent(in) :: Coord(3),ctime
      real(rp) :: Velocity(3)
      
      real(rp) :: distvector(3), distance,f, TimeCoord(3), auxTimeCoord(3)
      real(rp) :: DimensionlessDistance
      integer(ip) :: ieddy
      real(rp) :: Signal(3), TimeDistance
      
      Signal = 0.0_rp
      TimeDistance = mod(a%U0*ctime,2*a%EddyLengthScales(1)) 
      
      !Eddies go through the boundary at a velocity U0
      TimeCoord = Coord
      TimeCoord(1) = TimeCoord(1) - a%EddyLengthScales(1) + TimeDistance
      call AddContributionToSignal(TimeCoord)
      
      !We do not do YZ periodicity because of parallelism costs
      !call AddContributionYZPeriodicity(TimeCoord)
      
      !Periodicity of eddies in the x direction
      TimeCoord = Coord
      if (TimeDistance < a%EddyLengthScales(1) ) then
         TimeCoord(1) = TimeCoord(1) + a%EddyLengthScales(1) + TimeDistance
      else
         TimeCoord(1) = TimeCoord(1) -3*a%EddyLengthScales(1) + TimeDistance
      endif
      call AddContributionToSignal(TimeCoord)
      !We do not do YZ periodicity because of parallelism costs
      !call AddContributionYZPeriodicity(TimeCoord)
      
      Velocity = matmul(a%aij_matrix,Signal)
 
 
contains
      subroutine AddContributionToSignal(TimeCoord)
         real(rp) :: TimeCoord(3)
         
         class(LinkedListHeader), pointer  :: Head => NULL()
         class(LinkedListStorage), pointer :: Store => NULL()
         
         !The octree is 2dimensional, only for Y and Z (coord 2 and 3)
         call a%EddiesOctree%GetListForCoord(TimeCoord(2:3),Head,Store)
         call Store%GoToFirstOfList(Head)
         call Store%GetNext(ieddy)
         do while (ieddy /= -1)
            distvector = TimeCoord(:) - a%EddiesCenter(:,ieddy)
            distance = sqrt(dot_product(distvector,distvector))
         
            if (distance < a%EddyRadius) then
               
               DimensionlessDistance = distance/a%EddyRadius*3.09_rp
               
               call normal_01_pdf ( DimensionlessDistance, f)
               f = f*a%SignalNormalizationCoefficient
               
               !We divide by sqrtN outside the loop
               Signal = Signal+a%EddiesSign(:,ieddy)*f
            endif
         
            call Store%GetNext(ieddy)
         enddo
      end subroutine    
      
      subroutine AddContributionYZPeriodicity(TimeCoord)
         real(rp) :: TimeCoord(3)
         
         real(rp) :: auxTimeCoord(3)
         
         if (a%IsPeriodicY) then
            auxTimeCoord = TimeCoord
            if (TimeCoord(2) > a%gPlaneCoordinates(1,2) + (a%gPlaneCoordinates(2,2)-a%gPlaneCoordinates(1,2))/2.0_rp) then
               auxTimeCoord(2) = TimeCoord(2) - (a%gPlaneCoordinates(2,2)-a%gPlaneCoordinates(1,2))
            else
               auxTimeCoord(2) = TimeCoord(2) + (a%gPlaneCoordinates(2,2)-a%gPlaneCoordinates(1,2))
            endif
            call AddContributionToSignal(auxTimeCoord)
         endif   
         
         if (a%IsPeriodicZ) then
            auxTimeCoord = TimeCoord
            if (TimeCoord(3) > a%gPlaneCoordinates(1,3) + (a%gPlaneCoordinates(2,3)-a%gPlaneCoordinates(1,3))/2.0_rp) then
               auxTimeCoord(3) = TimeCoord(3) - (a%gPlaneCoordinates(2,3)-a%gPlaneCoordinates(1,3))
            else
               auxTimeCoord(3) = TimeCoord(3) + (a%gPlaneCoordinates(2,3)-a%gPlaneCoordinates(1,3))
            endif
            call AddContributionToSignal(auxTimeCoord)
         endif
      end subroutine
      
   end subroutine
   
   
   subroutine FindUTauFromWallLawParameters(a)
      class(TurbulentInletBoundaryConditions) :: a
      
      real(rp) :: uplus0,utau0,f0,f1,fderivative, yplus
      integer(ip) :: iiter
      
      !First guess is viscous part
      a%UTau = sqrt(a%U0*a%nu/a%delta)
      
      yplus = a%U0/a%UTau
      if (yplus <= 11) return
      
      
      !If it is not in the viscous part, then log part
      
      !Here we use a secant method in order to find u_tau
      !We assume we are in the logarithmic branch
      
      !Starting guess is y+ = 100
      uplus0 = 1/VonKarmanK*log(100.0_rp)+Cplus 
      utau0 = a%U0/uplus0
      a%UTau = utau0*0.9
      call EvaluateF(utau0,f0)
      call EvaluateF(a%UTau,f1)
      
      do iiter = 1,50
         call EvaluateDerivative(utau0,f0,a%UTau,f1,fderivative)
         
         if (abs(f1) < 1e-6) return
         
         utau0 = a%UTau
         f0 = f1
         
         a%UTau = a%UTau - f1/fderivative
         call EvaluateF(a%UTau,f1)
      enddo

      
      
      
contains      
      subroutine EvaluateF(utau,f)
         real(rp) :: utau,f
         
         f = a%U0/utau - 1/VonKarmanK*log(a%delta*utau/a%nu)-Cplus
      end subroutine   
      
      subroutine EvaluateDerivative(utau0,f0,utau1,f1,fderivative)
         real(rp) :: utau0,f0,utau1,f1,fderivative
         
         fderivative = (f1-f0)/(utau1-utau0)
      end subroutine
   end subroutine
   
   
   
   subroutine GetWallLawVelocity(a,y,ctime,Velocity)
      use Mod_Prob
      implicit none
      class(TurbulentInletBoundaryConditions) :: a
      real(rp), intent(in) :: y,ctime
      real(rp) :: Velocity(3)
      
      real(rp) :: yplus,uplus
      
      Velocity = 0.0_rp
      
      !y+
      yplus = y*a%UTau/a%nu
      
      if (yplus < 11.0_rp) then
         uplus = yplus
      else
         uplus = 1/VonKarmanK*log(yplus)+Cplus
      endif
      
      Velocity(1) = a%UTau*uplus
      
      if (Velocity(1) > a%U0) Velocity(1) = a%U0
   end subroutine
         
   !We need an octree for finding the ranges
    subroutine GiveRangeTIBC(a,ielem,range)
      implicit none
      class(TIBCRangeGiver) :: a
      integer(ip) :: ielem
      real(rp) :: range(*)
      
      integer(ip) :: icount, idime, j
      
      icount = 0
      do idime = 2,3   !only y and z
         icount = icount +1
         range(icount) = a%TIBC%EddiesCenter(idime,ielem)-a%TIBC%EddyLengthScales(idime)
         icount = icount+1
         range(icount) = a%TIBC%EddiesCenter(idime,ielem)+a%TIBC%EddyLengthScales(idime)
      enddo
   end subroutine
   
   subroutine InitializeTIBC(a,TIBC)
      implicit none
      class(TIBCRangeGiver) :: a
      type(TurbulentInletBoundaryConditions), target :: TIBC
      
      a%TIBC => TIBC
   end subroutine
   
   
   
   

end module
