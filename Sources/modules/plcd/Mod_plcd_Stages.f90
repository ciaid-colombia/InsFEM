module Mod_plcd_Stages
   use typre
   implicit none 


   type :: Substage
      real(rp) :: InitialLoadFactor = 0
      real(rp) :: FinalLoadFactor = 1
      real(rp) :: TimeStep = 1
      real(rp) :: TimeInterval = 10
      
      
      integer(ip) :: MaximumNonLinearIterations = 5
      real(rp)    :: IterationsTolerance = 1e-6
      integer(ip) :: NonLinearAlgorithm = 0

      
      real(rp) :: IniTime,EndTime
      real(rp) :: CurrentLoadFactor = 0.0_rp
      real(rp) :: PreviousLoadFactor = 0.0_rp
      real(rp) :: LoadFactorIncrement = 0.0_rp
   end type   

   type :: Stage 
      
      real(rp), allocatable :: NodalForces(:,:)
      integer(ip), allocatable :: kfl_fixno(:,:)
      real(rp), allocatable :: bvess(:,:),funno(:)
      integer(ip), allocatable :: kfl_fixbo(:)     !Element boundary fixity (nboun)
      type(r1p)  , allocatable :: bvnat(:)         !Natural bc values (nboun)
      integer(ip) :: bvnat_coun                    !Counts the dimension of the bvnat array
      integer(ip),  allocatable ::kfl_funbo(:)     !Functions for onboundaries bc (nboun)
      
      integer(ip) :: NumberOfSubstages = 1
      type(Substage), allocatable :: Substages(:)
      
      integer(ip) :: CurrentSubstage = 1
   end type

end module