subroutine plcd_ComputeForcesVectors(a)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem) :: a

   interface
      subroutine plcd_EnditeElmope(a)
         use typre
         use Mod_PLCD
         implicit none
         class(PLCDProblem) :: a
      end subroutine
      
      subroutine plcd_Enditebouope(a)
         use typre
         use MOD_PLCD
         implicit none
         class(PLCDProblem) :: a
      end subroutine
      
      subroutine plcd_bouope(a)
         use typre
         use Mod_PLCD
         implicit none
         class(PLCDProblem) :: a
      end subroutine
   end interface

   integer(ip) :: ndime

   !Set Internal and External Forces to Zero to check convergence
   a%InternalForcesVector = 0.0_rp
   a%ExternalForcesVector = 0.0_rp
   call a%Mesh%GetNdime(ndime)
   a%ExternalForcesVector(1:ndime,:) = a%css%CurrentLoadFactor*a%cs%NodalForces
      
   !Add Tractions in External Forces Vector
   call plcd_Enditebouope(a)

   !So that we know that we are not inside the endite, but in the begste
   call plcd_EnditeElmope(a)

   call a%Mesh%ArrayCommunicator%GhostCommunicate(size(a%InternalForcesVector,1),a%InternalForcesVector)
   call a%Mesh%ArrayCommunicator%GhostCommunicate(size(a%ExternalForcesVector,1),a%ExternalForcesVector)

   !Computation of the residual vector
   a%ResidualForcesVector = a%ExternalForcesVector-a%InternalForcesVector
end subroutine
