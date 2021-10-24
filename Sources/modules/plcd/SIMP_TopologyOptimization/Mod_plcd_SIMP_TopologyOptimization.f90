module Mod_plcd_SIMP_TopologyOptimization
   use typre
   use Mod_plcd_BaseElmope
   use Mod_plcd_SIMPMaterial
   use Mod_int2str
   use Mod_plcd_SIMP_ChiFilter
   use Mod_plcd_SIMP_PointWise
   use MPI
   implicit none
   private
   public SetPointersSIMP_TopologyOptimization, MemallSIMP,TurnofSIMP,AdaptiveSIMP

   type, extends(PointerSetter) :: SPSIMP_TopologyOptimization
contains
      procedure :: SpecificSet => SpecificSetSIMP_TopologyOptimization
   end type

   type(SPSIMP_TopologyOptimization) :: SetPointersSIMP_TopologyOptimization

contains

   !To be called from Memall
   subroutine MemallSIMP(a)
      class(PLCDProblem) :: a

      class(SIMPED), pointer :: mySIMPED
      class(PLCDMaterial), pointer :: myMat

      integer(ip) :: imaterial,ielem,nelem,igaus,npoin
      real(rp) :: chi

      class(FiniteElement), pointer :: e


      !First of all we set the p value for the SIMP material
      do imaterial = 1,size(a%Materials)
         myMat => a%Materials(imaterial)%p
         select type(mymat)
         type is (SIMPMaterial)
            call myMat%SetPValue(a%SIMPData%p)
         end select
      end do

      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','MemallSIMP')

      call a%Mesh%GetNelem(nelem)
      chi = a%SIMPData%VolumeFraction
      do ielem = 1,nelem
         !Only for elements of SIMP type
         ElementMatData => a%ElementMaterialsData(ielem)%p
         select type (ElementMatData)
         type is (SIMPED)
            mySIMPED => ElementMatData
            call a%Mesh%ElementLoad(ielem,e)

            do igaus = 1,e%pgaus
               e%igaus = igaus
               call mySIMPED%SetChiValue(e%igaus,chi)
            end do
         end select
      end do
      a%SIMPData%HistoryObjectiveVolumeFraction(1:3) = a%SIMPData%VolumeFraction
      a%SIMPData%HistoryMeanChi(1:3) = a%SIMPData%VolumeFraction

      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','MemallSIMP')

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%palloc(npoin,a%SIMPData%NodalChi,'NodalChi','SIMPMemall')

   end subroutine

   !To be called from Turnof
   subroutine TurnofSIMP(a)
      class(PLCDProblem) :: a

      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%pdealloc(npoin,a%SIMPData%NodalChi,'NodalChi','SIMPMemall')
   end subroutine

   !To be called from Adaptive
   subroutine AdaptiveSIMP(a)
      class(PLCDProblem) :: a

      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%pdealloc(size(a%SIMPData%NodalChi),a%SIMPData%NodalChi,'NodalChi','SIMPMemall')
      call a%Memor%palloc(npoin,a%SIMPData%NodalChi,'NodalChi','SIMPMemall')
   end subroutine

   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SpecificSetSIMP_TopologyOptimization(d)
      class(SPSIMP_TopologyOptimization) :: d
      !Element size required by EMDS
      if (a%kfl_TopologyOptimization == 1) then
         call ConcatenateProcedures(ProcHook%Initializations,InitializationsSIMP)
      endif
   end subroutine

   subroutine InitializationsSIMP

      procedure(), pointer :: GetStrain => NULL()
      procedure(), pointer :: AddSphericalComponent => NULL()
      integer(ip) :: vsize

      real(rp) :: strain(6)
      real(rp) :: lambda,totalVol,chi,newchi,compareValue,StrainEnergy,BK,correctionfactor,meanchi,ierr
      real(rp) :: gtotalvol, gmeanchi, glambda
      real(rp) :: lambdaGaussPoint
      real(rp) :: newObjectiveVolumeFraction,slope,weightfactor = 1.0_rp


      integer(ip) :: idime,jdime, imaterial

      class(SIMPED), pointer :: mySIMPED
      class(PLCDMaterial), pointer :: myMat

      class(ChiFilterInterface), pointer :: ChiFilter
      real(rp), pointer :: NodalChi(:)


      !Do not do this if this is the beginning of the step call
      if (a%itera == 0) return

      !ElementWise or PointWise SIMP
      if (a%SIMPData%ChispaceChoice == 0) then
         allocate(ElementCenterChiFilter::ChiFilter)
      elseif (a%SIMPData%ChispaceChoice == 1) then
         allocate(LumpedMassMatrixChiFilter::ChiFilter)
      endif

      !In this subroutine we perform a loop through the elements
      !We compute lambda and related data for each element
      !Finally we modify the chi value for each element according to
      !the corresponding correction

      call ChiFilter%Initialize(a)

      lambda = 0.0_rp
      totalVol = 0.0_rp
      !Element loop to compute lambda
      call a%Mesh%GetNelem(nelem)
      elements : do ielem = 1,nelem
         !Only for elements of SIMP type
         ElementMatData => a%ElementMaterialsData(ielem)%p
         select type (ElementMatData)
         type is (SIMPED)
            mySIMPED => ElementMatData

            !Load Element
            call a%Mesh%ElementLoad(ielem,e)
            call e%gather(e%ndime,eldisp,a%Displacement(:,:,1))
            call e%GetWeightFactor(weightfactor)

            !Default rule of integration, for getting data from Gauss Point
            do igaus = 1,e%pgaus
               e%igaus = igaus
               !Compute derivatives
               call e%elmder
               dvol = e%weigp(e%igaus)*e%detjm

               !ComputeDisplacementGradients
               call e%gradient(e%ndime,eldisp(:,:,1),gradDisp)

               call mySIMPED%ComputeHistoryAndConstitutiveTensor(e%igaus,gradDisp)
               call mySIMPED%ComputeStress(e%igaus,gradDisp)
               call mySIMPED%Material%CT%GetStrain(gradDisp,strain)
               call mySIMPED%GetStressTensorPointer(e%igaus,stress)

               !We assume only one value elementwise at this point
               call mySIMPED%GetChivalue(e%igaus,chi)

               !Original strain energy
               StrainEnergy = dot_product(stress,strain(1:size(stress)))/chi*a%SIMPData%p

               call ChiFilter%AddGaussPointContributionToStrainEnergy(ielem,e,StrainEnergy)
               call ChiFilter%AddGaussPointContributionToChi(ielem,e,chi)

               lambda = lambda + StrainEnergy*dvol*weightfactor
               totalVol = totalVol + dvol*weightfactor
            enddo

            call ChiFilter%PostGaussChi(ielem,e)
            call ChiFilter%PostGaussStrainEnergy(ielem,e)

         end select
      enddo   elements

      call ChiFilter%FinalizeStrainEnergyComputations
      call ChiFilter%FinalizeChiComputations

      glambda = 0.0_rp
      gtotalvol = 0.0_rp
      call MPI_REDUCE( lambda, glambda, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( totalvol, gtotalvol, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )

      lambda = glambda/gtotalvol

      CALL MPI_BCAST(lambda,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)

      if (lambda == 0.0_rp) then
         call ChiFilter%Finalize
         return
      endif



      call ChiFilter%UpdateChi(lambda)

      !Now we update the density function

      meanchi = 0.0_rp
      totalVol = 0.0_rp

      do ielem = 1,nelem
         !Only for elements of SIMP type
         ElementMatData => a%ElementMaterialsData(ielem)%p
         select type (ElementMatData)
         type is (SIMPED)
            mySIMPED => ElementMatData

            call a%Mesh%ElementLoad(ielem,e)
            call e%GetWeightFactor(weightfactor)

            do igaus = 1,e%pgaus
               e%igaus = igaus

               !Compute derivatives
               call e%elmder
               dvol = e%weigp(e%igaus)*e%detjm

               call ChiFilter%GetGaussPointChi(ielem,e,newchi)

               call mySIMPED%SetChiValue(e%igaus,newchi)

               meanchi = meanchi + newchi*dvol*weightfactor
            enddo
         end select
      enddo

      gmeanchi = 0.0_rp

      call MPI_REDUCE( meanchi, gmeanchi, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )

      meanchi = gmeanchi/gtotalvol

      CALL MPI_BCAST(meanchi,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)

      call ChiFilter%ResetChi

      !Correcting for the objective mass value
      correctionfactor = a%SIMPData%HistoryObjectiveVolumeFraction(1)/meanchi

      meanchi = 0.0_rp
      totalVol = 0.0_rp
      do ielem = 1,nelem
         !Only for elements of SIMP type
         ElementMatData => a%ElementMaterialsData(ielem)%p
         select type (ElementMatData)
         type is (SIMPED)
            mySIMPED => ElementMatData

            call a%Mesh%ElementLoad(ielem,e)
            call e%GetWeightFactor(weightfactor)

            do igaus = 1,e%pgaus
               e%igaus = igaus

               call mySIMPED%GetChivalue(e%igaus,chi)
               chi = max(min(chi * correctionfactor,a%SIMPData%chimax),a%SIMPData%chimin)

               call mySIMPED%SetChiValue(e%igaus,chi)

               call ChiFilter%AddGaussPointContributionToChi(ielem,e,chi)

               meanchi = meanchi + chi*dvol*weightfactor
               totalVol = totalVol + dvol*weightfactor
            enddo
            call ChiFilter%PostGaussChi(ielem,e)

         end select
      enddo

      call ChiFilter%FinalizeChiComputations
      call ChiFilter%PostprocessChi(a%FilePostpr)

      call ChiFilter%GetNodalChi(NodalChi)

      if (associated(NodalChi)) then
         a%SIMPData%NodalChi = NodalChi
      end if

      call ChiFilter%Finalize
      !-----------------------------------------------------------------------------------------------------------------
      !method for deciding which is the volume fraction value
      !which will give us the seeked mean chi
      !Correction based on the previous iterations
      gmeanchi = 0.0_rp
      call MPI_REDUCE( meanchi, gmeanchi, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      meanchi = gmeanchi/gtotalvol
      CALL MPI_BCAST(meanchi,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)

      a%SIMPData%HistoryMeanChi(1) = meanchi
      !I use a linear relationship
      newObjectiveVolumeFraction = a%SIMPData%VolumeFraction/meanchi*a%SIMPData%HistoryObjectiveVolumeFraction(1)

      a%SIMPData%HistoryMeanChi(2:3) = a%SIMPData%HistoryMeanChi(1:2)
      a%SIMPData%HistoryObjectiveVolumeFraction(2:3) = a%SIMPData%HistoryObjectiveVolumeFraction(1:2)
      a%SIMPData%HistoryObjectiveVolumeFraction(1) = newObjectiveVolumeFraction
      !-----------------------------------------------------------------------------------------------------------------
   end subroutine
end module


