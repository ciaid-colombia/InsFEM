module Mod_plcd_TD_TopologyOptimization
   use typre
   use Mod_plcd_BaseElmope
   use Mod_plcd_SIMPMaterial
   use Mod_int2str
   use Mod_plcd_SIMP_ChiFilter
   use Mod_plcd_SIMP_PointWise
   use MPI
   use Mod_plcd_TD_PolarizationTensor
   use Mod_plcd_TD_Derivative
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   private
   public SetPointersTD_TopologyOptimization, MemallTD,TurnofTD,AdaptiveTD

   type, extends(PointerSetter) :: SPTD_TopologyOptimization
contains
      procedure :: SpecificSet => SpecificSetTD_TopologyOptimization
   end type

   type(SPTD_TopologyOptimization) :: SetPointersTD_TopologyOptimization

contains

   !To be called from Memall
   subroutine MemallTD(a)
      use Mod_plcd_TD_StochasticTopologyOptimization
      implicit none
      class(PLCDProblem) :: a

      class(SIMPED), pointer :: mySIMPED
      class(PLCDMaterial), pointer :: myMat

      integer(ip) :: imaterial,ielem,nelem,igaus,npoin
      real(rp) :: chi,nu

      class(FiniteElement), pointer :: e


      !First of all we set the p value for the SIMP material
      do imaterial = 1,size(a%Materials)
         myMat => a%Materials(imaterial)%p
         select type(mymat)
         type is (SIMPMaterial)
            !P value is 1 for topological derivative
            call myMat%SetPValue(1.0_rp)
            call myMat%CT%GetPoissonRatio(nu)
            a%TDData%nu = nu
         end select
      end do

      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','MemallSIMP')


      call a%Mesh%GetNelem(nelem)
      chi = 1.0_rp
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

      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','MemallSIMP')

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%palloc(npoin,a%TDData%NodalChi,'NodalChi','SIMPMemall')

      !Initialize topological derivative
      !ElementWise or PointWise SIMP
      allocate(a%TopologicalDerivative)

      call a%TopologicalDerivative%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%TopologicalDerivative%SetMesh(a%Mesh)
      call a%TopologicalDerivative%SetMemor(a%Memor)
      call a%TopologicalDerivative%SetTDData(a%TDData)
      call a%TopologicalDerivative%SetLargeStrains(a%kfl_LargeStrains)
      call a%TopologicalDerivative%SetPTType(mymat%CT)
       if (a%kfl_StochasticTopologyOptimization == 1) then
         call plcd_TDSTO_Initial(a)
      endif

      !In this subroutine we perform a loop through the elements
      !We compute lambda and related data for each element
      !Finally we modify the chi value for each element according to
      !the corresponding correction
      call a%TopologicalDerivative%Initialize
      
     
      
      
      
      
      
   end subroutine

   !To be called from Turnof
   subroutine TurnofTD(a)
      class(PLCDProblem) :: a

      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%pdealloc(npoin,a%TDData%NodalChi,'NodalChi','SIMPMemall')

      !Finalize topological derivative
      call a%TopologicalDerivative%Finalize
      deallocate(a%TopologicalDerivative)

   end subroutine

   !To be called from Adaptive
   subroutine AdaptiveTD(a,itask)
      class(PLCDProblem) :: a
      character(6) :: itask

      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)
      call a%Memor%pdealloc(size(a%TDData%NodalChi),a%TDData%NodalChi,'NodalChi','SIMPMemall')
      call a%Memor%palloc(npoin,a%TDData%NodalChi,'NodalChi','SIMPMemall')

      call a%TopologicalDerivative%Adaptive(a,itask)
   end subroutine

   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SpecificSetTD_TopologyOptimization(d)
      class(SPTD_TopologyOptimization) :: d
      !Element size required by EMDS
      if (a%kfl_TopologyOptimization == 2 .and. a%kfl_InternalLineSearchLoop == 0 .and. a%TDData%kfl_PerformTopologyIteration == 1) then
            call ConcatenateProcedures(ProcHook%Initializations,InitializationsTD)
      endif
   end subroutine

   subroutine InitializationsTD
      procedure(), pointer :: GetStrain => NULL()
      procedure(), pointer :: AddSphericalComponent => NULL()
      integer(ip) :: vsize

      real(rp) :: strain(6)
      real(rp) :: lambda,totalVol,chi,compareValue,StrainEnergy,BK,correctionfactor,meanchi,ierr
      real(rp) :: gtotalvol, gmeanchi, glambda
      real(rp) :: lambdaGaussPoint
      real(rp) :: newObjectiveVolumeFraction,slope,weightfactor = 1.0_rp

      integer(ip) :: idime,jdime, imaterial
      class(SIMPED), pointer :: mySIMPED
      class(PLCDMaterial), pointer :: myMat
      real(rp), pointer :: NodalChi(:) => NULL()
      real(rp) :: Constraint
      
      real(rp) :: suma, gsuma,auxCompliance
      integer(ip) :: npoinLocal
      integer(ip) :: ielem2

      !Do not do this if this is the beginning of the step call
      if (a%itera == 0) return
      
      call a%TopologicalDerivative%Reset

      lambda = 0.0_rp
      totalVol = 0.0_rp
      meanchi = 0.0_rp
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

               call a%TopologicalDerivative%AddGaussPointData(ielem,e,strain,stress,chi)
               a%Stress(ielem)%a(:,igaus) = stress
            enddo

            call a%TopologicalDerivative%PostGaussData(ielem,e)
         end select
      enddo   elements
      
      !call a%FilePostpr%postgp(a%Stress,'StressPre',a%istep,a%ctime,a%Mesh,'voigtT')


      call a%TopologicalDerivative%FinalizeDataComputations
      
      !Now we should add the topological derivative of the various stochastic analysis
      if (a%kfl_StochasticTopologyOptimization == 1) then
         call plcd_TDSTO_AddStochasticTopologicalDerivative(a)
      endif
      
      !call PostprocessTopologicalDerivative(a%TopologicalDerivative,a%FilePostpr,a%itera,a%istep,a%ctime)
      
      call a%TopologicalDerivative%UpdateTopologicalDerivative(a%FilePostpr,a%itera,a%istep,a%ctime)

      call a%Mesh%GetNelem(nelem)
      do ielem = 1,nelem
         !Only for elements of SIMP type
         ElementMatData => a%ElementMaterialsData(ielem)%p
         select type (ElementMatData)
         type is (SIMPED)
            mySIMPED => ElementMatData

            !Load Element
            call a%Mesh%ElementLoad(ielem,e)

            !Default rule of integration, for getting data from Gauss Point
            do igaus = 1,e%pgaus
               e%igaus = igaus
               call a%TopologicalDerivative%GetGaussPointChi(ielem,e,chi)
               call mySIMPED%SetChiValue(e%igaus,chi)
            enddo
         end select
      enddo

      call a%TopologicalDerivative%GetNodalChi(NodalChi)
      
      if (associated(NodalChi)) then
         a%TDData%NodalChi = NodalChi
      end if
      
      !Postprocess Topological Derivative at each iteration
      if (a%kfl_PostprocessDisplacementAtEachIteration == 1) then
         call a%TopologicalDerivative%PostprocessChi(a%FilePostpr,a%itera,a%istep,a%ctime)
      endif
      
      call a%TopologicalDerivative%GetCompliance(auxCompliance)
      if (a%MPirank == a%MPIroot) then 
         write(a%lun_outph,100,advance='no') 'Istep: ', a%istep, ' Compliance: ', auxCompliance
         100 format (A11,I10,A11,E12.5)
         call a%TopologicalDerivative%WriteLogData(a%lun_outph)
         write(a%lun_outph,*) ' '
         call flush(a%lun_outph)
      endif
      
      if (a%TDData%kfl_WhenPerformTopologyOptimization /= 0) a%TDData%kfl_PerformTopologyIteration = 0
      
   end subroutine
end module


