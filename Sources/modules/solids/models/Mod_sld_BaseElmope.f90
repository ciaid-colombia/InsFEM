module Mod_sld_BaseElmope 
  use Mod_Mesh
  use Mod_php_SetTimeIntegrator
  use Mod_CauchyElement
  use Mod_sld_elmdir
  use Mod_Solids
  use Mod_UPSolids
  use Mod_UPSolids_NH
  use Mod_SUPSolids
  use Mod_SUPSolids_NH
  use Mod_SUPSolids_lin
  use Mod_SUPSolids_lin
  use Mod_sldup_calculateAU

  implicit none 
 
  real(rp), parameter :: zeroc = epsilon(0.0_rp)
 
  class(SolidsProblem),        pointer :: a       => NULL()
  class(SUPSolidsProblem),     pointer :: sup     => NULL()
  class(SUPSolidsProblem_lin), pointer :: sup_lin => NULL()
  class(SUPSolidsProblem_NH) , pointer :: sup_NH  => NULL()
  class(UPSolidsProblem)     , pointer :: up      => NULL()
  class(UPSolidsProblem_NH)  , pointer :: up_NH   => NULL()

  
  character(6) :: itask

  type :: PPointer  
      !Pointers
      procedure(),                NOPASS, pointer :: ExternalForces   => NULL()
      procedure(),                NOPASS, pointer :: PostGauss        => NULL()
      procedure(sld_elmk),        NOPASS, pointer :: sld_elmk         => NULL()
      procedure(sld_elmgeo),      NOPASS, pointer :: sld_elmgeo       => NULL()
      procedure(sld_elmrhu),      NOPASS, pointer :: sld_elmrhu       => NULL()
      procedure(getTauParameters),NOPASS, pointer :: getTauParameters => NULL()

  end type
  type(PPointer) :: ProcPointer

  !Hooks

  type :: PHook
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange   => NULL()
      procedure(), NOPASS, pointer :: ElmatsToZero    => NULL()
      procedure(), NOPASS, pointer :: Gathers         => NULL()
      procedure(), NOPASS, pointer :: PreGauss        => NULL()
      procedure(), NOPASS, pointer :: InGauss         => NULL()
      procedure(), NOPASS, pointer :: EndGauss        => NULL()
      procedure(), NOPASS, pointer :: PostGauss       => NULL()
      procedure(), NOPASS, pointer :: Dynamic         => NULL()
      procedure(), NOPASS, pointer :: DynamicForce    => NULL()
      procedure(), NOPASS, pointer :: DynamicMass     => NULL()
      procedure(), NOPASS, pointer :: calculateResidual    => NULL()
      procedure(), NOPASS, pointer :: Interpolates         => NULL()
      procedure(), NOPASS, pointer :: PrePostInterpolates  => NULL()
      procedure(), NOPASS, pointer :: PostInterpolates     => NULL()
      procedure(), NOPASS, pointer :: PostPostInterpolates => NULL()
      procedure(), NOPASS, pointer :: InGaussVec           => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats        => NULL()
      procedure(), NOPASS, pointer :: PreDirichlet         => NULL()
      procedure(), NOPASS, pointer :: Finalizations        => NULL()
      procedure(), NOPASS, pointer :: ResetArrays          => NULL()
      procedure(), NOPASS, pointer :: AllocateArrays       => NULL()
      procedure(), NOPASS, pointer :: DeallocateArrays     => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp         => NULL()        
      procedure(), NOPASS, pointer :: preAssemblyLhs       => NULL()
      procedure(), NOPASS, pointer :: AssemblyLhs          => NULL()
      procedure(), NOPASS, pointer :: postAssemblyLhs      => NULL()
      procedure(), NOPASS, pointer :: preAssemblyRhs       => NULL()
      procedure(), NOPASS, pointer :: AssemblyRhs          => NULL()
      procedure(), NOPASS, pointer :: postAssemblyRhs      => NULL()
      procedure(), NOPASS, pointer :: postAssembly         => NULL()
      procedure(), NOPASS, pointer :: ToLinearSystem       => NULL()
      procedure(), NOPASS, pointer :: ComputeResidual      => NULL()

  end type

  class(FiniteElement) , pointer :: e => NULL()
  type(PHook)                    :: ProcHook 
  type(TimeIntegratorDt2)        :: Integrator,IntegratorSGS
  integer(ip)                    :: ielty0 = 0 !Previous element type
  integer(ip)                    :: ielem = 0,nelem = 0,nsteps = 0,igaus = 0,nd = 0
  integer(ip)                    :: kfl_GoIteInGauss = 0,kfl_setPointForce = 0
  real(rp)                       :: det = 1.0_rp,densi = 0.0_rp,dvol = 0.0_rp
  real(rp)                       :: dvolt0 = 0.0_rp,divdisp
  real(rp)                       :: LHSdtinv2,LHSdtinv2SGS,tsch_deriv,gppress(1)
  real(rp), allocatable          :: elmat(:,:,:,:),elrhu(:,:),elrhd(:,:),elrhp(:,:),elrhs(:,:)
  real(rp), allocatable          :: c_elas(:,:),force(:,:)
  real(rp), allocatable          :: eldisp(:,:,:),elveloc(:,:,:), elaccel(:,:,:),elsou(:,:)
  real(rp), allocatable          :: elext(:) 
  real(rp), allocatable          :: gpdisp(:,:),gpsou(:)
  real(rp), allocatable          :: gradDisp(:,:),F_mat(:,:),F_matt(:,:),F_spat(:,:)
  real(rp), allocatable          :: piolaK2(:,:),tensor(:,:),grdFmat(:,:),grdFspt(:,:)
  real(rp), allocatable          :: ii(:,:) !identity matrix
  real(rp), allocatable          :: invar(:) !Tensor invariants
  real(rp), allocatable          :: B(:,:),gpaccel(:)
  real(rp), allocatable          :: strain(:)
  real(rp), allocatable          :: strain_aux(:)
  real(rp), allocatable          :: stress(:)
  real(rp), allocatable          :: elvmass(:),vmass(:)
  !---------------------Assembly matrices------------------------
  real(rp), allocatable          :: massmat(:,:,:,:)
  real(rp), allocatable          :: kmat(:,:,:,:)
  real(rp), allocatable          :: geomat(:,:,:,:)
  real(rp), allocatable          :: intForce(:,:)
  !--------FOR SLDSUP CALCULATIONS----------
  real(rp)                       :: K = 0.0_rp,G = 0.0_rp,lam=0.0_rp,Dv_scalar = 0.0_rp
  real(rp)                       :: tautau=0.0_rp
  real(rp), allocatable          :: AU(:),derJFki(:),divFmat(:),divFspt(:)
  real(rp), allocatable          :: C_dev(:,:),D_dev(:,:),D(:,:),P(:,:)
  real(rp), allocatable          :: elsigma(:,:),elpress(:,:) 
  real(rp), allocatable          :: gpsigma(:),grsig(:,:),divstr(:),grdpre(:)
  real(rp), allocatable          :: elmuv(:,:,:,:),elmut(:,:,:,:),elmuq(:,:,:,:) 
  real(rp), allocatable          :: elmsv(:,:,:,:),elmst(:,:,:,:),elmsq(:,:,:,:) 
  real(rp), allocatable          :: elmpv(:,:,:,:),elmpt(:,:,:,:),elmpq(:,:,:,:) 
  !--------SGS-----
  !--------Tau_u dynamic SGS stabilization-----
  real(rp), allocatable          :: elm_usgsdyn_VU(:,:,:,:)
  real(rp), allocatable          :: elm_usgsdyn_VP(:,:,:,:)
  !--------Tau_u dynamic SGS stabilization-----
  real(rp), allocatable          :: elm_usgs_VU(:,:,:,:)
  real(rp), allocatable          :: elm_usgs_VS(:,:,:,:)
  real(rp), allocatable          :: elm_usgs_VP(:,:,:,:)
  !--------Tau_u dynamic stabilization-----
  real(rp), allocatable          :: elm_usgs_EU(:,:,:,:)
  real(rp), allocatable          :: elm_usgs_QU(:,:,:,:)
  !--------Tau_u stabilization-----
  real(rp), allocatable          :: elm_usgs_ES(:,:,:,:),elm_usgs_EP(:,:,:,:)
  real(rp), allocatable          :: elm_usgs_QS(:,:,:,:),elm_usgs_QP(:,:,:,:)
  !--------Tau_s stabilization-----
  real(rp), allocatable          :: elm_ssgs_VU(:,:,:,:),elm_ssgs_VS(:,:,:,:)
  real(rp), allocatable          :: elm_ssgs_EU(:,:,:,:),elm_ssgs_ES(:,:,:,:)
  !--------Tau_p stabilization-----
  real(rp), allocatable          :: elm_psgs_VU(:,:,:,:),elm_psgs_VP(:,:,:,:)
  real(rp), allocatable          :: elm_psgs_QU(:,:,:,:),elm_psgs_QP(:,:,:,:)

  !Iterfaces for sldsup base procedures
  interface
      module subroutine getTauParametersUP(tau_u,tau_s,tau_p)
          real(rp), intent(inout):: tau_u,tau_s,tau_p
      end subroutine
      module subroutine getTauParameters(tau_u,tau_s,tau_p)
          real(rp), intent(inout):: tau_u,tau_s,tau_p
      end subroutine
      module subroutine getTauParameters_linear(tau_u,tau_s,tau_p)
          real(rp), intent(inout):: tau_u,tau_s,tau_p
      end subroutine

      module subroutine dynamicResidualUP
      end subroutine
      module subroutine calculateAU_U_dyn
      end subroutine
      module subroutine calculateAU_U
      end subroutine
      module subroutine calculateAU_P
      end subroutine
      module subroutine calculateGradients
      end subroutine
      module subroutine calculateFmatGradient
      end subroutine

      !----------------SUP------------
      module subroutine sigmaGather
      end subroutine 
      module subroutine pressGather
      end subroutine  
      module subroutine InterpolateGpSigma
      end subroutine
      module subroutine InterpolateGpPress
      end subroutine
      module subroutine ResetSUPBaseArrays 
      end subroutine
      module subroutine ResetSUPBaseArrays_SGS_dynsgs
      end subroutine
      module subroutine ResetSUPBaseArrays_SGS_dyn 
      end subroutine
      module subroutine calculatePressGradient
      end subroutine
      module subroutine calculatePressGradientsAndStressDiv
      end subroutine
      module subroutine AssembleRhsSUP
      end subroutine
      module subroutine AssembleRhsForcesSUP
      end subroutine
      module subroutine AssembleLhsSUP
      end subroutine
      module subroutine AssembleStiffnessMatrixSUP
      end subroutine
      module subroutine AssembleStiffnessMatrixSUP_SGS
      end subroutine
      module subroutine AssembleStiffnessMatrixSUP_SGS_dynsgs
      end subroutine
      module subroutine AssembleStiffnessMatrixSUP_SGS_dyn
      end subroutine
      module subroutine AllocateSUPAssemblyMatrices
      end subroutine
      module subroutine DeallocateSUPAssemblyMatrices
      end subroutine
      module subroutine AllocateBaseElmopeMatricesSUP
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesSUP
      end subroutine
      module subroutine AllocateBaseElmopeMatricesSUPLinear
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesSUPLinear
      end subroutine
      module subroutine ResetSUPNonlinearDerTerms
      end subroutine
      module subroutine ResetSUPBaseArrays_SGS
      end subroutine
      module subroutine AllocateHigherOrderDerivativeArrays
      end subroutine
      module subroutine DeallocateHigherOrderDerivativeArrays
      end subroutine
      module subroutine AllocateBaseElmopeMatricesSUPNonLinear
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesSUPNonLinear 
      end subroutine
      module subroutine AllocateBaseElmopeMatricesSUP_SGS
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesSUP_SGS
      end subroutine
      module subroutine AllocateBaseElmopeMatricesSUP_SGS_dynsgs
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesSUP_SGS_dynsgs
      end subroutine
      module subroutine AllocateBaseElmopeMatricesSUP_SGS_dyn
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesSUP_SGS_dyn
      end subroutine
      !-------------------UP-----------------
      module subroutine AssembleMassMatrixUP
      end subroutine 
      module subroutine calculateGPPress
      end subroutine
      module subroutine calculateGPStress
      end subroutine
      module subroutine calculateStressDivergence2D
      end subroutine
      module subroutine calculateStressDivergence3D
      end subroutine
      module subroutine calculateSigmaPress
      end subroutine
      module subroutine calculateUPGradients
      end subroutine
      module subroutine AssembleStiffnessMatrixUP
      end subroutine
      module subroutine AssembleRhsUP
      end subroutine
      module subroutine AssembleRhsForcesUP
      end subroutine
      module subroutine AssembleStiffnessMatrixUP_dynsgs
      end subroutine
      module subroutine AllocateBaseElmopeMatricesUP_dynsgs
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesUP_dynsgs
      end subroutine
      module subroutine ResetUPBaseArrays_dynsgs
      end subroutine
      module subroutine AssembleStiffnessMatrixUP_SGS
      end subroutine
      module subroutine ResetUPBaseArrays
      end subroutine
      module subroutine ResetUPBaseArrays_SGS
      end subroutine
      module subroutine AllocateBaseElmopeMatricesUP_SGS
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesUP_SGS
      end subroutine
      module subroutine AllocateUPAssemblyMatrices
      end subroutine
      module subroutine DeallocateUPAssemblyMatrices
      end subroutine
      module subroutine AllocateBaseElmopeMatricesUP
      end subroutine
      module subroutine DeallocateBaseElmopeMatricesUP
      end subroutine

  end interface

  
!#$COMPOSEPROCS 200
#include "COMPOSEPROCS_POINTERS_200.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_200.i90"

  subroutine SetPointersAndHooksToNULLSUB
  implicit none
     !External Procedures
     procedure() :: NULLSUB
  
     !Pointers
     ProcPointer%ExternalForces   => NULLSUB
     ProcPointer%PostGauss        => NULLSUB
     ProcPointer%getTauParameters => NULLSUB
  
     !Hooks
     ProcHook%Initializations      => NULLSUB
     ProcHook%OnIeltyChange        => NULLSUB
     ProcHook%ElmatsToZero         => NULLSUB
     ProcHook%Gathers              => NULLSUB
     ProcHook%PreGauss             => NULLSUB
     ProcHook%InGauss              => NULLSUB
     ProcHook%EndGauss             => NULLSUB
     ProcHook%PostGauss            => NULLSUB
     ProcHook%Dynamic              => NULLSUB
     ProcHook%DynamicForce         => NULLSUB
     ProcHook%DynamicMass          => NULLSUB
     ProcHook%calculateResidual    => NULLSUB
     ProcHook%Interpolates         => NULLSUB
     ProcHook%PrePostInterpolates  => NULLSUB
     ProcHook%PostInterpolates     => NULLSUB
     ProcHook%PostPostInterpolates => NULLSUB
     ProcHook%InGaussElmats        => NULLSUB
     ProcHook%PreDirichlet         => NULLSUB
     ProcHook%InGaussVec           => NULLSUB
     ProcHook%Finalizations        => NULLSUB
     ProcHook%ResetArrays          => NULLSUB
     ProcHook%AllocateArrays       => NULLSUB
     ProcHook%DeallocateArrays     => NULLSUB
     ProcHook%PhysicalProp         => NULLSUB
     ProcHook%preAssemblyLhs       => NULLSUB
     ProcHook%AssemblyLhs          => NULLSUB
     ProcHook%postAssemblyLhs      => NULLSUB
     ProcHook%preAssemblyRhs       => NULLSUB
     ProcHook%AssemblyRhs          => NULLSUB
     ProcHook%postAssemblyRhs      => NULLSUB
     ProcHook%postAssembly         => NULLSUB
     ProcHook%ToLinearSystem       => NULLSUB
     ProcHook%ComputeResidual      => NULLSUB
  end subroutine

  subroutine calculateGradientsAndDeter
     implicit none
     integer(ip) :: uf

     uf       = a%udofn
     gradDisp = 0.0_rp

     call e%gradient(uf,eldisp(:,:,1),gradDisp)

     call calculateMaterialTensor

  end subroutine

  subroutine calculateBoundaryGradientsAndDeter
     implicit none
     integer(ip) :: uf

     uf       = a%udofn
     gradDisp = 0.0_rp

     call e%gradientb(uf,eldisp(:,:,1),gradDisp)

     call calculateMaterialTensor

  end subroutine

  subroutine calculateMaterialTensor
    implicit none
    integer(ip) :: idime,jdime,nd

    call a%Mesh%GetNdime(nd)

    if(a%sld_type== 'NONLI' ) then

        det    = 0.0_rp
        F_mat  = 0.0_rp
        F_spat = 0.0_rp
        call get2ndIITensor(nd,ii)

        !we transform gradDisp into deformation tensor F
        !u = x_sp - X_mat; grad(u) = ii - F(x)^-1
        F_spat= ii - gradDisp 

        !F(X) = F^-1(x)
        call invmtx(F_spat,F_mat,det,nd)

        !calculate Jacobian of F
        call deter(F_mat,det,nd)

        !Tranpose of deformation tensor
        F_matt=transpose(F_mat)

        !calculate left cauchy green deformation tensor
        B = 0.0_rp
        B = matmul(F_mat,F_matt)
        !B = 2.0_rp*getSymGrad(nd,gradDisp)

        if(det .lt. zesld) then
            a%kfl_foldedElements = .true.
            det = 1.0_rp
            !call runend('Solid Base Elmope: Negative element determinant')
        end if
    else 

        gradDisp = getSymGrad(nd,gradDisp)

    end if

  end subroutine

  
  subroutine calculateStress(nd,elem)
    implicit none
    integer(ip),optional   :: elem
    integer(ip),intent(in) :: nd
    integer(ip)            :: sz
    real (rp)              :: traceaux,auxvar
    real(rp)               :: pushForwardaux(nd,nd)
    
    !Remember to call calculateGradientsAndDeter everytime before 
    !using this function!!!

    sz  = (nd*(nd+1))/2
    stress = 0.0_rp
    strain = 0.0_rp

    select case(a%sld_type)

    case('LINEA')

        call getVoigtStrain(sz,nd,gradDisp,strain)

        call a%GetPhysicalParameters(1.0_rp,sz,densi,c_elas,elem)

        stress = matmul(c_elas,strain)

    case('NONLI')

        tensor= 0.0_rp
        piolaK2 = 0.0_rp
        call get2ndIITensor(nd,ii)

        if(a%sld_model =='NEOHO') then

            !This model is given in spatial configuration naturally 

            !calculate cauchy stress tensor
            tensor = (1.0_rp/det)*(a%lambda*log(det)*ii + a%mu*(B-ii))

        else if(a%sld_model =='STVEN') then

            !This model is given in total lagrangian formulation
            !So first we have to calculate the PK2 and then pass it
            !To spatial configuration for updated lagrangian

            !we compute Green strain tensor
            tensor = 0.5*(matmul(F_matt,F_mat) - ii)
            call getVoigtStrain(sz,nd,tensor,strain)

            !calculate PK2
            call trace(nd,tensor,traceaux)
            piolaK2 = a%lambda*traceaux*ii + a%mu*2.0_rp*(tensor)

            !Pushforward of PK2
            call pushForward(nd,sz,F_mat,piolaK2,pushForwardaux)

            !calculate kirchoff stress and solve for cauchy => tau = J sigma
            tensor = (1.0_rp/det)*pushForwardaux

        else

            call runend('Solid Base Elmope: wrong type of material model')

        end if

        !Get stress tensor into voigt notation
        call getVoigtStress(sz,nd,tensor,stress)

        !Now we calculate strains for post process
        select case(a%sld_strain)

        case('GREEN')

            !-----Compute Green strain tensor
            tensor= 0.5*(matmul(F_matt,F_mat) - ii)

        case('ALMAN')
            !-----Compute almansi strain tensor
            call invmtx(B,B,det,nd)
            tensor= 0.5*(ii - B)

        case('SMALL')
            !-----Compute small strain tensor
            tensor = getSymGrad(nd,gradDisp)

        case default !default is almansi strains

            call invmtx(B,B,det,nd)
            tensor= 0.5*(ii - B)

        end select  
        !get into voigt notation
        call getVoigtStrain(sz,nd,tensor,strain)

    end select  

end subroutine

  function getSymGrad(nd,gradDisp) result(sym)
     implicit none
     integer(ip) :: idime,jdime
     integer(ip), intent(in)    :: nd
     real(rp), intent(in) :: gradDisp(nd,nd)
     real(rp) :: sym(nd,nd)

     sym=0.0_rp

     !symmetric gradient
     do idime=1,nd
         do jdime=1,nd
             sym(idime,jdime) = sym(idime,jdime) &
                     + (gradDisp(idime,jdime) + gradDisp(jdime,idime))
         end do
     end do

     sym=0.5_rp*sym

  end function

  function getPrincipalStrains(nd,sz,str) result(res)
    implicit none
    integer(ip), intent(in) :: nd,sz
    real(rp)   , intent(in) :: str(sz)
    real(rp)                :: val_1,val_2,PI,Q,R,theta
    real(rp)                :: res(nd)

    if (nd.eq.2) then

        val_1=0.5_rp*(str(1)+str(2)) + sqrt((0.5_rp*(str(1)-str(2)))**2 + (0.5_rp*str(3))**2)
        val_2=0.5_rp*(str(1)+str(2)) - sqrt((0.5_rp*(str(1)-str(2)))**2 + (0.5_rp*str(3))**2)


         if(val_1 >= val_2) then
             res(1)=val_1
             res(2)=val_2
         else 
             res(1)=val_2
             res(2)=val_1
         end if
        

    else 

        !-----calculate invariants----
        invar=0.0_rp
        call getTensorInvariants(nd,gradDisp,invar)

        !Solution for a cubic polynomial taken from wolfram
        !http://mathworld.wolfram.com/CubicFormula.html

        Q=(3.0_rp*invar(2) - invar(1)*invar(1))/9.0_rp

        if (Q .lt. 0.0_rp) then

            PI=4.D0*DATAN(1.D0)
            R=(2.0_rp*invar(1)*invar(1)*invar(1) - 9.0_rp*invar(1)*invar(2) + 27.0_rp*invar(3))/54.0_rp
            theta=acos(R/sqrt(-Q**3))

            res(1)= 2.0_rp*sqrt(-Q)*cos( theta             /3.0_rp) + (1.0_rp/3.0_rp)*invar(1)
            res(2)= 2.0_rp*sqrt(-Q)*cos((theta + 2.0_rp*PI)/3.0_rp) + (1.0_rp/3.0_rp)*invar(1)
            res(3)= 2.0_rp*sqrt(-Q)*cos((theta + 4.0_rp*PI)/3.0_rp) + (1.0_rp/3.0_rp)*invar(1)

        else

            res = 0.0_rp

        endif

    end if


  end function

!---------------------------------------------------------
  !Gathers and interpolates

  subroutine displacementGather
     implicit none
     integer(ip) :: nd,itime

     call a%Mesh%GetNdime(nd)
     
     call e%gather(nd,eldisp(:,:,1),a%disp(:,:,1))
     call e%gather(nd,eldisp(:,:,2),a%disp(:,:,3))
     do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(nd,eldisp(:,:,itime),a%disp(:,:,itime+1)) 
     enddo

  end subroutine  

  subroutine InterpolateGpAccel
     implicit none
     integer(ip) :: nd

     call a%Mesh%GetNdime(nd)
     
     call e%interpg(nd,elaccel(:,:,1),gpaccel)
  end subroutine
  
  subroutine InterpolateGpDisplacements
     implicit none
     integer(ip) :: nd,itime

     call a%Mesh%GetNdime(nd)
     
     call e%interpg(nd,eldisp(:,:,1),gpdisp(:,1))
     call e%interpg(nd,eldisp(:,:,2),gpdisp(:,2))
     do itime = 3,nsteps ! Time bdf2 and others
        call e%interpg(nd,eldisp(:,:,itime),gpdisp(:,itime))
     enddo

  end subroutine

  subroutine AssembleRhs
      implicit none
      integer(ip) :: nd,pn

      call a%Mesh%GetNdime(nd)
      pn  = e%pnode

      ! Assembly elrhu to elrhs
      elrhs(1:nd,1:pn) = elrhu(1:nd,1:pn) + elrhs(1:nd,1:pn) 
  end subroutine AssembleRhs

  subroutine AssembleMassMatrix
      implicit none
      integer(ip) :: nd,pn
      integer(ip) :: u1,uf,s1,sf,p1,bc

      call a%Mesh%GetNdime(nd)
      pn  = e%pnode

      call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      !Assembly of jacobian A(d): mass 
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn) + massmat(1:nd,1:pn,1:nd,1:pn)
      !-dNdt*L*(V)*tau + (r-tau_static*tau_dyn*r)*(V+tau*L*(V))
      !du_sgs/dt=(r-tau_static*tau_dyn*r)

      !du_sgs/dt,V = 0 para oss, se va todo (r-tau_static*tau_dyn*r)*(V+tau*L*(V))
      !du/dt en el residuo no va para oss
  end subroutine AssembleMassMatrix

  subroutine AssembleStiffnessMatrix
      implicit none
      integer(ip) :: nd,pn

      call a%Mesh%GetNdime(nd)
      pn  = e%pnode

      !Assembly of jacobian A(d): stiffness 
      elmat(1:nd,1:pn,1:nd,1:pn) = elmat(1:nd,1:pn,1:nd,1:pn) + kmat(1:nd,1:pn,1:nd,1:pn)
  end subroutine AssembleStiffnessMatrix

  subroutine AssembleGeometricMatrix
      implicit none
      integer(ip) :: nd,pn

      call a%Mesh%GetNdime(nd)
      pn  = e%pnode

      !Assembly of jacobian A(d): geometric mats 
      elmat(1:nd,1:pn,1:nd,1:pn) = elmat(1:nd,1:pn,1:nd,1:pn) + geomat(1:nd,1:pn,1:nd,1:pn)
  end subroutine AssembleGeometricMatrix

  !----------Resets-------------------

  subroutine ResetDynamicMatrix
     implicit none
      massmat  = 0.0_rp
  end subroutine ResetDynamicMatrix

  subroutine ResetSolidBaseArrays
      implicit none
      kmat     = 0.0_rp
  end subroutine ResetSolidBaseArrays

  subroutine ResetNonLinearSolidArrays
      geomat   = 0.0_rp
  end subroutine ResetNonLinearSolidArrays

  subroutine ResetNRSolidArrays
      intForce = 0.0_rp
  end subroutine ResetNRSolidArrays

  subroutine ResetAssemblyMatrices
     implicit none
      elrhu =  0.0_rp
      elrhs =  0.0_rp
      elmat =  0.0_rp
  end subroutine ResetAssemblyMatrices

!---------------------------------------------------------
  !Allocs and Deallocs

  subroutine AllocateSolidBase
      implicit none
      integer(ip) :: nd,mn,df,sz

      call a%Mesh%GetNdime(nd)
      sz  = (nd*(nd+1))/2
      mn  = e%mnode
      df  = a%ndofn

      call a%Memor%alloc(df,mn,df,mn,kmat      ,'kmat'      ,'sld_Elmope')
      call a%Memor%alloc(sz,sz,      c_elas    ,'c_elas'    ,'sld_Elmope')
      call a%Memor%alloc(sz,         stress    ,'stress'    ,'sld_Elmope')
      call a%Memor%alloc(sz,         strain    ,'strain'    ,'sld_Elmope')
      call a%Memor%alloc(sz,         strain_aux,'strain_aux','sld_Elmope')
  end subroutine AllocateSolidBase

  subroutine AllocateNonLinearSolidArrays
      implicit none
      integer(ip) :: nd,mn,df

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode
      df  = a%ndofn

      call a%Memor%alloc(nd,nd,      B       ,'B'       ,'sld_Elmope') 
      call a%Memor%alloc(nd,nd,      F_mat   ,'F_mat'   ,'sld_Elmope')
      call a%Memor%alloc(nd,nd,      F_matt  ,'F_matt'  ,'sld_Elmope')
      call a%Memor%alloc(nd,nd,      F_spat  ,'F_spat'  ,'sld_Elmope')
      call a%Memor%alloc(nd,nd,      tensor  ,'tensor'  ,'sld_Elmope')
      call a%Memor%alloc(df,mn,df,mn,geomat  ,'geomat'  ,'sld_Elmope')
      call a%Memor%alloc(nd,nd,      piolaK2 ,'piolaK2' ,'sld_Elmope')
  end subroutine AllocateNonLinearSolidArrays

  subroutine AllocateNRSolidArrays
      implicit none
      integer(ip) :: nd,mn

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      call a%Memor%alloc(nd,mn,      intForce,'intForce','sld_Elmope')
  end subroutine AllocateNRSolidArrays

  subroutine AllocateDynamicArray
      implicit none
      integer(ip) :: uf,mn

      uf  = a%udofn
      mn  = e%mnode

      call a%Memor%alloc(uf,mn,uf,mn,massmat,'massmat','sld_Elmope')
  end subroutine AllocateDynamicArray

  subroutine AllocateAssemblyArrays
      implicit none
      integer(ip) :: nd,df,mn

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode
      df  = a%ndofn

      call a%Memor%alloc(df,mn,df,mn,elmat,'elmat','sld_Elmope')
      call a%Memor%alloc(df,mn,      elrhs,'elrhs','sld_Elmope')
      call a%Memor%alloc(nd,mn,      elrhu,'elrhu','sld_Elmope')
  end subroutine AllocateAssemblyArrays

  subroutine AllocateBaseElmopeArrays
      implicit none
      integer(ip) :: nd,nc,mn

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode
      nc  = a%ncomp - 1
      det = 1.0_rp
      nsteps = 2

      call a%Memor%alloc(nd,nd,   ii      ,'ii'      ,'sld_Elmope')
      call a%Memor%alloc(nd,      invar   ,'invar'   ,'sld_Elmope')
      call a%Memor%alloc(nd,      elext   ,'elext'   ,'sld_Elmope')
      call a%Memor%alloc(nd,mn,   elsou   ,'elsou'   ,'sld_Elmope')
      call a%Memor%alloc(nd,      gpsou   ,'gpsou'   ,'sld_Elmope')
      call a%Memor%alloc(nd,mn,nc,eldisp  ,'eldisp'  ,'sld_Elmope')
      call a%Memor%alloc(nd,nc,   gpdisp  ,'gpdisp'  ,'sld_Elmope')
      call a%Memor%alloc(nd,nd,   gradDisp,'gradDisp','sld_Elmope')

  end subroutine
  
  subroutine DeallocateSolidBase
      implicit none
      integer(ip) :: nd,mn,df,sz

      call a%Mesh%GetNdime(nd)
      sz  = (nd*(nd+1))/2
      mn  = e%mnode
      df  = a%ndofn

      call a%Memor%dealloc(df,mn,df,mn,kmat      ,'kmat'      ,'sld_Elmope')
      call a%Memor%dealloc(sz,sz,      c_elas    ,'c_elas'    ,'sld_Elmope')
      call a%Memor%dealloc(sz,         stress    ,'stress'    ,'sld_Elmope')
      call a%Memor%dealloc(sz,         strain    ,'strain'    ,'sld_Elmope')
      call a%Memor%dealloc(sz,         strain_aux,'strain_aux','sld_Elmope')
  end subroutine DeallocateSolidBase

  subroutine DeallocateNonLinearSolidArrays
      implicit none
      integer(ip) :: nd,mn,df

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode
      df  = a%ndofn

      call a%Memor%dealloc(nd,nd,      B       ,'B'       ,'sld_Elmope')
      call a%Memor%dealloc(nd,nd,      F_mat   ,'F_mat'   ,'sld_Elmope')
      call a%Memor%dealloc(nd,nd,      F_matt  ,'F_matt'  ,'sld_Elmope')
      call a%Memor%dealloc(nd,nd,      F_spat  ,'F_spat'  ,'sld_Elmope')
      call a%Memor%dealloc(nd,nd,      tensor  ,'tensor'  ,'sld_Elmope')
      call a%Memor%dealloc(df,mn,df,mn,geomat  ,'geomat'  ,'sld_Elmope')
      call a%Memor%dealloc(nd,nd,      piolaK2 ,'piolaK2' ,'sld_Elmope')
  end subroutine DeallocateNonLinearSolidArrays

  subroutine DeallocateNRSolidArrays
      implicit none
      integer(ip) :: nd,mn

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      call a%Memor%dealloc(nd,mn,      intForce,'intForce','sld_Elmope')
  end subroutine DeallocateNRSolidArrays

  subroutine DeallocateDynamicArray
      implicit none
      integer(ip) :: uf,mn

      uf  = a%udofn
      mn  = e%mnode

      call a%Memor%dealloc(uf,mn,uf,mn,massmat,'massmat','sld_Elmope')
  end subroutine DeallocateDynamicArray

  subroutine DeallocateAssemblyArrays
      implicit none
      integer(ip) :: nd,df,mn

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode
      df  = a%ndofn

      call a%Memor%dealloc(df,mn,df,mn,elmat,'elmat','sld_Elmope')
      call a%Memor%dealloc(df,mn,      elrhs,'elrhs','sld_Elmope')
      call a%Memor%dealloc(nd,mn,      elrhu,'elrhu','sld_Elmope')
  end subroutine DeallocateAssemblyArrays

  subroutine DeallocateBaseElmopeArrays
      implicit none
      integer(ip) :: nd,nc,mn

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode
      nc  = a%ncomp - 1

      call a%Memor%dealloc(nd,nd,   ii      ,'ii'      ,'sld_Elmope')
      call a%Memor%dealloc(nd,      invar   ,'invar'   ,'sld_Elmope')
      call a%Memor%dealloc(nd,      elext   ,'elext'   ,'sld_Elmope')
      call a%Memor%dealloc(nd,mn,   elsou   ,'elsou'   ,'sld_Elmope')
      call a%Memor%dealloc(nd,      gpsou   ,'gpsou'   ,'sld_Elmope')
      call a%Memor%dealloc(nd,mn,nc,eldisp  ,'eldisp'  ,'sld_Elmope')
      call a%Memor%dealloc(nd,nc,   gpdisp  ,'gpdisp'  ,'sld_Elmope')
      call a%Memor%dealloc(nd,nd,   gradDisp,'gradDisp','sld_Elmope')

  end subroutine

  subroutine AllocGpAcceleration
     implicit none
     integer(ip) :: nd
      call a%Memor%alloc(nd,  gpaccel,'gpaccel'  ,'sld_Elmope')
  end subroutine

  subroutine DeallocGpAcceleration
     implicit none
     integer(ip) :: nd
      call a%Memor%dealloc(nd,  gpaccel,'gpaccel'  ,'sld_Elmope')
  end subroutine

  subroutine AllocVelandAccel
     implicit none
     integer(ip) :: nd,mn

     call a%Mesh%GetNdime(nd)
     mn  = e%mnode

     call a%Memor%alloc(nd,mn,2,elveloc,'elveloc','sld_Elmope')
     call a%Memor%alloc(nd,mn,2,elaccel,'elaccel','sld_Elmope')
  end subroutine

  subroutine GatherDynamic
     implicit none
     integer(ip) :: nd

     call a%Mesh%GetNdime(nd)

     call e%gather(nd,elveloc(:,:,1),a%veloc(:,:,1))
     call e%gather(nd,elaccel(:,:,1),a%accel(:,:,1))
     call e%gather(nd,elveloc(:,:,2),a%veloc(:,:,2))
     call e%gather(nd,elaccel(:,:,2),a%accel(:,:,2))
  end subroutine

  subroutine DeallocVelandAccel
     implicit none
     integer(ip) :: nd,mn

     call a%Mesh%GetNdime(nd)
     mn  = e%mnode

     call a%Memor%dealloc(nd,mn,2,elveloc,'elveloc','sld_Elmope')
     call a%Memor%dealloc(nd,mn,2,elaccel,'elaccel','sld_Elmope')
  end subroutine
!---------------------------------------------------------
  !Newmark
  subroutine InterpolateNewmark
     implicit none
     integer(ip) :: nd

     call a%Mesh%GetNdime(nd)

     call e%interpg(nd,elveloc(:,:,2),gpdisp(:,3))!gpveloc)
     call e%interpg(nd,elaccel(:,:,2),gpdisp(:,4))!gpaccel)

  end subroutine

end module Mod_sld_BaseElmope

#include "sldsup/Elmopes/Mod_sldsup_BaseElmope.f90"   
#include "sldup/Elmopes/Mod_sldup_BaseElmope.f90"   
