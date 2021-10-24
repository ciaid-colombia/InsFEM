module Mod_sld_begiteElmopeSUP
   use Mod_CauchyElement
   use Mod_SUPSolids
   use Mod_SUPSolids_NH
   use Mod_sld_BaseElmopeRoutines
   use Mod_sldsup_elmope_NH
   implicit none   
   
   contains

  ! subroutine SetPointersBegiteSUP
  !    use typre
  !    implicit none
  !    integer(ip) :: kfl_nonlinear,nelty

  !    call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateSUPAssemblyMatrices)
  !    call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateSUPAssemblyMatrices)

  !    call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesSUP)
  !    call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP)

  !    call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesSUPNonLinear)
  !    call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUPNonLinear)

  !    call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays)

  !    call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp)

  !    call ConcatenateProcedures(ProcHook%Gathers,pressGather)
  !    call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpPress)
  !          
  !    call ConcatenateProcedures(ProcHook%PostInterpolates,calculateForcesSGS)

  !    call ConcatenateProcedures(ProcHook%AssemblyRhs,AssembleForceRhsSGS)

  !end subroutine SetPointersBegiteSUP

  !subroutine calculateForcesSGS
  !    implicit none
  !    integer(ip) :: nd,tn
  !    real(rp)    :: tau_u,tau_s,tau_p

  !    call sup_NH%Mesh%GetNdime(nd)
  !    tn = (nd*(nd+1))/2

  !    call getTauParameters(tau_u,tau_s,tau_p)

  !    call ProcHook%PhysicalProp
  !    call calculateGradientsAndDeter
  !    call calculatePressGradientsAndStressDiv

  !    !Compute contributions to RHS :
  !    call sup_rhs_usgs_S_NH_extforces(e,nd,tn,dvol,G,det,tau_u,F_mat,F_spat,grsig,gpsigma,elrhd,elext)
  !    call sup_rhs_usgs_P_NH_extforces(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,elrhp,elext)

  !end subroutine

  ! subroutine AssembleForceRhsSGS
  !     implicit none
  !     integer(ip) :: u1,uf,s1,sf,p1,bc

  !     call sup_NH%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

  !     sup_NH%extForce(ielem)%a(u1:uf,:) = elrhu(:,:)
  !     sup_NH%extForce(ielem)%a(s1:sf,:) = elrhd(:,:) 
  !     sup_NH%extForce(ielem)%a(p1,:)    = elrhp(1,:) 

  ! end subroutine

end module Mod_sld_begiteElmopeSUP

subroutine sldsup_begite(SldProblem)
!-----------------------------------------------------------------------
! NAME 
!    sldsup_begite
! DESCRIPTION
!    This routine starts an internal iteration for the elastic solid prob. 
!-----------------------------------------------------------------------
   use typre
   use Mod_sld_begiteElmope
   !use Mod_sld_begiteElmopeSUP
   implicit none
   class(SUPSolidsProblem_NH), target :: SldProblem

   sup_NH => SldProblem
   sup    => SldProblem
   a      => SldProblem

   !Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
   sup_NH%disp(:,:,1)  = sup_NH%disp(:,:,2)
   sup_NH%sigma(:,:,1) = sup_NH%sigma(:,:,2)
   sup_NH%press(:,1)   = sup_NH%press(:,2)

   !call SetPointersBegite
   !call SetPointersBegiteSUP

   !call sld_elemLoop

end subroutine sldsup_begite
