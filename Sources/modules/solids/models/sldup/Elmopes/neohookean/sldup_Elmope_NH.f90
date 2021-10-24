module Mod_sldup_elmope_NH
  use Mod_sld_BaseElmopeRoutines
  use Mod_CauchyElement_tools
  use Mod_CauchyElement
  use Mod_UPSolids_NH
  use Mod_sldup_InterpolateResidualProjection
  implicit none   
  real(rp)    :: tau_u,tau_s,tau_p
  integer(ip) :: u1,uf,s1,sf,p1,bc

   
  contains

  subroutine SetPointersElmopeUP_NH
     use typre
     implicit none
     procedure() :: NULLSUB
     integer(ip) :: kfl_nonlinear,nelty
           
     call up_NH%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)

     call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)

     call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateAssemblyArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateAssemblyArrays)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUP)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUPNonLinear)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUPNonLinear)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesUP_SGS)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesUP_SGS)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateUPAssemblyMatrices)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateUPAssemblyMatrices)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateHigherOrderDerivativeArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateHigherOrderDerivativeArrays)

     call ConcatenateProcedures(ProcHook%ResetArrays,ResetUPBaseArrays)
     call ConcatenateProcedures(ProcHook%ResetArrays,ResetUPBaseArrays_SGS)

     call ConcatenateProcedures(ProcHook%ResetArrays,ResetForceVector)

     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpDisplacements)
     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpPress)

     call SetPointersExternalForces(0)
     call SetPointersExternalForces(1)
     call SetPointersExternalForces(100)
     call ConcatenateProcedures(ProcHook%Initializations,calculateExternalForcesUP)

     call SetPointers_NH

     call SetNRPointers

     call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixUP)
     call ConcatenateProcedures(ProcHook%preAssemblyRhs,AssembleRhsUP)
     call ConcatenateProcedures(ProcHook%AssemblyRhs,AssembleRhsForcesUP)

     !--------------------------------------------------------------------------
     if(up_NH%kfl_timei==1) then
         !Dynamic allocations

         call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateDynamicArray)
         call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateDynamicArray)

         !Dynamic resets
         call ConcatenateProcedures(ProcHook%ResetArrays,ResetDynamicMatrix)

         !Dynamic elemental build
         call ConcatenateProcedures(ProcHook%DynamicMass,buildMassMat)
     
         !Dynamic build to LHS
         call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleMassMatrixUP)

         call ConcatenateProcedures(ProcHook%Initializations,AllocVelandAccel)
         call ConcatenateProcedures(ProcHook%Initializations,InitTimeIntegration)

         call ConcatenateProcedures(ProcHook%Gathers        ,GatherDynamic)

         if (up%kfl_tacsg == 1) then

             call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesUP_dynsgs )
             call ConcatenateProcedures(ProcHook%DeallocateArrays  ,DeallocateBaseElmopeMatricesUP_dynsgs )
             call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixUP_dynsgs)

             call ConcatenateProcedures(ProcHook%ResetArrays,ResetUPBaseArrays_dynsgs)

             !Dynamic subscales terms for LHS
             call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH_ASGS_dynsgs)
             call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH_dynsgs)
             if(kfl_nonlinear==1) then 
                 !call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH_dynsgs_nonli)
             endif

         endif

         call ConcatenateProcedures(ProcHook%Finalizations  ,DeallocVelandAccel)

         !Dynamic matrix routines
         call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicMass)
         call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicForce)

     endif

     call ConcatenateProcedures(ProcHook%ToLinearSystem,AssembleLinearSystem)
     !-------------------------------------------------------
     !If more than one element type, then pointers should be reset when ielty changes!
     call a%Mesh%GetNelty(nelty)
     if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

    !Mass matrix
    call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%Dynamic)

  end subroutine SetPointersElmopeUP_NH
     
  subroutine SetPointers_NH
     use typre
     implicit none
     integer(ip) :: kfl_nonlinear
     call up_NH%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)

     call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp)

     call ConcatenateProcedures(ProcHook%InGauss,calculateUPGradients)
     if (kfl_nonlinear == 1) then 
         call ConcatenateProcedures(ProcHook%InGauss,calculateFmatGradient)
     endif

     call ConcatenateProcedures(ProcHook%Gathers,displacementGather)
     call ConcatenateProcedures(ProcHook%Gathers,pressGather)

     !GALERKIN
     call ConcatenateProcedures(ProcHook%PhysicalProp,calculateTauAndMatrixOrganization)

     call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_UP_NH)
     call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH)

     ProcPointer%sld_elmrhu => sld_elmrhu
     call ConcatenateProcedures(ProcHook%EndGauss,integrateForces)

     if (up%kfl_repro == 0) then
         !ASGS
         call ConcatenateProcedures(ProcHook%EndGauss,integrateForces_SGS)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH_ASGS)
         call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_UP_NH_ASGS)

         if(a%kfl_timei==1) then

             call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpAccel)
             call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH_ASGS_dyn)
             call ConcatenateProcedures(ProcHook%PostInterpolates,calculateAU_U_dyn)
             if (kfl_nonlinear == 1) then
                 call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH_ASGS_dyn_nonli)
             endif
         else
             call ConcatenateProcedures(ProcHook%PostInterpolates,calculateAU_U)
         endif

         call ConcatenateProcedures(ProcHook%PostPostInterpolates,completeAU)
         call ConcatenateProcedures(ProcHook%PostInterpolates,calculateAU_P)

     elseif (up%kfl_repro == 1) then

         !OSGS
         call SetPointersInterpolateResidualProjection(0)
         call SetPointersInterpolateResidualProjection(1)
         call SetPointersInterpolateResidualProjection(100)

         if(a%kfl_timei==1) then
             call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpAccel)
         endif

         call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_UP_NH_OSGS)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH_OSGS)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH_OSGS_repro)

         call ConcatenateProcedures(ProcHook%PrePostInterpolates,calculateAU_U)
         call ConcatenateProcedures(ProcHook%PrePostInterpolates,calculateAU_P_OSGS)
         call ConcatenateProcedures(ProcHook%PostInterpolates,completeAU)
         call ConcatenateProcedures(ProcHook%PostPostInterpolates,addReproToAU)

     endif

     !if (kfl_nonlinear == 1) then
     !    call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPNonlinearDerTerms)
     !    call ConcatenateProcedures(ProcHook%InGauss,calculateInGaussNonLinearTerms)
     !    call ConcatenateProcedures(ProcHook%EndGauss,PostGaussElmats_NH_SGS_nonlinear)
     !    call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_UP_NH_ASGS_nonli)
     !endif

  end subroutine SetPointers_NH

   subroutine SetNRPointers

     !Residual
     if(a%kfl_timei==1) then
         call ConcatenateProcedures(ProcHook%calculateResidual,dynamicResidualUP)
     endif

     !Elemental assemblies
     ProcPointer%getTauParameters => getTauParametersUP
     call ConcatenateProcedures(ProcHook%preAssemblyLhs,residual_and_jacobian)

     call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixUP_SGS)
     call ConcatenateProcedures(ProcHook%postAssemblyLhs,completeForcingTermUP)


   end subroutine SetNRPointers

   subroutine residual_and_jacobian
       implicit none

       call ProcHook%calculateResidual

       if(a%kfl_timei==1) then
           !We complete the mass matrix Jacobian
          massmat = tsch_deriv*a%dtinv*a%dtinv*massmat
      endif

   end subroutine
   
   subroutine addReproToAU

       AU = AU + gprep

   end subroutine addReproToAU

   subroutine completeForcingTermUP
      implicit none
      integer(ip)          :: pn,nd,tn,df,idime,i,count
      real(rp),allocatable :: var(:,:),aux(:,:)
      !Very useful formulation for ROM, instead of calculating residual 
      !we calculate the variable
      !by solving A*U_i+1 = r + A*U_i, ------> U_i+1 - U_i = d_u
      !Finally, the problem to solve is : U_i+1 = A^-1*(r + A*U_i)
      !Usual Newton Raphson assembly is: d_u = A^-1 r

      !We add the previous value so we do not calculate residual 
      !but the total variable

      call up%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode
      df  = a%ndofn

      call up%Memor%alloc(df,pn,var,'var','sldup_Elmope')
      call up%Memor%alloc(df,pn,aux,'aux','sldup_Elmope')

      do i=1,pn
          count = 1
          do idime=u1,uf
              var(idime,i) = eldisp(count,i,1)
              count = count + 1
          enddo
          var(p1,i)    = elpress(1,i)
      enddo

      call matvecmult(df,pn,elmat,var,aux)

      elrhs = elrhs + aux

      call up%Memor%dealloc(df,pn,var,'var','sldup_Elmope')
      call up%Memor%dealloc(df,pn,aux,'aux','sldup_Elmope')

   end subroutine completeForcingTermUP

  subroutine calculateTauAndMatrixOrganization
      implicit none

      call up_NH%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

  end subroutine


  subroutine calculateAU_P_OSGS
     implicit none
     integer(ip) :: nd,tn

     call up_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     call getAU_P_OSGS(nd,G,lam,det,B,gppress,AU(p1))

  end subroutine calculateAU_P_OSGS


  subroutine completeAU

      AU(u1:uf) = AU(u1:uf)/(2.0_rp*up%mu)

  end subroutine completeAU

  !----------------------------------------------------------
  !PostGauss Matrices

  subroutine GaussElmats_UP_NH
     implicit none
     integer(ip) :: nd,tn

     call up_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     call  up_elmVU_NH(e,nd,tn,dvol,det,G,gpsigma   ,F_mat,F_spat ,elmuv)
     call  up_elmVP_NH(e,nd,   dvol,G,                             elmpv)
     call sup_elmQU_NH(e,nd,   dvol,det,G,lam,F_spat,F_mat,gppress,elmuq)
     call sup_elmQP_NH(e,nd,   dvol,det,  lam,                     elmpq)

  end subroutine

   subroutine GaussElmats_NH_ASGS_dynsgs
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_u,tau_s,tau_p,dvol_mod

       call up_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       dvol_mod = dvol/(2.0_rp*up%mu)
       !-------------------Tau_u matrix calculation----------
       call sup_elm_usgs_VU_NH_dynsgs(e,nd,tn,dvol_mod,tautau,densi,a%dtinv2   ,elm_usgsdyn_VU)
       call  up_elm_usgs_VS_NH_dynsgs(e,nd,tn,dvol_mod,tautau,    divstr,F_spat,elm_usgsdyn_VU)
       call sup_elm_usgs_VP_NH_dynsgs(e,nd,tn,dvol_mod,tautau                  ,elm_usgsdyn_VP)

   end subroutine GaussElmats_NH_ASGS_dynsgs

   subroutine GaussElmats_NH_ASGS_dyn
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_u,tau_s,tau_p,dvol_mod

       call up_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------Tau_u matrix calculation----------
       dvol_mod = dvol/(2.0_rp*up%mu)
       call sup_elm_usgs_QU_NH_dyn(e,nd,tn,dvol_mod,det,G,lam,grdpre,gppress,F_mat,F_spat,densi,a%dtinv2,elm_usgs_QU)

   end subroutine GaussElmats_NH_ASGS_dyn

   subroutine GaussElmats_NH_ASGS_dyn_nonli
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_u,tau_s,tau_p

       call up_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------Tau_u matrix calculation----------
       !This term cancels out when using OSS
       !call up_elm_usgs_QU_NH_dyn_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,densi,a%dtinv2,elm_usgs_QU)

   end subroutine GaussElmats_NH_ASGS_dyn_nonli

  subroutine GaussElmats_UP_NH_ASGS
     implicit none
     integer(ip) :: nd,tn
     real(rp)    :: dvol_mod

     call up_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     dvol_mod = dvol/(2.0_rp*up%mu)

     !-------------------Tau_u matrix calculation----------
     call sup_elm_usgs_QP_NH(e,nd,tn,dvol_mod,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QP)

     !!-------------------Tau_p matrix calculation----------
     call sup_elm_psgs_VU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_VU)
     call sup_elm_psgs_VP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_VP)
     call sup_elm_psgs_QU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_QU)
     call sup_elm_psgs_QP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_QP)


  end subroutine GaussElmats_UP_NH_ASGS

  subroutine GaussElmats_UP_NH_OSGS
     implicit none
     integer(ip) :: nd,tn
     real(rp)    :: dvol_mod

     call up_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     dvol_mod = dvol/(2.0_rp*up%mu)

     !-------------------Tau_u matrix calculation----------
     call sup_elm_usgs_QP_NH(e,nd,tn,dvol_mod,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QP)

     !-------------------Tau_p matrix calculation----------
     call sup_elm_psgs_VU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_VU)
     !call sup_elm_psgs_VP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_VP)
     call sup_elm_psgs_QU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_QU)
     !call sup_elm_psgs_QP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_QP)

  end subroutine GaussElmats_UP_NH_OSGS

   subroutine calculateRHS_UP_NH
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: dvol_mod

       call up_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       dvol_mod = dvol/(2.0_rp*up%mu)
       call sup_rhs_U_NH(e,nd,tn,dvol_mod,               gpsigma,gppress,elrhu)
       call sup_rhs_P_NH(e,nd,tn,dvol,G,lam,det,B,        gppress,elrhp,divdisp)

   end subroutine

   subroutine calculateRHS_UP_NH_ASGS
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_pmod

       call up_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------U RHS calculation----------
       tau_pmod = tau_p/(2.0_rp*up%mu)
       call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_pmod,AU(p1),elrhu)
       !!-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,AU(u1:uf),elrhp)
       call sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,det,tau_p,AU(p1),elrhp)

   end subroutine calculateRHS_UP_NH_ASGS

  subroutine calculateRHS_UP_NH_OSGS
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: tau_pmod

      call up_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      !-------------------U RHS calculation----------
      tau_pmod = tau_p/(2.0_rp*up%mu)
      call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_pmod,AU(p1),elrhu)
      !-------------------P RHS calculation----------
      call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,AU(u1:uf),elrhp)
      call sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,det,tau_p,AU(p1),elrhp)

  end subroutine calculateRHS_UP_NH_OSGS

  subroutine calculateRHS_UP_NH_OSGS_repro
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: tau_pmod

      call up_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      !-------------ADJOINT TESTING OF PROJECTION -----------

      !-------------------U RHS calculation----------
      tau_pmod = tau_p/(2.0_rp*up%mu)
      call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_pmod,gprep(p1),elrhu)
      !-------------------P RHS calculation----------
      call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,gprep(u1:uf),elrhp)
      call sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,det,tau_p,gprep(p1),elrhp)

  end subroutine calculateRHS_UP_NH_OSGS_repro

  subroutine calculateRHS_UP_NH_dynsgs
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: dvol_mod,AU_aux(e%ndime+1)

      call up_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      dvol_mod = dvol/(2.0_rp*up%mu)

      !-------------------U RHS calculation----------
      !(Vh*(1-tau_k^-1*tau_t),elext-AU)
      call sup_rhs_usgs_U_NH_AU_dynsgs(e,nd,tn,dvol,tautau,elext,AU(u1:uf),elrhu)

      !subscale past acceleration
      call sup_rhs_usgs_U_NH_SGSACCEL_dynsgs(e,nd,tn,dvol_mod,a%dtinv2,densi,tautau,up_NH%u_sgs(ielem)%a(:,1,igaus),elrhu)
      !-------------------P RHS calculation----------
      call sup_rhs_usgs_P_NH_dynsgs(e,nd,tn,dvol_mod,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,a%dtinv2,densi,up_NH%u_sgs(ielem)%a(:,1,igaus),elrhp)

  end subroutine calculateRHS_UP_NH_dynsgs

  subroutine calculateRHS_UP_NH_dynsgs_nonli
      implicit none
      integer(ip) :: nd,tn

      call up_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      !-------------------U RHS calculation----------
      !call up_rhs_usgs_U_NH(e,nd,tn,dvol,tau_u,AU(u1:uf),elrhu)
      !-------------------P RHS calculation----------
      !call up_rhs_usgs_P_NH_dynsgs_nonli(e,nd,tn,dvol,G,lam,det,tau_u,gppress,derJFki,divFmat,grdFspt,a%dtinv2,densi,up_NH%u_sgs(ielem)%a(:,1,igaus),elrhp)

  end subroutine calculateRHS_UP_NH_dynsgs_nonli

  subroutine calculateInGaussNonLinearTerms
      implicit none
      integer(ip) :: nd,tn

      call up_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call sld_calculateJacobianDer_Fspat_ki(e,nd,tn,det,F_spat,grdFmat,derJFki)

  end subroutine calculateInGaussNonLinearTerms

  subroutine PostGaussElmats_NH_SGS_nonlinear
     implicit none
     integer(ip) :: nd,tn

     call up_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     !-------------------Tau_u matrix calculation----------
     !TODO
     !call up_elm_usgs_VU_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QS)
     !call up_elm_usgs_VP_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QP)

     !call sup_elm_usgs_QU_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QU)
     call sup_elm_usgs_QP_NH_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,elm_usgs_QP)

  end subroutine

   subroutine calculateRHS_UP_NH_ASGS_nonli
       implicit none
       integer(ip) :: nd,tn

       call up_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------U RHS calculation----------
       !TODO

       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,tau_u,-elext,up_NH%extForce(ielem)%a(p1,:))

   end subroutine calculateRHS_up_NH_ASGS_nonli

  subroutine calculateExternalForcesUP
      implicit none

      !Compute contributions to RHS :
      elext = 0.0_rp
      call ProcPointer%ExternalForces   
      elext = elext/(2.0_rp*up%mu)

  end subroutine

  subroutine integrateForces
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: dvol_mod

      call up_NH%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2

      !Compute contributions to RHS :
      call ProcPointer%sld_elmrhu(e,dvol,elext,up_NH%extForce(ielem)%a(u1:uf,:))

  end subroutine

  subroutine integrateForces_SGS
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: dvol_mod

      call up_NH%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2

      call sup_rhs_usgs_P_NH_extforces(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,up_NH%extForce(ielem)%a(p1,:),elext)

  end subroutine

   subroutine ResetForceVector
         implicit none
   
      up_NH%extForce(ielem)%a(:,:)=0.0_rp
   
   end subroutine

  subroutine PhysicalProp
      implicit none
      integer(ip) :: nd,tn

      call up_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call up_NH%UPGetPhysicalParameters(nd,tn,densi,lam,G)
  end subroutine

  subroutine OnIeltyChange
     implicit none
      
     if (e%ielty /= ielty0) then
        call SetPointersElmopeUP_NH
        ielty0 = e%ielty
     endif
  end subroutine

   subroutine lhsrhscheck_up
      implicit none
      integer(ip)          :: pn,nd,tn,df,idime,i,count,mn
      real(rp),allocatable :: var(:,:),aux(:,:),elLhs(:,:),realRhs(:,:),realForcesRhs(:,:)
      !Very useful formulation for ROM, instead of calculating residual 
      !we calculate the variable
      !by solving A*U_i+1 = r + A*U_i, ------> U_i+1 - U_i = d_u
      !Finally, the problem to solve is : U_i+1 = A^-1*(r + A*U_i)
      !Usual Newton Raphson assembly is: d_u = A^-1 r

      !We add the previous value so we do not calculate residual 
      !but the total variable

      call up%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode
      df  = a%ndofn
      mn  = e%mnode

      call up%Memor%alloc(df,pn,var,'var','sldup_Elmope')
      call up%Memor%alloc(df,pn,aux,'aux','sldup_Elmope')
      call up%Memor%alloc(df,pn,elLhs,'elLhs','sldup_Elmope')
      call up%Memor%alloc(df,pn,realRhs,'realRhs','sldup_Elmope')
      call up%Memor%alloc(df,pn,realForcesRhs,'realForcesRhs','sldup_Elmope')

      do i=1,pn
          count = 1
          do idime=u1,uf
              var(idime,i) = eldisp(count,i,1)
              count = count + 1
          enddo
          count = 1
          var(p1,i)    = elpress(1,i)
      enddo

      call matvecmult(df,pn,elmat,var,aux)

      !Compare with Totalview, ideally both elLhs and realRhs should be the same
      realForcesRhs(u1:uf,1:pn) = realForcesRhs(u1:uf,1:pn) + a%extForce(ielem)%a(u1:uf,:) 
      realForcesRhs(p1,1:pn)    = realForcesRhs(p1,1:pn)    + a%extForce(ielem)%a(p1,:)    

      realForcesRhs = realForcesRhs + elrhs
      elLhs         = -aux
      realRhs       =  elrhs

      call up%Memor%dealloc(df,pn,var,'var','sldup_Elmope')
      call up%Memor%dealloc(df,pn,aux,'aux','sldup_Elmope')
      call up%Memor%dealloc(df,pn,elLhs,'elLhs','sldup_Elmope')
      call up%Memor%dealloc(df,pn,realRhs,'realRhs','sldup_Elmope')
      call up%Memor%dealloc(df,pn,realForcesRhs,'realForcesRhs','sldup_Elmope')

   end subroutine lhsrhscheck_up

end module   

!Solids_ELMOPE subroutine   
subroutine sldup_Elmope_NH(SldProblem)
   use Mod_sldup_elmope_NH

   implicit none
   class(UPSolidsProblem_NH), target :: SldProblem
   integer(ip):: ndi
   
   !Setting sld problem pointer in sld_BaseElmope
   up_NH  => SldProblem
   up     => SldProblem
   a      => SldProblem

   call up%Mesh%GetNdime(ndi)

   call SetPointersGeneral
   call SetPointersPointForces
   call SetPointersElmopeUP_NH
   
   !TODO: there is a memory fault somewhere, it jumps from the AllocateArrays to gpaccel
   call a%Memor%alloc(ndi,  gpaccel,'gpaccel'  ,'sld_Elmope')
   call sld_elemLoop
   call a%Memor%dealloc(nd,  gpaccel,'gpaccel'  ,'sld_Elmope')

end subroutine sldup_Elmope_NH
