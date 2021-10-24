module Mod_sldsup_elmope_NH
  use Mod_sldsup_calculateAU
  use Mod_CauchyElement_tools
  use Mod_CauchyElement
  use Mod_SUPSolids_NH
  use Mod_sldsup_InterpolateResidualProjection
  use Mod_sld_BaseElmopeRoutines
  implicit none   
  real(rp)    :: tau_u,tau_s,tau_p
  integer(ip) :: u1,uf,s1,sf,p1,bc
   
  contains

  subroutine SetPointersElmopeSUP_NH
     use typre
     implicit none
     procedure() :: NULLSUB
     integer(ip) :: kfl_nonlinear,nelty
           
     call sup_NH%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)

     call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)

     call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateAssemblyArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateAssemblyArrays)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUP)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUPNonLinear)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUPNonLinear)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUP_SGS)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP_SGS)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateSUPAssemblyMatrices)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateSUPAssemblyMatrices)

     call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateHigherOrderDerivativeArrays)
     call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateHigherOrderDerivativeArrays)

     call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays)
     call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays_SGS)

     call ConcatenateProcedures(ProcHook%ResetArrays,ResetForceVector)

     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpDisplacements)
     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpSigma)
     call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpPress)

     call SetPointersExternalForces(0)
     call SetPointersExternalForces(1)
     call SetPointersExternalForces(100)
     call ConcatenateProcedures(ProcHook%Initializations,calculateExternalForces)

     call SetPointers_NH

     call SetNRPointers

     call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixSUP)
     call ConcatenateProcedures(ProcHook%preAssemblyRhs,AssembleRhsSUP)
     call ConcatenateProcedures(ProcHook%AssemblyRhs,AssembleRhsForcesSUP)

     !--------------------------------------------------------------------------
     if(sup_NH%kfl_timei==1) then
         !Dynamic allocations

         call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesSUP_SGS_dyn)
         call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP_SGS_dyn)

         call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateDynamicArray)
         call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateDynamicArray)

         !call ConcatenateProcedures(ProcHook%AllocateArrays,AllocGpAcceleration)
         !call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocGpAcceleration)

         !Dynamic resets
         call ConcatenateProcedures(ProcHook%ResetArrays,ResetDynamicMatrix)
         call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays_SGS_dyn)

         call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixSUP_SGS_dyn)

         !Dynamic elemental build
         call ConcatenateProcedures(ProcHook%DynamicMass,buildMassMat)
     
         !Dynamic build to LHS
         call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleMassMatrix)

         call ConcatenateProcedures(ProcHook%Initializations,AllocVelandAccel)
         call ConcatenateProcedures(ProcHook%Initializations,InitTimeIntegration)

         call ConcatenateProcedures(ProcHook%Gathers        ,GatherDynamic)

         if (sup%kfl_tacsg == 1) then

             !Dynamic subscales terms for LHS
             call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesSUP_SGS_dynsgs)
             call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP_SGS_dynsgs)
             call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPBaseArrays_SGS_dynsgs)
             call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH_ASGS_dynsgs)
             call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH_dynsgs)
             if(kfl_nonlinear==1) then 
                 call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH_dynsgs_nonli)
             endif

             call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixSUP_SGS_dynsgs)

         endif

         call ConcatenateProcedures(ProcHook%Finalizations  ,DeallocVelandAccel)

         !Dynamic matrix routines
         call ConcatenateProcedures(ProcHook%Dynamic,ProcHook%DynamicMass)

     endif

     call ConcatenateProcedures(ProcHook%ToLinearSystem,AssembleLinearSystem)
     !-------------------------------------------------------
     !If more than one element type, then pointers should be reset when ielty changes!
     call a%Mesh%GetNelty(nelty)
     if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

    !Mass matrix
    call ConcatenateProcedures(ProcHook%EndGauss,ProcHook%Dynamic)

  end subroutine SetPointersElmopeSUP_NH
     
  subroutine SetPointers_NH
     use typre
     implicit none
     integer(ip) :: kfl_nonlinear
     call sup_NH%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)

     call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp)


     call ConcatenateProcedures(ProcHook%InGauss,calculateGradients)
     if (kfl_nonlinear == 1) then 
         call ConcatenateProcedures(ProcHook%InGauss,calculateFmatGradient)
     endif

     call ConcatenateProcedures(ProcHook%Gathers,displacementGather)
     call ConcatenateProcedures(ProcHook%Gathers,sigmaGather)
     call ConcatenateProcedures(ProcHook%Gathers,pressGather)

     !GALERKIN
     call ConcatenateProcedures(ProcHook%PhysicalProp,calculateTauAndMatrixOrganization)

     call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH)
     call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH)

     ProcPointer%sld_elmrhu => sld_elmrhu
     call ConcatenateProcedures(ProcHook%EndGauss,integrateForces)

     if (sup%kfl_repro == 0) then
         !ASGS
         call ConcatenateProcedures(ProcHook%EndGauss,integrateForces_SGS)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH_ASGS)
         call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH_ASGS)

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

         call ConcatenateProcedures(ProcHook%PostInterpolates,calculateAU_SP)

     elseif (sup%kfl_repro == 1) then
         !OSGS
         call SetPointersInterpolateResidualProjection(0)
         call SetPointersInterpolateResidualProjection(1)
         call SetPointersInterpolateResidualProjection(100)

         if(a%kfl_timei==1) then
             call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpAccel)
         endif

         call ConcatenateProcedures(ProcHook%PostInterpolates,calculateAU_U)
         call ConcatenateProcedures(ProcHook%PostInterpolates,calculateAU_SP_OSGS)

         call ConcatenateProcedures(ProcHook%EndGauss,GaussElmats_NH_OSGS)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH_OSGS)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH_OSGS_repro)

         call ConcatenateProcedures(ProcHook%PostInterpolates,addReproToAU)

     endif


     if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%ResetArrays,ResetSUPNonlinearDerTerms)
         call ConcatenateProcedures(ProcHook%InGauss,calculateInGaussNonLinearTerms)
         call ConcatenateProcedures(ProcHook%EndGauss,PostGaussElmats_NH_SGS_nonlinear)
         call ConcatenateProcedures(ProcHook%EndGauss,calculateRHS_SUP_NH_ASGS_nonli)
     endif

  end subroutine SetPointers_NH

   subroutine SetNRPointers

     !Residual
     if(a%kfl_timei==1) then
         call ConcatenateProcedures(ProcHook%calculateResidual,dynamicResidual)
     endif

     !Elemental assemblies
     ProcPointer%getTauParameters => getTauParameters
     call ConcatenateProcedures(ProcHook%preAssemblyLhs,residual_and_jacobian)

     call ConcatenateProcedures(ProcHook%AssemblyLhs,AssembleStiffnessMatrixSUP_SGS)
     call ConcatenateProcedures(ProcHook%postAssemblyLhs,completeForcingTermSUP)


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

   subroutine completeForcingTermSUP
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

      call sup%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode
      df  = a%ndofn

      call sup%Memor%alloc(df,pn,var,'var','sldsup_Elmope')
      call sup%Memor%alloc(df,pn,aux,'aux','sldsup_Elmope')

      do i=1,pn
          count = 1
          do idime=u1,uf
              var(idime,i) = eldisp(count,i,1)
              count = count + 1
          enddo
          count = 1
          do idime=s1,sf
              var(idime,i) = elsigma(count,i)
              count = count + 1
          enddo
          var(p1,i)    = elpress(1,i)
      enddo

      call matvecmult(df,pn,elmat,var,aux)

      elrhs = elrhs + aux

      call sup%Memor%dealloc(df,pn,var,'var','sldsup_Elmope')
      call sup%Memor%dealloc(df,pn,aux,'aux','sldsup_Elmope')

   end subroutine completeForcingTermSUP

  !----------------------------------------------------------
  !PostGauss Matrices

  subroutine calculateTauAndMatrixOrganization
      implicit none

      call sup_NH%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

  end subroutine

  subroutine calculateAU_SP
     implicit none
     integer(ip) :: nd,tn

     call sup_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     call getAU_S(nd,tn,G,det,B,gpsigma,AU(s1:sf))
     call getAU_P(nd,G,lam,det,B,gppress,AU(p1))

  end subroutine calculateAU_SP

  subroutine GaussElmats_NH
     implicit none
     integer(ip) :: nd,tn

     call sup_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     call sup_elmVS_NH(e,nd,tn,dvol,                               elmsv)
     call sup_elmVP_NH(e,nd,   dvol,                               elmpv)
     call sup_elmEU_NH(e,nd,tn,dvol,det,G,    F_spat,F_mat,gpsigma,elmut)
     call sup_elmES_NH(e,nd,tn,dvol,det,G,                         elmst)
     call sup_elmQU_NH(e,nd,   dvol,det,G,lam,F_spat,F_mat,gppress,elmuq)
     call sup_elmQP_NH(e,nd,   dvol,det,  lam,                     elmpq)

  end subroutine

   subroutine GaussElmats_NH_ASGS_dynsgs
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_u,tau_s,tau_p

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------Tau_u matrix calculation----------
       call sup_elm_usgs_VU_NH_dynsgs(e,nd,tn,dvol,tautau,densi,a%dtinv2,elm_usgs_VU)
       call sup_elm_usgs_VS_NH_dynsgs(e,nd,tn,dvol,tautau,elm_usgs_VS)
       call sup_elm_usgs_VP_NH_dynsgs(e,nd,tn,dvol,tautau,elm_usgs_VP)

   end subroutine GaussElmats_NH_ASGS_dynsgs

   subroutine GaussElmats_NH_ASGS_dyn
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_u,tau_s,tau_p

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------Tau_u matrix calculation----------
       !This term cancels out when using OSS
       call sup_elm_usgs_EU_NH_dyn(e,nd,tn,dvol,det,G,    grsig,gpsigma,F_mat,F_spat,densi,a%dtinv2,elm_usgs_EU)
       call sup_elm_usgs_QU_NH_dyn(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,densi,a%dtinv2,elm_usgs_QU)

   end subroutine GaussElmats_NH_ASGS_dyn

   subroutine GaussElmats_NH_ASGS_dyn_nonli
       implicit none
       integer(ip) :: nd,tn
       real(rp)    :: tau_u,tau_s,tau_p

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------Tau_u matrix calculation----------
       !This term cancels out when using OSS
       call sup_elm_usgs_EU_NH_dyn_nonli(e,nd,tn,dvol,det,G,gpsigma,derJFki,divFmat,grdFspt,densi,a%dtinv2,elm_usgs_EU)
       call sup_elm_usgs_QU_NH_dyn_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,densi,a%dtinv2,elm_usgs_QU)

   end subroutine GaussElmats_NH_ASGS_dyn_nonli

  subroutine GaussElmats_NH_ASGS
     implicit none
     integer(ip) :: nd,tn

     call sup_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     !-------------------Tau_u matrix calculation----------
     call sup_elm_usgs_ES_NH(e,nd,tn,dvol,det,G,    grsig,gpsigma,F_mat,F_spat,elm_usgs_ES)
     call sup_elm_usgs_EP_NH(e,nd,tn,dvol,det,G,    grsig,gpsigma,F_mat,F_spat,elm_usgs_EP)
     call sup_elm_usgs_QS_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QS)
     call sup_elm_usgs_QP_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QP)

     !-------------------Tau_s matrix calculation----------
     call sup_elm_ssgs_VU_NH(e,nd,tn,dvol,det,G,    gpsigma,F_mat,F_spat,elm_ssgs_VU)
     call sup_elm_ssgs_VS_NH(e,nd,tn,dvol,det,G                         ,elm_ssgs_VS)
     call sup_elm_ssgs_EU_NH(e,nd,tn,dvol,det,G,    gpsigma,F_mat,F_spat,elm_ssgs_EU)
     call sup_elm_ssgs_ES_NH(e,nd,tn,dvol,det,G                         ,elm_ssgs_ES)

     !-------------------Tau_p matrix calculation----------
     call sup_elm_psgs_VU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_VU)
     call sup_elm_psgs_VP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_VP)
     call sup_elm_psgs_QU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_QU)
     call sup_elm_psgs_QP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_QP)


  end subroutine GaussElmats_NH_ASGS

  subroutine GaussElmats_NH_OSGS
     implicit none
     integer(ip) :: nd,tn

     call sup_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     !-------------------Tau_u matrix calculation----------
     call sup_elm_usgs_ES_NH(e,nd,tn,dvol,det,G,    grsig,gpsigma,F_mat,F_spat,elm_usgs_ES)
     call sup_elm_usgs_EP_NH(e,nd,tn,dvol,det,G,    grsig,gpsigma,F_mat,F_spat,elm_usgs_EP)
     call sup_elm_usgs_QS_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QS)
     call sup_elm_usgs_QP_NH(e,nd,tn,dvol,det,G,lam,grdpre,gppress,F_mat,F_spat,elm_usgs_QP)

     !-------------------Tau_s matrix calculation----------
     call sup_elm_ssgs_VU_NH(e,nd,tn,dvol,det,G,    gpsigma,F_mat,F_spat,elm_ssgs_VU)
     !call sup_elm_ssgs_VS_NH(e,nd,tn,dvol,det,G                         ,elm_ssgs_VS)
     call sup_elm_ssgs_EU_NH(e,nd,tn,dvol,det,G,    gpsigma,F_mat,F_spat,elm_ssgs_EU)
     !call sup_elm_ssgs_ES_NH(e,nd,tn,dvol,det,G                         ,elm_ssgs_ES)

     !-------------------Tau_p matrix calculation----------
     call sup_elm_psgs_VU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_VU)
     !call sup_elm_psgs_VP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_VP)
     call sup_elm_psgs_QU_NH(e,nd,   dvol,det,G,lam,gppress,F_mat,F_spat,elm_psgs_QU)
     !call sup_elm_psgs_QP_NH(e,nd,   dvol,det,  lam                     ,elm_psgs_QP)

  end subroutine GaussElmats_NH_OSGS

   subroutine calculateRHS_SUP_NH
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       call sup_rhs_U_NH(e,nd,tn,dvol,            gpsigma,gppress,elrhu)
       call sup_rhs_S_NH(e,nd,tn,dvol,G,    det,B,gpsigma        ,elrhd)
       call sup_rhs_P_NH(e,nd,tn,dvol,G,lam,det,B,        gppress,elrhp,divdisp)

   end subroutine

  subroutine calculateAU_SP_OSGS
     implicit none
     integer(ip) :: nd,tn

     call sup_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     call getAU_S_OSGS(nd,tn,G,B,gpsigma,AU(s1:sf))
     call getAU_P_OSGS(nd,G,lam,det,B,gppress,AU(p1))

  end subroutine calculateAU_SP_OSGS

   subroutine calculateRHS_SUP_NH_OSGS
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------U RHS calculation----------
       call sup_rhs_ssgs_U_NH(e,nd,tn,dvol,tau_s,AU(s1:sf),elrhu)
       call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_p,AU(p1),elrhu)
       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH(e,nd,tn,dvol,G,det,tau_u,grsig,gpsigma,F_mat,F_spat,AU(u1:uf),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,AU(u1:uf),elrhp)

   end subroutine calculateRHS_SUP_NH_OSGS

   subroutine calculateRHS_SUP_NH_ASGS
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------U RHS calculation----------
       call sup_rhs_ssgs_U_NH(e,nd,tn,dvol,tau_s,AU(s1:sf),elrhu)
       call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_p,AU(p1),elrhu)
       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH(e,nd,tn,dvol,G,det,tau_u,grsig,gpsigma,F_mat,F_spat,AU(u1:uf),elrhd)
       call sup_rhs_ssgs_S_NH(e,nd,tn,dvol,G,det,tau_s,AU(s1:sf),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,AU(u1:uf),elrhp)
       call sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,det,tau_p,AU(p1),elrhp)

   end subroutine calculateRHS_SUP_NH_ASGS

   subroutine calculateRHS_SUP_NH_dynsgs
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------U RHS calculation----------
       !(Vh*(1-tau_k^-1*tau_t),elext-AU)
       call sup_rhs_usgs_U_NH_AU_dynsgs(e,nd,tn,dvol,tautau,elext,AU(u1:uf),elrhu)
       !subscale past acceleration
       call sup_rhs_usgs_U_NH_SGSACCEL_dynsgs(e,nd,tn,dvol,a%dtinv2,densi,tautau,sup_NH%u_sgs(ielem)%a(:,1,igaus),elrhu)
       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH_dynsgs(e,nd,tn,dvol,G,det,tau_u,grsig,gpsigma,F_mat,F_spat,a%dtinv2,densi,sup_NH%u_sgs(ielem)%a(:,1,igaus),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH_dynsgs(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,a%dtinv2,densi,sup_NH%u_sgs(ielem)%a(:,1,igaus),elrhp)

   end subroutine

   subroutine calculateRHS_SUP_NH_dynsgs_nonli
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH_dynsgs_nonli(e,nd,tn,dvol,G,det,tau_u,gpsigma,derJFki,divFmat,grdFspt,a%dtinv2,densi,sup_NH%u_sgs(ielem)%a(:,1,igaus),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH_dynsgs_nonli(e,nd,tn,dvol,G,lam,det,tau_u,gppress,derJFki,divFmat,grdFspt,a%dtinv2,densi,sup_NH%u_sgs(ielem)%a(:,1,igaus),elrhp)

   end subroutine calculateRHS_SUP_NH_dynsgs_nonli

   subroutine calculateRHS_SUP_NH_ASGS_repro
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------ADJOINT TESTING OF PROJECTION -----------

       !-------------------U RHS calculation----------
       call sup_rhs_ssgs_U_NH(e,nd,tn,dvol,tau_s,gprep(s1:sf),elrhu)
       call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_p,gprep(p1),elrhu)
       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH(e,nd,tn,dvol,G,det,tau_u,grsig,gpsigma,F_mat,F_spat,gprep(u1:uf),elrhd)
       call sup_rhs_ssgs_S_NH(e,nd,tn,dvol,G,det,tau_s,gprep(s1:sf),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,gprep(u1:uf),elrhp)
       call sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,det,tau_p,gprep(p1),elrhp)

   end subroutine calculateRHS_SUP_NH_ASGS_repro

   subroutine calculateRHS_SUP_NH_OSGS_repro
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------ADJOINT TESTING OF PROJECTION -----------

       !!-------------------U RHS calculation----------
       !call sup_rhs_ssgs_U_NH_repro(e,nd,tn,dvol,tau_s,gprep(s1:sf),elrhu)
       !call sup_rhs_psgs_U_NH_repro(e,nd,tn,dvol,tau_p,gprep(p1),elrhu)
       !!-------------------S RHS calculation----------
       !call sup_rhs_usgs_S_NH_repro(e,nd,tn,dvol,G,det,tau_u,grsig,F_mat,F_spat,gpsigma,gprep(u1:uf),elrhd)
       !!call sup_rhs_ssgs_S_NH_repro(e,nd,tn,dvol,G,det,tau_s,gprep(s1:sf),elrhd)
       !!-------------------P RHS calculation----------
       !call sup_rhs_usgs_P_NH_repro(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,gprep(u1:uf),elrhp)
       !!call sup_rhs_psgs_P_NH_repro(e,nd,tn,dvol,G,lam,det,B,tau_p,gprep(p1),elrhp)

       !-------------------U RHS calculation----------
       call sup_rhs_ssgs_U_NH(e,nd,tn,dvol,tau_s,gprep(s1:sf),elrhu)
       call sup_rhs_psgs_U_NH(e,nd,tn,dvol,tau_p,gprep(p1),elrhu)
       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH(e,nd,tn,dvol,G,det,tau_u,grsig,gpsigma,F_mat,F_spat,gprep(u1:uf),elrhd)
       call sup_rhs_ssgs_S_NH(e,nd,tn,dvol,G,det,tau_s,gprep(s1:sf),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,gprep(u1:uf),elrhp)
       call sup_rhs_psgs_P_NH(e,nd,tn,dvol,lam,det,tau_p,gprep(p1),elrhp)

   end subroutine calculateRHS_SUP_NH_OSGS_repro

   subroutine calculateInGaussNonLinearTerms
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       call sld_calculateJacobianDer_Fspat_ki(e,nd,tn,det,F_spat,grdFmat,derJFki)

   end subroutine calculateInGaussNonLinearTerms

  subroutine PostGaussElmats_NH_SGS_nonlinear
     implicit none
     integer(ip) :: nd,tn

     call sup_NH%Mesh%GetNdime(nd)
     tn=(nd*(nd+1))/2

     !-------------------Tau_u matrix calculation----------
     call sup_elm_usgs_ES_NH_nonli(e,nd,tn,dvol,det,G,gpsigma,derJFki,divFmat,grdFspt,elm_usgs_ES)
     call sup_elm_usgs_EP_NH_nonli(e,nd,tn,dvol,det,G,gpsigma,derJFki,divFmat,grdFspt,elm_usgs_EP)
     call sup_elm_usgs_QS_NH_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,elm_usgs_QS)
     call sup_elm_usgs_QP_NH_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,elm_usgs_QP)

  end subroutine

   subroutine calculateRHS_SUP_NH_ASGS_nonli
       implicit none
       integer(ip) :: nd,tn

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       !-------------------S RHS calculation----------
       call sup_rhs_usgs_S_NH_nonli(e,nd,tn,dvol,det,G,gpsigma,derJFki,divFmat,grdFspt,tau_u,AU(u1:uf),elrhd)
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,tau_u,AU(u1:uf),elrhp)

       !-------------------nonlinear derivative part of Forces----------
       call sup_rhs_usgs_S_NH_nonli(e,nd,tn,dvol,det,G,gpsigma,derJFki,divFmat,grdFspt,tau_u,-elext,sup_NH%extForce(ielem)%a(s1:sf,:))
       !-------------------P RHS calculation----------
       call sup_rhs_usgs_P_NH_nonli(e,nd,tn,dvol,det,G,lam,gppress,derJFki,divFmat,grdFspt,tau_u,-elext,sup_NH%extForce(ielem)%a(p1,:))

   end subroutine calculateRHS_SUP_NH_ASGS_nonli

  subroutine calculateExternalForces
      implicit none

      !Compute contributions to RHS :
      elext = 0.0_rp
      call ProcPointer%ExternalForces   

  end subroutine

  subroutine integrateForces
      implicit none

      !Compute contributions to RHS :
      call ProcPointer%sld_elmrhu(e,dvol,elext,sup_NH%extForce(ielem)%a(u1:uf,:))

  end subroutine

  subroutine integrateForces_SGS
      implicit none
      integer(ip) :: nd,tn

      call sup_NH%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2

      call sup_rhs_usgs_S_NH_extforces(e,nd,tn,dvol,G,det,tau_u,F_mat,F_spat,grsig,gpsigma,sup_NH%extForce(ielem)%a(s1:sf,:),elext)
      call sup_rhs_usgs_P_NH_extforces(e,nd,tn,dvol,G,lam,det,tau_u,grdpre,gppress,F_mat,F_spat,sup_NH%extForce(ielem)%a(p1,:),elext)

  end subroutine

   subroutine ResetForceVector
         implicit none
   
      sup_NH%extForce(ielem)%a(:,:)=0.0_rp
   
   end subroutine

  subroutine PhysicalProp
      implicit none
      integer(ip) :: nd,tn

      call sup_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call sup_NH%SUPGetPhysicalParameters(nd,tn,densi,lam,G)
  end subroutine

  subroutine OnIeltyChange
     implicit none
      
     if (e%ielty /= ielty0) then
        call SetPointersElmopeSUP_NH
        ielty0 = e%ielty
     endif
  end subroutine

   subroutine lhsrhscheck
      implicit none
      integer(ip)          :: pn,nd,tn,u1,uf,s1,sf,p1,bc,df,idime,i,count,mn
      real(rp),allocatable :: var(:,:),aux(:,:),elLhs(:,:),realRhs(:,:),realForcesRhs(:,:)
      !Very useful formulation for ROM, instead of calculating residual 
      !we calculate the variable
      !by solving A*U_i+1 = r + A*U_i, ------> U_i+1 - U_i = d_u
      !Finally, the problem to solve is : U_i+1 = A^-1*(r + A*U_i)
      !Usual Newton Raphson assembly is: d_u = A^-1 r

      !We add the previous value so we do not calculate residual 
      !but the total variable

      call sup%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode
      df  = a%ndofn
      mn  = e%mnode

      call sup%Memor%alloc(df,pn,var,'var','sldsup_Elmope')
      call sup%Memor%alloc(df,pn,aux,'aux','sldsup_Elmope')
      call sup%Memor%alloc(df,pn,elLhs,'elLhs','sldsup_Elmope')
      call sup%Memor%alloc(df,pn,realRhs,'realRhs','sldsup_Elmope')
      call sup%Memor%alloc(df,pn,realForcesRhs,'realForcesRhs','sldsup_Elmope')

      do i=1,pn
          count = 1
          do idime=u1,uf
              var(idime,i) = eldisp(count,i,1)
              count = count + 1
          enddo
          count = 1
          do idime=s1,sf
              var(idime,i) = elsigma(count,i)
              count = count + 1
          enddo
          var(p1,i)    = elpress(1,i)
      enddo

      call matvecmult(df,pn,elmat,var,aux)

      !Compare with Totalview, ideally both elLhs and realRhs should be the same
      realForcesRhs(u1:uf,1:pn) = realForcesRhs(u1:uf,1:pn) + a%extForce(ielem)%a(u1:uf,:) 
      realForcesRhs(s1:sf,1:pn) = realForcesRhs(s1:sf,1:pn) + a%extForce(ielem)%a(s1:sf,:) 
      realForcesRhs(p1,1:pn)    = realForcesRhs(p1,1:pn)    + a%extForce(ielem)%a(p1,:)    

      realForcesRhs = realForcesRhs + elrhs
      elLhs         = -aux
      realRhs       =  elrhs

      call sup%Memor%dealloc(df,pn,var,'var','sldsup_Elmope')
      call sup%Memor%dealloc(df,pn,aux,'aux','sldsup_Elmope')
      call sup%Memor%dealloc(df,pn,elLhs,'elLhs','sldsup_Elmope')
      call sup%Memor%dealloc(df,pn,realRhs,'realRhs','sldsup_Elmope')
      call sup%Memor%dealloc(df,pn,realForcesRhs,'realForcesRhs','sldsup_Elmope')

   end subroutine lhsrhscheck

end module   

!Solids_ELMOPE subroutine   
subroutine sldsup_Elmope_NH(SldProblem)
   use Mod_sldsup_elmope_NH

   implicit none
   class(SUPSolidsProblem_NH), target :: SldProblem
   
   !Setting sld problem pointer in sld_BaseElmope
   up     => SldProblem
   sup_NH => SldProblem
   sup    => SldProblem
   a      => SldProblem

   call SetPointersGeneral
   call SetPointersPointForces
   call SetPointersElmopeSUP_NH
   
   !TODO: there is a memory fault somewhere, it jumps from the AllocateArrays to gpaccel
   call a%Memor%alloc(nd,  gpaccel,'gpaccel'  ,'sld_Elmope')
   call sld_elemLoop
   call a%Memor%dealloc(nd,  gpaccel,'gpaccel'  ,'sld_Elmope')

   !!----------------------------------------------
   !ielem = 1
   !call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_ElmLoop')

   !!Allocate Arrays in BaseElmope
   !call ProcHook%AllocateArrays

   !call ProcHook%PhysicalProp
   !call ProcHook%Initializations

   !call a%Mesh%GetNelem(nelem)

   !elements : do ielem = 1,nelem
   !    !Load Element
   !    a%kfl_foldedElements=.false.
   !    call a%Mesh%ElementLoad(ielem,e)  

   !    call ProcHook%OnIeltyChange

   !    !ElmatsToZero
   !    call ProcHook%ResetArrays

   !    call ProcHook%Gathers

   !    !if(sup_NH%itera==1) then
   !    !    eldisp  = rand()
   !    !    elsigma = rand()
   !    !    elpress = rand()
   !    !endif

   !    !Cartesian derivatives and Jacobian at center of gravity
   !    call e%elmdcg

   !    !Element length at center of gravity
   !    call e%elmlen

   !    dvolt0 = 0.0_rp

   !    call ProcHook%PhysicalProp

   !    !Build Jacobian, A(d) matrix for the latest state of the body
   !    gauss_points: do igaus=1,e%pgaus

   !        e%igaus = igaus

   !        call ProcHook%PreGauss

   !        call ProcHook%InGauss

   !        !Interpolate         
   !        call ProcHook%Interpolates

   !        call ProcHook%PostInterpolates

   !        call ProcHook%EndGauss

   !    enddo gauss_points

   !    call ProcPointer%PostGauss

   !    call ProcHook%PostGauss

   !    call ProcHook%preAssemblyLhs

   !    call ProcHook%AssemblyLhs

   !    call ProcHook%postAssemblyLhs

   !    call ProcHook%preAssemblyRhs

   !    call lhsrhscheck

   !    call ProcHook%AssemblyRhs

   !    call ProcHook%postAssemblyRhs

   !    call ProcHook%ToLinearSystem

   !    call lhsrhscheck


   !enddo elements

   !!Hook
   !call ProcHook%Finalizations
   !call ProcHook%DeallocateArrays

   !!DeallocateElement
   !call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_ElmLoop')

end subroutine sldsup_Elmope_NH
