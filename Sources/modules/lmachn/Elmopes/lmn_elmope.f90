module Mod_lmn_elmope
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_Element
   use Mod_LowMachElement
   use Mod_ConvectiveElement
   use Mod_LowMach
   use Mod_lmn_elmdir

   !For setting pointers
   use Mod_lmn_BaseElmope
   use Mod_lmn_ComputeAdvectionVelocity
   use Mod_lmn_ComputeTemperature
   use Mod_lmn_ComputeDensity
   use Mod_lmn_NonLinearDerivatives
   use Mod_lmn_ExternalForces
   use Mod_lmn_InterpolateGradients
   use Mod_lmn_ComputeGpResidual
   use Mod_lmn_ComputeResidualProjection
   use Mod_lmn_InterpolateResidualProjection
   use Mod_lmn_SubgridSpaceResidual
   use Mod_lmn_ComputeTaus
   use Mod_lmn_ComputeTestf
   use Mod_lmn_ComputeSubscales
   use Mod_lmn_ComputeNonLinearSubscales
   use Mod_lmn_DynSubsElmope
   use Mod_lmn_OrthogonalSubscales
   use Mod_lmn_Linearization
   use Mod_lmn_HangingNodes
   implicit none   
   
contains

   !----------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      use typre
      implicit none

      !External Procedures
      procedure() :: NULLSUB

      integer(ip) :: kfl_nonlinear,nelty
      
      call ResetProcedureComposition
      
      !-----------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
     
      !----------------------------------------------------------
      !Defaults

      !Pointers
      ProcPointer%PostGaussElmats => PostGaussElmats
      ProcPointer%lmn_elmbuv      => lmn_elmbuv
      ProcPointer%lmn_elmrhu      => lmn_elmrhu
      ProcPointer%lmn_elmrhp      => lmn_elmrhp
      ProcPointer%lmn_elmbuq      => lmn_elmbuq
      ProcPointer%lmn_elmbtw      => lmn_elmbtw
      ProcPointer%lmn_elmrht      => lmn_elmrht
      ProcPointer%TimeIntegrationToElext => lmn_TimeIntegrationToElext

      !Hooks
      !Allocate-Deallocate
      call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeArrays)
      call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)
      !ToZero
      call ConcatenateProcedures(ProcHook%ElmatsToZero,ToZeroArrays)
      !Assembly
      call ConcatenateProcedures(ProcHook%Assembly,AssemblyElmats)
      
      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersAdvectionVelocity(0)     
      call SetPointersTemperature(0)     
      call SetPointersDensity(0)     
      call SetPointersExternalForces(0)        
      call SetPointersInterpolateGradients(0)         
      call SetPointersComputeTaus(0)                  
      call SetPointersComputeTestf(0)                 
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersComputeSubgridSpaceResidual(0)  
      call SetPointersNonLinearDerivatives(0)         
      call SetPointersGetSubscales(0)   
      call SetPointersDynSubsElmope(0)
      call SetPointersOrthogonalSubscales(0)
      call SetPointersLinearization(0)
      call SetPointersHangingNodes(0)

      !-------------------------------------------------------------------
      !Now we set the required pointers 
      
      call SetPointersAdvectionVelocity(1)
      call SetPointersTemperature(1)
      call SetPointersDensity(1)

      !Confined flow
      if (a%kfl_confi == 1) call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussConfined)

      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call SetPointersNonLinearDerivatives(1)
         call ConcatenateProcedures(ProcHook%InGauss,InGaussVolumesNonLinear)
         call ConcatenateProcedures(ProcHook%InGaussElmats,ProcPointer%PostGaussElmats)
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsNonLinear)
         ProcPointer%PostGaussElmats => NULLSUB
      endif
      
      call SetPointersExternalForces(1)
      call SetPointersComputeTaus(1)
      call SetPointersComputeTestf(1)
      call SetPointersOrthogonalSubscales(1)
      
      !-----------------------------------------------------------
      !Tracking of the subscales
      if (a%kfl_tacsg == 1) call SetPointersDynSubsElmope(1)
    
      !-------------------------------------------------------
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

      !Linearization pointers
      call SetPointersLinearization(1)
      call SetPointersHangingNodes(1)

      !--------------------------------------------------------------------------
      !We deallocate the procedure flags, so that they can be set the next time
      call SetPointersAdvectionVelocity(100)     
      call SetPointersTemperature(100)     
      call SetPointersDensity(100)     
      call SetPointersExternalForces(100)        
      call SetPointersInterpolateGradients(100)         
      call SetPointersComputeTaus(100)                  
      call SetPointersComputeTestf(100)                 
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersComputeSubgridSpaceResidual(100)  
      call SetPointersNonLinearDerivatives(100)         
      call SetPointersGetSubscales(100)   
      call SetPointersOrthogonalSubscales(100)
      call SetPointersDynSubsElmope(100)
      call SetPointersLinearization(100)
      call SetPointersHangingNodes(100)
   end subroutine

   !--------------------------------------------------------
   !ToZero
   subroutine ToZeroArrays
      implicit none

      elmat=0.0_rp
      elrhs=0.0_rp
      elmuv=0.0_rp
      elrhu=0.0_rp
      wrmat1=0.0_rp
      wrmat2=0.0_rp
      elmpq=0.0_rp
      elmpv=0.0_rp
      elmuq=0.0_rp
      elrhp=0.0_rp
      elmtw=0.0_rp
      elmtv=0.0_rp
      elrht=0.0_rp
   end subroutine
 
   !----------------------------------------------------------
   !PostGauss Matrices
   subroutine PostGaussElmats
      implicit none
      
      !BLOCK U,V
      !Viscosity terms : mu*(grad v, gradt u)
      call elmvis(e,dvolt0,acvis,wrmat1)

      !Viscosity terms : mu*(gradt v, grad u)
      call lmn_elmvis_div(e,dvolt0,acvis,elmuv)
      call lmn_elmbuv_div(e,acvis,dvolt0,elmuv)

      !Compute contribution to elemental matrix : BLOCK P,Q
      call lmn_elmbpq(e,dvolt1,acden,elmpq)
         
      !Conduction term
      call elmvis(e,dvolt0,actco,wrmat2)
   end subroutine
   
   !----------------------------------------------------------
   !Assembly Matrices
   subroutine AssemblyElmats
      implicit none
      integer(ip) :: idime
      
      !Assembly wrmat1  in elmuv
      forall (idime = 1:e%ndime)
         elmuv(idime,1:e%pnode,idime,1:e%pnode) = elmuv(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      end forall   
  
      !Assembly wrmat2 in elmtw
      elmtw(1,1:e%pnode,1,1:e%pnode) = elmtw(1,1:e%pnode,1,1:e%pnode) + wrmat2(1:e%pnode,1:e%pnode)
  
      !Assembly elmuv to elmat
      elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) + elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      
      !Assembly elrhu to elrhs
      elrhs(1:e%ndime,1:e%pnode) = elrhu(1:e%ndime,1:e%pnode) + elrhs(1:e%ndime,1:e%pnode)
   
      !Assembly elmuq to elmat
      elmat(e%ndime+2,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(e%ndime+2,1:e%pnode,1:e%ndime,1:e%pnode) + elmuq(1,1:e%pnode,1:e%ndime,1:e%pnode)

      !Assembly elmpv & elmpq to elmat
      elmat(1:e%ndime,1:e%pnode,e%ndime+2,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,e%ndime+2,1:e%pnode) + elmpv(1:e%ndime,1:e%pnode,1,1:e%pnode)
      elmat(e%ndime+2,1:e%pnode,e%ndime+2,1:e%pnode) = elmat(e%ndime+2,1:e%pnode,e%ndime+2,1:e%pnode) + elmpq(1,1:e%pnode,1,1:e%pnode)

      !Assembly elmtv & elmtw to elmat
      elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,e%ndime+1,1:e%pnode) + elmtv(1:e%ndime,1:e%pnode,1,1:e%pnode)
      elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) = elmat(e%ndime+1,1:e%pnode,e%ndime+1,1:e%pnode) + elmtw(1,1:e%pnode,1,1:e%pnode)

      !Assembly elrhp to elrhs
      elrhs(e%ndime+2,1:e%pnode) = elrhp(1,1:e%pnode) + elrhs(e%ndime+2,1:e%pnode)
     
      !Assembly elrht to elrhs
      elrhs(e%ndime+1,1:e%pnode) = elrht(1,1:e%pnode) + elrhs(e%ndime+1,1:e%pnode)

   end subroutine
   
   !----------------------------------------------------------
   !NonLinear Elements   
   subroutine InGaussVolumesNonLinear
      implicit none
      dvolt0=0.0_rp
      dvolt1=0.0_rp
   end subroutine

   subroutine InGaussElmatsNonLinear
      implicit none
      call lmn_elmbuv_nln(e,dvol,acvis,testf_mom,elmuv)
      call lmn_elmbuq_nln(e,dvol,acvis,acden,timom,elmuq)
      call lmn_elmbtw_nln(e,dvol,actco,testf_ene,elmtw)
   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine

   !----------------------------------------------------------
   !Confined   
   subroutine InGaussConfined
      implicit none
      call lmn_elmbpq_cnf(e,a%epspe,gpden(1),acvis,dvol,elmpq)
      call lmn_elmrhp_cnf(e,a%epspe,gpden(1),acvis,dvol,gppre(1),elrhp)
   end subroutine

end module   
     
!LMN_ELMOPE subroutine   
subroutine lmn_elmope(LMProblem)
   use Mod_lmn_elmope
   implicit none
   class(LowMachProblem), target :: LMProblem
   a=>LMProblem
  
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lmn_elmope')
   
   !Allocate Arrays in BaseElmope
   call ProcHook%AllocateArrays

   
   !Physical Parameters
   a%dvolt = 0.0_rp
 
   !Initialize Statistics
   call a%InitStats
   
   !Hook
   call ProcHook%Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)  
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
      
      !ElmatsToZero
      call ProcHook%ElmatsToZero

      !Gathers
      call ElementGathers 
      !Hook
      call ProcHook%Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      !Compute the characteristic length scale
      call elmchl(e,a%kfl_advec,elvel,chale)
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         call ElementInterpolates
         !Hook
         call ProcHook%Interpolates
         
         !Temperature
         call ProcPointer%ComputeTemperature

         !Density
         call ProcPointer%ComputeDensity
         call vecnor(gpden,1,acden,2)

         !Advection velocity      
         call ProcPointer%ComputeAdvectionVelocity

         !Advection velocity norm
         call vecnor(gpadv,e%ndime,gpvno,2)
         
         !Compute a·grad(N) for v and T
         call ComputeAGradV(e,gpadv,AGradN)
         
         !Compute the stability parameters 
         call ProcHook%ComputeTaus
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         call ProcHook%ComputeTestf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                !w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                !w(gp)*detjm*tau1
         
         !Compute Elext, Temporal Derivatives
         elext_mom=0.0_rp
         elext_pre=0.0_rp
         elext_ene=0.0_rp

         !Compute vector of external forces
         call ProcPointer%ExternalForces
         
         !Time integration
         call ProcPointer%TimeIntegrationToElext
      
         !Gradient
         call e%gradient(1,eltem,grtem)
         call e%gradient(e%ndime,elvel,grvel)

         !InGaussElmats

         !RHS
         !Compute contributions to RHS : Block U
         call ProcPointer%lmn_elmrhu(e,Integrator,dvol,a%dtinv,ticon,gpden,gpvel,testf_mom,elext_mom,elrhu)

         !Compute contributions to RHS : Block T
         call ProcPointer%lmn_elmrht(e,Integrator,dvol,a%dtinv,testf_ene,elext_ene,actex,acden,gptem,gtemp,acpth,elrht)
            
         !Compute contributions to RHS : Block P
         call ProcPointer%lmn_elmrhp(e,Integrator,timom,dvol,a%dtinv,gpden,elext_mom,elrhp)

         !LHS
         !Compute contributions to elemental matrix : Block U,V
         call ProcPointer%lmn_elmbuv(e,dvol,acden,LHSdtinv,AGradN,testf_mom,elmuv)
         call lmn_elmbuv_con(e,acden,ticon,dvol,elmuv)
         
         !Compute contributions to elemental matrix : Block P,V
         call lmn_elmbpv(e,dvol,testf_mom,elmpv)
            
         !Compute contribution to elemental matrix : BLOCK T,V
         call lmn_elmbtv(e,acden,ticon,dvol,aGradN,actex,elmtv)

         !Compute contributions to elemental matrix : Block U,Q
         call ProcPointer%lmn_elmbuq(e,timom,dvol,acden,LHSdtinv,AGradN,elmuq)     
        
         !Compute contribution to elemental matrix : BLOCK T,W
         call ProcPointer%lmn_elmbtw(e,dvol,acden,LHSdtinv,AgradN,testf_ene,elmtw)

         !Hook
         call ProcHook%InGaussElmats
         
         !Statistics
         call a%InGaussStats(acden,acvis,gpvno,chale,timom)
         
      enddo gauss_points
            
      !Post Gauss Matrices    
      call ProcPointer%PostGaussElmats
      
      !Matrix composition
      call ProcHook%Assembly

      !Total volume
      a%dvolt = a%dvolt +dvolt0
   enddo elements
   
   call a%FinalizeStats
   
   !Hook
   call ProcHook%Finalizations
   
   !Arrays Deallocations
   call ProcHook%DeallocateArrays

   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','lmn_elmope')

end subroutine


