module Mod_nsm_ComputeDissipation
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_ComputeSubscales
   use Mod_nsm_SubgridSpaceResidual
   use Mod_nsm_ComputeTestf
   use Mod_nsm_ComputeGpResidual
   use Mod_nsm_InterpolateResidualProjection
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_ComputeTaus
   implicit none
   private
   public SetPointersComputeDissipation
   
   type, extends(PointerSetter) :: SPComputeDissipation
contains
      procedure :: SpecificSet => SpecificSetComputeDissipation
   end type
   type(SPComputeDissipation) :: SetPointersComputeDissipation
   
   !Dissipation
   real(rp), allocatable :: eldis(:)
   real(rp), allocatable :: gpDisTest(:)
   real(rp) :: gpdis
   logical  :: logical_aux_nodal, logical_aux_gp
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SpecificSetComputeDissipation(d)
      implicit none
      class(SPComputeDissipation) :: d
      
      !First we decide if we need to compute the dissipations
      logical_aux_nodal = .false.
      logical_aux_gp = .false.
      !Dissipation as a nodal array
      if (a%npp_stepi(15) /= 0) then   
         if (mod(a%istep,a%npp_stepi(15))==0) then
            logical_aux_nodal = .true.
         endif
      endif
      if (a%kfl_dispa > 1 ) logical_aux_nodal = .true.
      !Dissipation as a gauss point array
      if (a%npp_stepi(20) /= 0) then   
         if (mod(a%istep,a%npp_stepi(20))==0) then
            logical_aux_gp = .true.
         endif
      endif
      
      !If nodal or gp dissipations need to be computed
      if ((logical_aux_gp .eqv. .true.) .or. (logical_aux_nodal .eqv. .true.)) then
      
         call ConcatenateProcedures(ProcHook_PreLoop,PreLoopDis)
         call ConcatenateProcedures(ProcHook_Initializations,InitializationsDis)
         call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroDis)
         
         !We need the gradients
         call SetPointersInterpolateGradients%Set
         
         !Put dissipation to zero
         call ConcatenateProcedures(ProcHook_Interpolates,GpDisToZero)
         
         !Smagorinsky
         if (a%kfl_cotur == -1) then
         
         !Numerical dissipation
         else
            !We will need to compute the Tau Values
            call SetPointersComputeTaus%Set
            
            !We will need to compute Testf
            call SetPointersComputeTestf%Set
         
            !-L*(u_h)
            if (a%kfl_repro == 1 .and. (a%kfl_repro_SkipFE == 0 .or. e%linea == 0) )then
               call runend('Dissipation: For OSS, only ready for SKIP FE and linear elements')
            elseif (a%kfl_repro > 1) then
               call runend('Dissipation: For OSS, only ready for SKIP FE and linear elements')
            elseif (a%kfl_repro == 1 .and. a%kfl_repro_SkipFE == 1 .and. e%linea == 1) then
               !We need the residual at the Gauss point
               call SetPointersComputeGpResidual%Set
               
               !We need to interpolate the residual projection at the gauss point
               call SetPointersInterpolateResidualProjection%Set
            
               !We will compute the dissipation as the product R_ort*Subscales (positive for static subscales : R_ort*Tau*R_ort)
               call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,ComputeGpDisTestOSSLinearSkipFE)
            elseif (a%kfl_repro == 0) then
               call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,ComputeGpDisTest)
            else
               call runend('Dissipation: option not ready')
            endif
            
            !Get the subscale at the Gauss Point
            call SetPointersGetSubscales%Set
         
            !Add numerical dissipation
            call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,ComputeGpDisNumerical)
            
            !Dynamic subgrid scales, additional dissipation only if ASGS 
            !corresponding to the term (u_h,|partial_t |tilde u)
            if (a%kfl_tacsg == 1 .and. a%kfl_repro == 0) then
               call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,AdditionalDissipation_DSS_ASGS)
            endif
            
         endif   
         
         !Add molecular and Smagorinsky dissipation
         call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,ComputeGpDis)
         
         !For Nodal point dissipation
         if (logical_aux_nodal .eqv. .true.) then
            call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,GpdisToEldis)
            call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyDis)
         endif
         
         !For Gauss point dissipation
         if (logical_aux_gp .eqv. .true.) then
            call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,GpdisToGPDissipation)
         endif
         
         call ConcatenateProcedures(ProcHook_Finalizations,FinalizationsDis)
         call ConcatenateProcedures(ProcHook_PostLoop,ProjectDissipation)
      endif  
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !-------------------------------------------------------------------
   !Dissipation
   subroutine PreLoopDis
      implicit none
      
      a%Dissipation = 0.0_rp
   end subroutine
   
   subroutine InitializationsDis
      implicit none
      
      !Matrices alloc
      call a%Memor%alloc(e%mnode,eldis,'eldis','nsm_EndsteElmope')
      call a%Memor%alloc(e%ndime+1,gpdisTest,'gpdisTest','nsm_EndsteElmope')
   end subroutine

   subroutine FinalizationsDis
      implicit none
      
      !Matrices dealloc
      call a%Memor%dealloc(e%mnode,eldis,'eldis','nsm_EnditeElmope')
      call a%Memor%dealloc(e%ndime+1,gpdisTest,'gpdisTest','nsm_EndsteElmope')
   end subroutine

   subroutine ElmatsToZeroDis
      implicit none
      
      eldis = 0.0_rp
   end subroutine
   
   subroutine GpDisToZero
      implicit none
      
      gpdis = 0.0_rp
   end subroutine

   !For Numerical Dissipation
   subroutine ComputeGpDisTest
      implicit none
      
      integer(ip) :: idime,jdime

      GpDisTest = 0.0_rp
      
      !Velocity subgrid scales
      !We use the test function here, which is -L*(v_h)*timom
      !As a consequence we need to use -testf/timom
      do idime = 1,e%ndime
         GpDisTest(idime) = -dot_product(testf(1:e%pnode),elvel(idime,1:e%pnode,1))/timom
      enddo
      
      !We also need to add the contribution of the pressure gradient to the GpDisTest
      GpDisTest(1:e%ndime) = GpDisTest(1:e%ndime) - grpre(1,1:e%ndime)
         
      !Pressure subgrid scales
      !div div
      GpDisTest(e%ndime+1) = divvel
   end subroutine
   
   subroutine ComputeGpDisTestOSSLinearSkipFE
      implicit none
      
      !Here we are cheating
      !We use the residual projection because we know it coincides with the projection of the test functions
      !For more general cases, which include the full residual, or non-linear elements
      !this won't work and needs to be specifically coded (an additional projection is required)
      GpDisTest(1:e%ndime) = -( gpres(1:e%ndime) - gprep(1:e%ndime) )
      GpDisTest(e%ndime+1) = gpres(e%ndime+1) - gprep(e%ndime+1)
   end subroutine
   
   subroutine ComputeGpDisNumerical
      implicit none
      
      !Numerical Diffusion
      gpdis = gpdis + dot_product(gpdisTest(1:e%ndime+1),GpSGS(1:e%ndime+1))
   end subroutine   
   
   subroutine ComputeGpDis
      implicit none
      real(rp) :: VGradNorm2
      
      !Molecular and turbulent (acvis smago) dissipation
      call dot(grvel,grvel,e%ndime*e%ndime,VGradNorm2)
      gpdis = gpdis + acvis*VGradNorm2
   end subroutine
   
   !Dissipation due to transient subscales if ASGS, (u_h,\partial \tilde u)
   subroutine AdditionalDissipation_DSS_ASGS
      implicit none
      
      !First order subscales
      gpdis = gpdis + dot_product(gpvel(1:e%ndime,1),(a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)-a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)))*ReferenceDtinv
   end subroutine
   
   subroutine GpdisToEldis
      implicit none
      
      eldis(1:e%pnode) = eldis(1:e%pnode) + e%shape(1:e%pnode,e%igaus)*gpdis*dvol
   end subroutine
   
   subroutine AssemblyDis
      implicit none
      
      !a%Dissipation(e%lnods(1:e%pnode)) = a%Dissipation(e%lnods(1:e%pnode)) + eldis(1:e%pnode)
      call a%Mesh%AssemblyToArray(e,1_ip,eldis,a%Dissipation) 
   end subroutine
   
   subroutine GpdisToGPDissipation
      implicit none
      
      !We store the gauss point dissipation in the array
      a%GPDissipation(ielem)%a(igaus) = gpdis
   end subroutine
   
   subroutine ProjectDissipation
      implicit none
      
      call a%Project(1,a%Dissipation) 
   end subroutine

end module
 
