module Mod_nsf_elmope_3rd
   use typre
   use Mod_NSFractionalStep
   use Mod_TimeIntegrator
   use Mod_Memor
   use Mod_Element
   use Mod_ConvectiveElement
   use Mod_php_SetTimeIntegrator
   use Mod_nsm_elmdir
   use Mod_php_elmdir

   !For setting pointers
   use Mod_nsm_BaseElmope
   use Mod_nsm_TemperatureCoupling
   use Mod_nsm_ComputeAdvectionVelocity
   use Mod_nsm_NonLinearDerivatives
   use Mod_nsm_PhysicalProperties
   use Mod_nsm_ExternalForces
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_TurbulenceModel
   use Mod_nsm_ComputeGpResidual
   use Mod_nsm_ComputeResidualProjection
   use Mod_nsm_InterpolateResidualProjection
   use Mod_nsm_SubgridSpaceResidual
   use Mod_nsm_ComputeTauSmoothing
   use Mod_nsm_ComputeTaus
   use Mod_nsm_ComputeTestf
   use Mod_nsm_ComputeSubscales
   use Mod_nsm_ComputeNonLinearSubscales
   use Mod_nsm_ComputeDissipation
   use Mod_nsf_DynSubsElmope
   use Mod_nsm_SwitchOff
   use Mod_nsm_LevelSetCoupling   
   use Mod_nsm_FreeSurfaceMatrices
   use Mod_nsm_HangingNodes
   use Mod_nsf_FreeSurfaceStabilization_3rd
   use Mod_nsm_ElasticBoundaryDir
   implicit none

   class(NSFractionalStepProblem), pointer :: b => NULL()
   
   real(rp), allocatable :: elpreInc(:,:)  !Pressure Increment
   real(rp)              :: gppreInc(1,1)
   real(rp)              :: val_aux
   integer(ip)           :: inode,ipoin
   
contains
   
   !--------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      
      !External Procedures
      procedure() :: NULLSUB

      integer(ip) :: kfl_nonlinear, nelty
      
      call ResetProcedureComposition
      
      !-------------------------------------------------------------------
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      
      !------------------------------------------------------------------
      !Initialize the ProcedureFlags
      call SetPointersPhysicalProperties%Initialize 
      call SetPointersSwitchOff%Initialize
      call SetPointersLevelSetCoupling%Initialize 
      call SetPointersFreeSurfaceMatrices%Initialize 
      call SetPointersHangingNodes%Initialize
      call SetPointersFSS_3rd%Initialize
      call SetPointersElasticBoundaryDir%Initialize

      !Now we set the pointers
      !------------------------------------------------------------------
      !Disconnection of elements
      if (a%kfl_SwitchOff == 1) call SetPointersSwitchOff%Set
         
      !HangingNodes
      call SetPointersHangingNodes%Set
         
      call SetPointersLevelSetCoupling%Set
      call SetPointersFreeSurfaceMatrices%Set 
      call SetPointersPhysicalProperties%SetTask('Elmope')
      call SetPointersPhysicalProperties%Set
      call SetPointersFSS_3rd%Set
      call SetPointersElasticBoundaryDir%Set

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange
      
      !Finalizations
      call SetPointersPhysicalProperties%Finalize
      call SetPointersSwitchOff%Finalize
      call SetPointersLevelSetCoupling%Finalize
      call SetPointersFreeSurfaceMatrices%Finalize 
      call SetPointersHangingNodes%Finalize
      call SetPointersFSS_3rd%Finalize
      call SetPointersElasticBoundaryDir%Finalize
   end subroutine

   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine

end module

subroutine nsf_elmope_3rd(NSProblem)
   !Navier-Stokes elemental operations (Fractional step - end-of-step u):
   !1. Compute elemental matrix and RHS 
   !2. Impose Dirichlet boundary conditions
   use typre
   use Mod_nsf_elmope_3rd
   implicit none

   class(NSFractionalStepProblem), target :: NSProblem
   
   interface
      subroutine nsm_elmrhf_eov(e,dvolu,denac,dtinv,gpveln,gppre,elrhs)
         use typre
         use Mod_Element
         implicit none
         class(FiniteElement)        :: e
         real(rp),    intent(in)    :: gpveln(e%ndime)
         real(rp),    intent(in)    :: dvolu,dtinv,gppre,denac
         real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      end subroutine   
   end interface
   
   a=>NSProblem%NavierStokesProblem
   b=>NSProblem

   !Properties (rho, nu, sigma.)
   acvis=a%MatProp(imat)%visco
   acden=a%MatProp(imat)%densi

   call SetPointers
   
   !Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   !If it is Crank Nicolson the third step is solved as BDF1
   if (a%kfl_tsche_1st_current == 'CN   ') then
      !The third step of the fractional step is solved as BDF1
      a%kfl_tsche_1st_current = 'BDF1 '
      call Integrator%Init(a%kfl_tsche_1st_current)
      call Integrator%GetLHSDtinv(a%dtinv,LHSdtinv)
      call Integrator%GetNumberOfTimeSteps(nsteps)
   endif

   !Memory Allocations
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsf_elmope_3rd')
   call a%Memor%alloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','nsf_elmope_3rd')
   
   call a%Memor%alloc(e%ndime,e%mnode,1,elvel,'elvel','nsf_elmope_3rd')
   call a%Memor%alloc(e%mnode,2,elpre,'elpre','nsf_elmope_3rd')
   call a%Memor%alloc(e%mnode,1,elpreInc,'elpreInc','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,1,gpvel,'gpvel','nsf_elmope_3rd')

   call ProcHook_Initializations

   !Loop over elements
   call a%Mesh%GetNelem(nelem)
   elements: do ielem=1,nelem
      call a%Mesh%ElementLoad(ielem,e)

      !Hook
      call ProcHook_OnIeltyChange      
      
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      wrmat1=0.0_rp

      !Gathering operations
      call e%gather(e%ndime,elvel(:,:,1),a%veloc(:,:,1)) ! int. u_n+1
      call e%gather(1,elpre(:,1),a%press(:,1)) ! p_n+1
      call e%gather(1,elpre(:,2),a%press(:,2)) ! p_n or previous pressure iteration
      
      !ElpreInc is the difference between pressure at n and pressure at n+1
      
      elpreInc(:,1) = elpre(:,1) - elpre(:,2)
      
      !Cartesian derivatives and Jacobian for linear elements
      call e%elmdcg
      
      !PreGauss
      call ProcHook_PreGauss

      !Loop on gauss points
      gauss_points: do igaus=1,e%pgaus
         e%igaus = igaus
         
         call e%elmder

         call ProcHook_InGauss
         
         dvol = e%weigp(e%igaus)*e%detjm

         ! Interpolation (Gauss point values)
         call e%interpg(e%ndime,elvel(:,:,1),gpvel(:,1)) ! u_int
         call e%interpg(1,elpreInc(:,1),gppreInc)

         call ProcHook_Interpolates
         
         !Physical Properties
         call a%GetPhysicalParameters(imat,acden,acvis)
         !Hook
         call ProcHook_PhysicalProp 

         !Advection velocity norm
         call vecnor(gpvel(:,1),e%ndime,gpvno,2)         

         call ProcHook_ComputeTaus

         call ProcHook_InGaussElmats
         
         !Compute contributions to RHS : Block U
         call nsm_elmrhf_eov(e,dvol,acden,LHSdtinv,gpvel(:,1),gppreInc(1,1),elrhs)
         
         !Compute contributions to elemental matrix : Block U,V
         val_aux = LHSdtinv*acden*dvol
         
         !if(a%kfl_colev==1 .and. a%kfl_fsurf==1)then
         !   call elmmas_Lumped(e,val_aux,wrmat1)
         !else
            call elmmas(e,val_aux,wrmat1)
         !end if

      end do gauss_points

      !Matrix composition
      do idime = 1,e%ndime
         elmat(idime,1:e%pnode,idime,1:e%pnode) = elmat(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      enddo
      
      !PreDirichlet
      call ProcHook_PreDirichlet

      !Prescribe Dirichlet boundary conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
      call php_elmdir(a,e,e%ndime,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
      
      do idime=1,e%ndime
         do inode=1,e%pnode
         if(abs(elmat(idime,inode,idime,inode)) <1e-14 .and. abs(elmat(idime,inode,idime,inode)) > 0.0_rp) then
            ipoin=e%lnods(inode)
            write(*,*) ielem,inode,ipoin,'3rd','pivot_zero'
         endif
         end do
      end do      
      
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   end do elements
   
   !Finalizations
   call ProcHook_Finalizations
   
   !Memory Deallocations
   call a%Memor%Dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsf_elmope_3rd')
   call a%Memor%Dealloc(e%mnode,e%mnode,wrmat1,'wrmat1','nsf_elmope_3rd')
   call a%Memor%Dealloc(e%ndime,e%mnode,elrhs,'elrhs','nsf_elmope_3rd')
   
   call a%Memor%Dealloc(e%ndime,e%mnode,1,elvel,'elvel','nsf_elmope_3rd')
   call a%Memor%Dealloc(e%mnode,2,elpre,'elpre','nsf_elmope_3rd')
   call a%Memor%Dealloc(e%mnode,1,elpreInc,'elpreInc','nsf_elmope_3rd')
   call a%Memor%Dealloc(e%ndime,1,gpvel,'gpvel','nsf_elmope_3rd')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsf_elmope_3rd')
end subroutine nsf_elmope_3rd



subroutine nsm_elmrhf_eov(e,dvolu,denac,dtinv,gpveln,gppre,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for Fractional step : end-of-step u
    !    (v, u_int/dt) + (p_n+1-p_n, div v)
    !
    !-----------------------------------------------------------------------
    use typre
    use Mod_Element
    implicit none
    class(FiniteElement)        :: e
    real(rp),    intent(in)    :: gpveln(e%ndime)
    real(rp),    intent(in)    :: dvolu,dtinv,gppre,denac
    real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)

    integer(ip)                :: inode,idime
    real(rp)                   :: aux

    do inode=1,e%pnode
       aux = e%shape(inode,e%igaus)*dvolu*denac*dtinv
       elrhs(:,inode) = gpveln(:)*aux + gppre*e%cartd(:,inode)*dvolu + elrhs(:,inode)
       
    end do
end subroutine nsm_elmrhf_eov

