module Mod_supm_elmope   
   use Mod_supm_BaseElmope
   use Mod_supm_TemperatureCoupling
   use Mod_supm_ExternalForces
   use Mod_supm_ComputeAdvectionVelocity
   use Mod_supm_PressTempSubscale
   use Mod_supm_PhysicalProperties
   use Mod_supm_TemporalDerivatives
   use Mod_supm_LevelSetCoupling
   use Mod_supm_StokesProblem
   use Mod_supm_EnrichElement
   use Mod_supm_DiscontCapturing
   use Mod_supm_PostGauss
   use Mod_supm_NonLinearElem
   use Mod_supm_OrthogonalSplitTerms
   use Mod_supm_ComputeGpSplitTerms
   use Mod_supm_InGaussNonLinear
   use Mod_supm_ElemVElastic
   use Mod_supm_ElemNoVElastic
   use Mod_supm_ElemLogarithmic
   use Mod_supm_LogarithmicProblem
   use Mod_supm_ComputeTaus
   use Mod_supm_ComputeTauSmoothing
   use Mod_supm_ComputeResidualProjection
   use Mod_supm_ComputeGpResidual
   use Mod_supm_InterpolateResidualProjection
   use Mod_supm_SubgridSpaceResidual
   use Mod_supm_DynSubsElmope
   use Mod_supm_OrthogonalSubscales
   use Mod_supm_InterpolateSplitTermsProjection
   use Mod_supm_ComputeSplitTermsProjection
   use Mod_supm_ComputeTestf
   implicit none
  
contains

   !*********************
   !SetPointers
   !*********************
   subroutine SetPointers
      use typre
      implicit none
      integer(ip) :: kfl_nonlinear,nelty
      
      call ResetProcedureComposition
      !Set All Pointers To NULL() (in Mod_supm_BaseElmope)
      call SetPointersAndHooksToNULL
      imat=1
      
      !Pointers
      ProcPointer%PreAssembly_sup => NULLSUB

      !----------------------------------------------------------
      !Initialize the ProcedureFlags
      !----------------------------------------------------------
      call SetPointersNoVElastic(0,'Elmope')
      call setPointersElemVElastic(0,'Elmope')
      call SetPointersDynSubsElmopeTF(0)
      call SetPointersInGaussNonLinearDerivatives(0)

      call SetPointersPostGauss(0)
      call SetPointersNonLinearElem(0,'Elmope')
      call SetPointersDiscontCapturing(0,'Elmope')
      call SetPointersEnrichElement(0)
      call SetPointersStokesProblem(0)
      call SetPointersLevelSetCoupling(0,'Elmope')
      call SetPointersComputeTaus(0)
      call SetPointersTemporalDerivatives(0)
      call SetPointersPhysicalPropertiesSUP(0,'Elmope')
      call SetPointersExternalForces(0)
      call SetPointersTemperatureCoupling(0)
      call SetPointersAdvectionVelocity(0)
      call SetPointersComputeGpResidual(0) 
      
      call SetPointersLogarithmicComponents(0)
      call SetPointersElemLCR(0,'elmope')
      
      call SetPointersComputeResidualProjection(0)
      call SetPointersInterpolateResidualProjection(0)
      call SetPointersOrthogonalSubscales(0)
      call SetPointersComputeSubgridSpaceResidual(0)
      
      
      call SetPointersComputeSplitTermsProjection(0)
      call SetPointersInterpolateSplitTermsProjection(0)
      call SetPointersOrthogonalSplitTherms(0)
      call SetPointersComputeTestf(0)
      call SetPointersPressTempSubscale(0)


      !----------------------------------------------------------
      !Set the required pointers
      !----------------------------------------------------------
      call SetPointersNoVElastic(1,'Elmope')
      call SetPointersPostGauss(1)
      call setPointersElemVElastic(1,'Elmope')
      call SetPointersComputeTestf(1)
      
      !Penalty
      if(a%kfl_penal ==1)  ProcPointer%PenaltyTerm => elmatpenalty
      
      auxdim=e%ndime
      !Dynamic subscales
      if (a%kfl_tacsg == 1) call SetPointersDynSubsElmopeTF(1)
      
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call SetPointersInGaussNonLinearDerivatives(1)
         call SetPointersNonLinearElem(1,'Elmope') 
      end if   
 
      if(kfl_nonlinear==1 .or. a%kfl_colev==1)then
          call ConcatenateProcedures(ProcHook_InGauss,InGaussVolumesNonLinear)
          call ConcatenateProcedures(ProcHook_InGaussElmats, ProcPointer%PostGaussElmats_sup)
          ProcPointer%PostGaussElmats_sup => NULLSUB      
      end if
      
      !ResidualProjection
      !call SetPointersComputeResidualProjection(1)
      call SetPointersOrthogonalSubscales(1)
      
      !Discontinuity Capturing
      call SetPointersDiscontCapturing(1,'Elmope')
      
      !Split-Oss
      !call SetPointersSplitOss(1,'Elmope')  
      call SetPointersOrthogonalSplitTherms(1)

      !Compute Taus
      call SetPointersComputeTaus(1)
      
      !PhysicalProperties
      call SetPointersPhysicalPropertiesSUP(1,'Elmope')
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange
      
      !Level Set
      call SetPointersLevelSetCoupling(1,'Elmope')  
      !Advection Velocity
      call SetPointersAdvectionVelocity(1)
      !Compute External forces
      call SetPointersExternalForces(1)
      !Temperature Coupling model
      call SetPointersTemperatureCoupling(1)
      !Logarithm conformation tensor model
      call SetPointersLogarithmicComponents(1)
      call SetPointersElemLCR(1,'elmope')
      call SetPointersPressTempSubscale(1)

      
      !----------------------------------------------------------
      !Deallocate the procedure flags, so that they can be set the next time
      !----------------------------------------------------------

      call SetPointersNoVElastic(100,'Elmope')
      call setPointersElemVElastic(100,'Elmope')
      call SetPointersDynSubsElmopeTF(100)
      call SetPointersInGaussNonLinearDerivatives(100)
      call SetPointersOrthogonalSplitTherms(100)
      call SetPointersNonLinearElem(100,'Elmope')
      call SetPointersPostGauss(100)
      call SetPointersDiscontCapturing(100,'Elmope')
      call SetPointersComputeResidualProjection(100)
      call SetPointersComputeSplitTermsProjection(100)
      call SetPointersComputeGpResidual(100)  
      call SetPointersInterpolateSplitTermsProjection(100)
      call SetPointersInterpolateResidualProjection(100)
      call SetPointersOrthogonalSubscales(100)
      call SetPointersComputeSubgridSpaceResidual(100) 
      call SetPointersEnrichElement(100)
      call SetPointersStokesProblem(100)
      call SetPointersLevelSetCoupling(100,'Elmope')
      call SetPointersComputeTaus(100)
      call SetPointersTemporalDerivatives(100)
      call SetPointersPhysicalPropertiesSUP(100,'Elmope')
      call SetPointersAdvectionVelocity(100)
      call SetPointersExternalForces(100)
      call SetPointersTemperatureCoupling(100)
      call SetPointersLogarithmicComponents(100)
      call SetPointersElemLCR(100,'elmope')
      call SetPointersComputeTestf(100)
      call SetPointersPressTempSubscale(100)


   end subroutine   
   
   !Nonlinear elements
   subroutine InGaussVolumesNonLinear
      implicit none
      call e%elmder
      call e%elmhes
      
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp
   end subroutine
   
   !Change of type of element
   subroutine OnIeltyChange
      implicit none
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
   !Penalty term
   subroutine elmatpenalty
      call supm_elmbpqpena(e,acvis,a%penal,dvol,elmpq) 
   end subroutine
   
   
   !*******************************************
   !General Subroutines called by SUPM_ELMOPE
   !*******************************************
   subroutine ComputeElmlen
      call e%elmlen
   end subroutine
   
   subroutine ComputeCGDerivatives
      call e%elmdcg
   end subroutine
   
   subroutine ComputeDerivatives
      call e%elmder
   end subroutine  

   subroutine ElextEstabToZero
      elextSEstab=0.0_rp
      elextSEstab2=0.0_rp
      elextSEstab3=0.0_rp
      elextSEstab4=0.0_rp
      elextSMat=0.0_rp
      elextSEstabMat=0.0_rp
      elextEstab=0.0_rp   
      elextEstab2=0.0_rp  
      elextEstab3=0.0_rp
      elextEstab4=0.0_rp
      elextEstab5=0.0_rp
   end subroutine
   
   subroutine PrintNodalValues 
      !This subroutine print in a Matlab File elmat (left hand side) 
      !and elrhs (right hand side) of a particular element.
      use Mod_iofile
      implicit none
      integer :: i,j,k,l, ipoin,cont, lunfile
      real    :: coordinates(e%ndime,e%pnode)
    
      if(ielem==3718) then
            call iofile(0,lunfile,'/home/lmoreno/Documentos/Matlab_files/LCR_element/elmsv.m','Restrict')
            
            cont=1
            do j=1, e%pnode
               do k=1, auxtens+e%ndime+1
                  write(lunfile,*) 'elmat(:,',cont,')=['!, j, k
                  do i=1, e%pnode
                     do l=1,auxtens+e%ndime+1
                        write(lunfile,*)  elmat(l,i,k,j)
                     end do  
                  end do
                  write(lunfile,*) '];'
               cont=cont+1
               end do
            end do

         write(lunfile,*) 'elrhs(:)=['
         do i=1,e%pnode
            do j=1,auxtens+e%ndime+1
            write(lunfile,*) elrhs(j,i)
            end do
         end do
         write(lunfile,*) '];'
      
               coordinates= e%elcod
            do ipoin=1,e%pnode
               write(lunfile,*) 'nodei(',ipoin,',:)=['
               do j=1,e%ndime
                  write(lunfile,*) coordinates(j,ipoin)
               end do
               write(lunfile,*) '];'
            end do
            
         coordinates= e%elcod
            do ipoin=1,e%pnode
               write(*,*) 'Nodo:', ipoin
               do j=1,e%ndime
                  write(*,*) coordinates(j,ipoin)
               end do
            end do
      end if
      call iofile(2,lunfile,'/home/lmoreno/Documentos/Matlab_files/LCR_element/elmsv.m','Restrict')
   end subroutine
   
end module   


!*****************************
!SUPM_ELMOPE subroutine  
!*****************************

subroutine supm_elmope(TFNSProblem)
   use Mod_supm_elmope
   implicit none
   class(ThreeFieldNSProblem) , target :: TFNSProblem
   real(rp) :: lamb
   interface
      subroutine supm_elmdir(a,e,elmat,elrhs)
         use typre
         use Mod_nsm_elmdir
         use Mod_ThreeField
         use Mod_Element
         use Mod_Mesh
         use Mod_php_elmdir
         implicit none
         class(ThreeFieldNSProblem)             :: a
         class(FiniteElement), intent(in)       :: e
         real(rp), intent(inout)                :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)
         integer(ip)                            :: aux,aux1
      end subroutine
   end interface
   
   a=>TFNSProblem
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)  

   ielem = 1
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sub_elem')

   !Set the Pointers for execution
   call SetPointers
   
   !Allocate matrices and variables
   call ProcHook_PreAllocate
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','sub_elem')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','sub_elem')
   call AllocateBaseElmopeArrays
   call AllocateBaseElmopeMatrices
   
   !Initialize Statistics
   call a%InitStats
   call ProcHook_Initializations

   call a%Mesh%GetNelem(nelem)
   
   do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)

      call ComputeCGDerivatives
      call ComputeElmlen
      
      call ProcHook_OnIeltyChange
   
      call ProcHook_PreGauss
            
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      elmuv=0.0_rp
      elrhu=0.0_rp
      wrmat1=0.0_rp
      elmpq=0.0_rp
      elmpv=0.0_rp
      elmuq=0.0_rp
      elrhp=0.0_rp
      !Constitutive Case
      wrmat2=0.0_rp
      elmst=0.0_rp
      elmut=0.0_rp
      elmpt=0.0_rp
      elmsv=0.0_rp
      elmsq=0.0_rp
      elrhc=0.0_rp
      
      elmsv2=0.0_rp
      elrhu2=0.0_rp

      call VelocityAndSigmaGathers
      
      call ProcHook_Gathers
      
      !Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale,a%kfl_hdifumin)
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      dvolt2=0.0_rp
      dvolt3=0.0_rp
      
      do igaus = 1,e%pgaus
         e%igaus = igaus
         call ProcHook_InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Velocity gradient 
         call e%gradient(e%ndime,elvel,grvel) 
         
         !Interpolate         
         call InterpolateGpVelocityAndSigma
         !Stress gradient 
         call e%gradient(auxtens,elsig,grsig)  
         !Interpolate hooks   
         call ProcHook_Interpolates
            
         !PostInterpolate
         call ProcHook_PostInterpolate
            
         !Physical Parameters
         call a%GetPhysicalParameters(imat,acden,acvis)        
         
         !Set physical options
         call ProcHook_PhysicalProp    

         !Advection velocity      
         call  ProcPointer%ComputeAdvectionVelocity_sup
            
         !Advection velocity norm
         call ProcPointer%ComputeAdvectionVelocityNorm_sup
            
         !Compute a·grad(V)
         call ProcPointer%ComputeAGradV_sup
            
         !Compute Taus
         call  ProcPointer%TauConstitutive 
   
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         !call nsm_ComputeTestf(e,acden,timom,AGradV,testf)
         call ProcHook_ComputeTestf
            
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*timom + dvolt1                ! w(gp)*detjm*tau1
         dvolt2=dvol*tidiv + dvolt2                ! w(gp)*detjm*tau2
            
         !Compute Elext, Temporal Derivatives
         elext = 0.0_rp
         elextC= 0.0_rp
         elextS= 0.0_rp
         call ElextEstabToZero
         
         !Compute vector of external forces
         call  ProcPointer%ExternalForces_sup          
            
         !Time integration
         call nsi_TimeIntegrationToEltemp(e,Integrator,acden,a%dtinv,gpvel,elext)
            
         !Viscoelastic time integration
         call ProcPointer%TemporalDerivatives   
            
         !InGaussElmats
         Call ProcPointer%TwoFieldTerms
            
         !Constitutive equation elemental components
         call ProcPointer%ConstitutiveComponents   
         
         !Laplacian components in viscoelastic fluid
         call ProcHook_HessianComponents    
         call ProcHook_InGaussElmats
            
         !Penalty
         call ProcPointer%PenaltyTerm

         !Statistics
         call a%InGaussStats_sup(acden,acvis,gpvno,chale,timom,tisig)
      enddo

      !PostGauss   
      call  ProcPointer%PostGaussElmats_sup 
      
      !Matrix composition
      ! Assembly wrmat1  in elmuv
      forall (idime = 1:e%ndime)
         elmuv(idime,1:e%pnode,idime,1:e%pnode) = elmuv(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
      end forall 
      !adition of (grad_sym(u),grad_sym(v)) => wrmat2 in elmuv
       elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) =elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) &
            + wrmat2(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)

      !Constitutive equation
      ! Assembly elmst to elmat
      elmat(1:auxtens,1:e%pnode,1:auxtens,1:e%pnode) = elmat(1:auxtens,1:e%pnode,1:auxtens,1:e%pnode) + elmst(1:auxtens,1:e%pnode,1:auxtens,1:e%pnode)
      ! Assembly elmut to elmat
      elmat(1:auxtens,1:e%pnode,auxtens+1:auxtens+e%ndime,1:e%pnode) = elmat(1:auxtens,1:e%pnode,auxtens+1:auxtens+e%ndime,1:e%pnode) + &
            elmut(1:auxtens,1:e%pnode,1:e%ndime,1:e%pnode)
      ! Assembly elmpt to elmat
      elmat(1:auxtens,1:e%pnode,auxtens+e%ndime+1,1:e%pnode) = elmat(1:auxtens,1:e%pnode,auxtens+e%ndime+1,1:e%pnode) + elmpt(1:auxtens,1:e%pnode,1,1:e%pnode)
      ! Assembly elmrhc to elrhs
      elrhs(1:auxtens,1:e%pnode) = elrhc(1:auxtens,1:e%pnode) + elrhs(1:auxtens,1:e%pnode)    
   
      ! Momentum equation
      ! Assembly elmsv to elmat
      elmat(auxtens+1:auxtens+e%ndime,1:e%pnode,1:auxtens,1:e%pnode) = elmat(auxtens+1:auxtens+e%ndime,1:e%pnode,1:auxtens,1:e%pnode)  &
         + elmsv(1:e%ndime,1:e%pnode,1:auxtens,1:e%pnode)
      ! Assembly elmuv to elmat
      elmat(auxtens+1:auxtens+e%ndime,1:e%pnode,auxtens+1:auxtens+e%ndime,1:e%pnode) = elmat(auxtens+1:auxtens+e%ndime,1:e%pnode,auxtens+1:auxtens+e%ndime,1:e%pnode) &
         + elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      ! Assembly elmpv to elmat
      elmat(auxtens+1:auxtens+e%ndime,1:e%pnode,auxtens+e%ndime+1,1:e%pnode) = elmat(auxtens+1:auxtens+e%ndime,1:e%pnode,auxtens+e%ndime+1,1:e%pnode) &
         + elmpv(1:e%ndime,1:e%pnode,1,1:e%pnode)
      ! Assembly elrhu to elrhs
      elrhs(auxtens+1:auxtens+e%ndime,1:e%pnode) = elrhs(auxtens+1:auxtens+e%ndime,1:e%pnode) + elrhu(1:e%ndime,1:e%pnode)

      ! Continuity Equation
      !Assembly elmsq to elmat
      elmat(auxtens+1+e%ndime,1:e%pnode,1:auxtens,1:e%pnode) = elmat(auxtens+1+e%ndime,1:e%pnode,1:auxtens,1:e%pnode) + elmsq(1,1:e%pnode,1:auxtens,1:e%pnode)   
      ! Assembly elmuq to elmat
      elmat(auxtens+e%ndime+1,1:e%pnode,auxtens+1:auxtens+e%ndime,1:e%pnode) = elmat(auxtens+1+e%ndime,1:e%pnode,auxtens+1:auxtens+e%ndime,1:e%pnode) &
         + elmuq(1,1:e%pnode,1:e%ndime,1:e%pnode)
      ! Assembly elmpq to elmat
      elmat(auxtens+e%ndime+1,1:e%pnode,auxtens+e%ndime+1,1:e%pnode) = elmat(auxtens+e%ndime+1,1:e%pnode,auxtens+e%ndime+1,1:e%pnode) + elmpq(1,1:e%pnode,1,1:e%pnode)
      ! Assembly  elrhp to elrhs
      elrhs(auxtens+e%ndime+1,1:e%pnode) = elrhs(auxtens+e%ndime+1,1:e%pnode) + elrhp(1,1:e%pnode) 
      
      !call PrintNodalValues
      
      call  ProcPointer%PreAssembly_sup       
      
      !Dirichlet Boundary Conditions
      call supm_elmdir(a,e,elmat,elrhs)
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo      

   !Finalizations
   call a%FinalizeStats
   call ProcHook_Finalizations 
   
   !DeallocateMatrices
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','sub_elem')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','sub_elem')
   call DeallocateBaseElmopeMatrices
   call DeallocateBaseElmopeArrays
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','sub_elem')
   
end subroutine

