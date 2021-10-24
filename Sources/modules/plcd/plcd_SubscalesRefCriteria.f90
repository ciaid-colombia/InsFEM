subroutine plcd_SubscalesRefCriteria(b,error,TotalEstimatedError)
   use typre
   use Mod_PLCD
   use Mod_Element
   use Mod_ZZErrorEstimator
   use Mod_plcd_BaseElmope
   use Mod_plcd_StrainGenerator
   use Mod_plcd_ExternalForces
   use MPI
   implicit none
   class(PLCDProblem), target :: b
   real(rp) :: error(:)
   
   real(rp) :: dispnorm,gpdisp(1),subscalesnorm,grdisp(9),graddispnorm,smoothStress(6)
   
   integer(ip) :: vsize,ndime,ierr,ielem2
   real(rp) ::  ortStress(6), ortStrain(6),symnorm
   
   real(rp), allocatable :: elSmoothGradient(:,:,:), ortGradDisp(:,:), GradDispSmooth(:,:), symortGradDisp(:,:)
   
   class(PLCDMaterial), pointer :: Material
   real(rp) :: G, ortStressNorm, ortDivu, tau2, taue, tau, TotalEstimatedError,TotalEstimatedError0
   integer(ip) :: idime,npoinLocal
   real(rp) :: weightfactor, ContinuityResidual,DisplacementSubscales(3)
   real(rp), allocatable :: elpre(:)
   real(rp) :: gppre(1), gradp(3),invk,gpMomentumResidual(3)
   real(rp), allocatable :: elMomentumResidual(:,:)
   real(rp) :: elext(4) = 0.0_rp
   
   a=> b


   if (.not. a%UseUPFormulation) call runend('Subscales cannot be used as error estimator if UPFormulation is not used')
   if (.not. a%UseSmoothedDisplacementGradient) call runend('Subscales cannot be used as error estimator if SmoothedDisplacementGradient is not used')
   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
   call a%Memor%alloc(e%ndime,e%mnode,1,eldisp,'eldisp','plcd_Elmope')
   call a%Memor%alloc(e%mnode,elpre,'elpre','plcd_Elmope')
   call a%Memor%alloc(e%ndime,e%mnode,elMomentumResidual,'elMomentumResidual','UPFormulation')
   
   call a%Mesh%GetNdime(ndime)
   call GetVoigtsize(ndime,vsize)

   call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elSmoothGradient,'elSmoothGradient','plcd_Elmope')
   call a%Memor%alloc(e%ndime,e%ndime,ortGradDisp,'ortGradDisp','plcd_Elmope')
   call a%Memor%alloc(e%ndime,e%ndime,symortGradDisp,'symortGradDisp','plcd_Elmope')
   call a%Memor%alloc(e%ndime,e%ndime,GradDisp,'GradDisp','plcd_Elmope')
   call a%Memor%alloc(e%ndime,e%ndime,GradDispSmooth,'GradDispSmooth','plcd_Elmope')         

   
   call a%Mesh%GetNelem(nelem)
   do ielem = 1,nelem
   
      ielem2 = ielem
      
      call a%Mesh%ElementLoad(ielem,e)  
      
      ElementMatData => a%ElementMaterialsData(ielem)%p
      
      call e%gather(e%ndime,eldisp(:,:,1),a%Displacement(:,:,1))
      call e%gather(1,elpre,a%Pressure(:,1))
      call e%gather(e%ndime,elMomentumResidual,a%UPResidualProjection)
      call e%elmdel
      call e%elmlen
      
      if (a%ErrorEstimatorTypeOfSubscales == 0) then
         !Orthogonal subscales
         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%elmder
            dvol = e%weigp(e%igaus)*e%detjm
            
            !Stabilization contribution 
            call e%gradient(1,elpre,gradP)
            call ElementMatData%GetSecantCNorm(e%igaus,G)
            call e%elmlen
            tau = a%staco(1)*e%hleng(1)**2/(2*G)
            
            !Split OSS, velocity subscales
            call e%interpg(e%ndime,elMomentumResidual,gpMomentumResidual)
            !gpMomentumResidual = 0.0_rp
            
            call GetExternalForces(e,a,elext)
            
            DisplacementSubscales(1:e%ndime) = tau*(gradP(1:e%ndime)-elext(1:e%ndime)-gpMomentumResidual(1:e%ndime))
            
            !call vecnor(a%UPSubscales(ielem)%a(1:e%ndime,igaus),e%ndime,subscalesnorm,2)
            
            call vecnor(DisplacementSubscales,e%ndime,subscalesnorm,2)
            
            call ElementMatData%GetSecantCNorm(e%igaus,G)
         
            !tau = 1*(e%hleng(1)**2)/(2*G)
            error(ielem) = error(ielem) + ((subscalesnorm**2)/tau)*dvol
            
         enddo
      elseif (a%ErrorEstimatorTypeOfSubscales == 1) then
         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%elmder
            dvol = e%weigp(e%igaus)*e%detjm
            
            !ASGS subscales as error indicator
            !Pressure gradient at gauss point
            call e%gradient(1,elpre,gradP)
            
            call ElementMatData%GetSecantCNorm(e%igaus,G)
            call e%elmlen
            tau = a%staco(1)*e%hleng(1)**2/(2*G)
            
            elext = 0.0_rp   
            call GetEXternalForces(e,a,elext)

            
            DisplacementSubscales(1:e%ndime) = tau*(gradP(1:e%ndime)-elext(1:e%ndime))
            call vecnor(DisplacementSubscales,e%ndime,subscalesnorm,2)
            error(ielem) = error(ielem) + ((subscalesnorm**2)/tau)*dvol
         enddo
      endif
         
      
      !If using smoothed stresses, the difference between smoothed and non-smoothed stresses is 
      !also an error indicator, which contributes as the "subscales on the element boundaries" term
      
      !Also we use the orthogonal projection of divu as pressure subscales
   
      call e%gather(e%ndime*e%ndime,elSmoothGradient,a%SmoothedDisplacementGradient)
   
      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         
         !ComputeDisplacementGradients
         call e%gradient(e%ndime,eldisp(:,:,1),gradDisp)
         !Real smooth computation
         call e%interpg(e%ndime*e%ndime,elSmoothGradient,GradDispSmooth)
         
         ortGradDisp = gradDisp-GradDispSmooth
         
         call ElementMatData%GetMaterialPointer(Material)
         call Material%CT%GetStrain(ortGradDisp,ortStrain)            
         call ElementMatData%GetConstitutiveTensorPointer(e%igaus,C)
         ortStress(1:vsize) =  matmul(C,ortStrain(1:vsize))
         
         call Material%CT%ComputeStressNorm(ortStress,ortStressNorm)
         
         call ElementMatData%GetSecantCNorm(e%igaus,G)
         
         !compute taue
         taue = (e%hleng(1)*0.5)/(2*G)
         error(ielem) = error(ielem) + 4*(ortStressNorm**2)*taue/(e%hleng(1))*dvol
         
         !Contribution from the pressure subscales
         !Orthogonal subscales
         if (a%ErrorEstimatorTypeOfSubscales == 0) then
            
            ortDivu = 0.0_rp
            do idime = 1,e%ndime
               ortDivu = ortDivu + ortGradDisp(idime,idime)
            enddo
            
            elext = 0.0_rp   
            call GetEXternalForces(e,a,elext)
            
            
            tau2 = 1*2*G
            error(ielem) = error(ielem) + (((ortDivu)**2)/tau2)*dvol
         
         !ASGS   
         elseif (a%ErrorEstimatorTypeOfSubscales == 1) then
            call e%interpg(1,elpre,gppre)
            
            call ElementMatData%GetInverseVolumetricDeformationModulus(e%igaus,invK)
         
            elext = 0.0_rp   
            call GetEXternalForces(e,a,elext)
         
            !Divergence term
            ContinuityResidual = 0.0_rp
            do idime = 1,e%ndime
               ContinuityResidual = ContinuityResidual + GradDisp(idime,idime)
            enddo
            ContinuityResidual = ContinuityResidual + invK*gppre(1) - elext(e%ndime+1)
            tau2 = 1*2*G
            error(ielem) = error(ielem) + ((ContinuityResidual**2)/tau2)*dvol
         endif   
      enddo
   enddo
   
   
   call a%Memor%dealloc(e%ndime,e%mnode,elMomentumResidual,'elMomentumResidual','UPFormulation')
   call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elSmoothGradient,'elSmoothGradient','plcd_Elmope')
   call a%Memor%dealloc(e%ndime,e%ndime,ortGradDisp,'ortGradDisp','plcd_Elmope')
   call a%Memor%dealloc(e%ndime,e%ndime,symortGradDisp,'symortGradDisp','plcd_Elmope')
   call a%Memor%dealloc(e%ndime,e%ndime,GradDisp,'GradDisp','plcd_Elmope')
   call a%Memor%dealloc(e%ndime,e%ndime,GradDispSmooth,'GradDispSmooth','plcd_Elmope')      
   
   call a%Memor%dealloc(e%ndime,e%mnode,1,eldisp,'eldisp','plcd_Elmope')
   call a%Memor%dealloc(e%mnode,elpre,'elpre','plcd_Elmope')
   
   call a%FilePostpr%postgp(error,'ErrorEstimator',a%istep,a%ctime,a%Mesh)
   
   !Compute the total error
   TotalEstimatedError0 = 0.0_rp
   call a%Mesh%GetNpoinLocal(npoinLocal)
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      weightfactor = 1.0-real(count(e%lnods(1:e%pnode)>npoinLocal))/real(e%pnode)
   
   
      TotalEstimatedError0 = TotalEstimatedError0 + error(ielem)*weightfactor
   enddo
   call MPI_ALLREDUCE(TotalEstimatedError0,TotalEstimatedError,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
   
   TotalEstimatedError = sqrt(TotalEstimatedError)
   if (a%MPIrank == a%MPIroot .and. a%kfl_adap) write(a%lun_adapt,*) ' Estimated error: ',TotalEstimatedError
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
   
end subroutine   
