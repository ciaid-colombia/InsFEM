module Mod_plcd_UPFormulation
   use typre  
   use Mod_plcd_BaseElmope
   use Mod_PointerSetter
   use Mod_plcd_ExternalForces
   implicit none
   private
   public SetPointersUPFormulation, SetPointersUPFormulationEndite 
   
   type, extends(PointerSetter) :: SPUPFormulation
contains
      procedure :: SpecificSet => SpecificSetPointersUPFormulation
   end type 
   type(SPUPFormulation) :: SetPointersUPFormulation
   
   type, extends(PointerSetter) :: SPUPFormulationEndite
contains
      procedure :: SpecificSet => SpecificSetPointersUPFormulationEndite
   end type 
   type(SPUPFormulationEndite) :: SetPointersUPFormulationEndite
      
      
   real(rp), allocatable :: elpre(:)
   
   
   real(rp), allocatable :: elresp(:,:),UPResidualProjection(:,:),gradP(:), elMomentumResidual(:,:)
   real(rp) ::  gppre(1)   
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SpecificSetPointersUPFormulation(d)
      implicit none
      class(SPUPFormulation) :: d
            !-------------------------------------------------------
            !UPFormulation
            if (a%UseUpFormulation) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocElmope)
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocElmope)
               call ConcatenateProcedures(ProcHook%PreGauss,GathersElmope)
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyGMatrix)
            endif
   end subroutine   
   
   
   !For actual computations
   !------------------------------------------------------
   subroutine AllocElmope
      call a%Memor%alloc(e%ndime,e%mnode,elMomentumResidual,'elMomentumResidual','UPFormulation')
   
   
   end subroutine
   
   
   
   subroutine DeallocElmope
      call a%Memor%dealloc(e%ndime,e%mnode,elMomentumResidual,'elMomentumResidual','UPFormulation')
   
   
   
   end subroutine
   
   subroutine GathersElmope
      call e%gather(e%ndime,elMomentumResidual,a%UPResidualProjection)
   
   
   
   end subroutine
   
   
   
   subroutine ComputeAndAssemblyGMatrix
      real(rp) :: gelmat(e%ndime,e%pnode,1,e%pnode)
      
      integer(ip) :: idime,inode,jnode
      real(rp) :: invK,gpMomentumResidual(e%ndime)
      
      real(rp) ::   tau, G
      
      !(div v, p)
      forall(idime = 1:e%ndime, inode = 1:e%pnode, jnode = 1:e%pnode)
         gelmat(idime,inode,1,jnode) = e%cartd(idime,inode)*e%shape(jnode,e%igaus)*dvol
      end forall   
      !G mat
      GaussElmat(1:e%ndime,1:e%pnode,e%ndime+1:e%ndime+1,1:e%pnode) = GaussElmat(1:e%ndime,1:e%pnode,e%ndime+1:e%ndime+1,1:e%pnode) - gelmat
      
      !Gt mat
      forall(idime = 1:e%ndime, inode = 1:e%pnode, jnode = 1:e%pnode)
         GaussElmat(e%ndime+1,jnode,idime,inode) = GaussElmat(e%ndime+1,jnode,idime,inode) + gelmat(idime,inode,1,jnode)
      end forall
      
      
      !p p mat, without stabilization
      call ElementMatData%GetInverseVolumetricDeformationModulus(e%igaus,invK)
      
      forall (inode = 1:e%pnode,jnode = 1:e%pnode)
         GaussElmat(e%ndime+1,inode,e%ndime+1,jnode) = GaussElmat(e%ndime+1,inode,e%ndime+1,jnode) + invK*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol
      end forall
   
      !Stabilization contribution 
      call ElementMatData%GetSecantCNorm(e%igaus,G)
      call e%elmlen
      tau = a%staco(1)*e%hleng(1)**2/(2*G)
      
      !qp
      forall (inode = 1:e%pnode,jnode = 1:e%pnode)
         GaussElmat(e%ndime+1,inode,e%ndime+1,jnode) = GaussElmat(e%ndime+1,inode,e%ndime+1,jnode) + tau*dot_product(e%cartd(:,inode),e%cartd(:,jnode))*dvol
      end forall
       
      !Orthogonal subscales split oss
      !Nothing to be done here because a Newton Raphson scheme is used, it is already in the external forces
      
   end subroutine
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
   !ENDITE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   subroutine SpecificSetPointersUPFormulationEndite(d)
      implicit none
      class(SPUPFormulationEndite) :: d
            !-------------------------------------------------------
            !UPFormulation
            if (a%UseUpFormulation) then
               call ConcatenateProcedures(ProcHook%Initializations,Alloc)
               call ConcatenateProcedures(ProcHook%Finalizations,Dealloc)
               call ConcatenateProcedures(ProcHook%InGauss,PressureInterpolatePAndPassItToEMD)
               call ConcatenateProcedures(ProcHook%PreGauss,Gathers)
               call ConcatenateProcedures(ProcHook%InGaussElmats,PressureContributionToInternalForces)
            endif
   end subroutine   
   
   !For actual computations
   !------------------------------------------------------
   subroutine Alloc
      call a%Memor%alloc(e%mnode,elpre,'elpre','UPFormulation')
      call a%Memor%alloc(e%ndime,e%mnode,elMomentumResidual,'elMomentumResidual','UPFormulation')
      call a%Memor%alloc(e%ndime,gradP,'gradP','AllocOSS')
      
   end subroutine
   
   subroutine DeAlloc
      implicit none
      integer(ip) :: ipoin,npoin
      
      !Modify external forces with pressure boundary conditions
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         if (ipoin == a%nodpr) then
            if (a%UseUPFormulation .and. (a%kfl_confi == 1)) then
               a%ExternalForcesVector(a%ndofbc+1,ipoin) = a%InternalForcesVector(a%ndofbc+1,ipoin)
            endif
         endif
      enddo
   
      call a%Memor%dealloc(e%mnode,elpre,'elpre','UPFormulation')
      call a%Memor%dealloc(e%ndime,e%mnode,elMomentumResidual,'elMomentumResidual','UPFormulation')
      call a%Memor%dealloc(e%ndime,gradP,'gradP','AllocOSS')
   end subroutine
   
   subroutine Gathers
      call e%gather(1,elpre,a%Pressure(:,1))
      call e%gather(e%ndime,elMomentumResidual,a%UPResidualProjection)
   end subroutine
   
   subroutine PressureInterpolatePAndPassItToEMD
      call e%interpg(1,elpre,gppre)
      
      call ElementMatData%SetPressureForComputeHistoryAndConstitutiveTensor(gppre(1))
   end subroutine
   
   subroutine PressureContributionToInternalForces
      real(rp) :: gelmat(e%ndime,e%pnode,1,e%pnode), elext(a%ndofn)
      
      integer(ip) :: idime,inode,jnode
      real(rp) :: invK
      
      real(rp) :: divu,gpMomentumResidual(e%ndime),DisplacementSubscales(e%ndime)
      
      real(rp) ::  tau, G
      
      !Pressure gradient at gauss point
      call e%gradient(1,elpre,gradP)
      
      !Momentum equation G
      GaussElInternalForces(1:e%ndime,1:e%pnode) = GaussElInternalForces(1:e%ndime,1:e%pnode) - e%cartd(:,1:e%pnode)*gppre(1)*dvol
      
      !Continuity equation G
      divu = 0.0_rp
      do idime = 1,e%ndime
         divu = divu + gradDisp(idime,idime)
      enddo
      GaussElInternalForces(e%ndime+1,1:e%pnode) = GaussElInternalForces(e%ndime+1,1:e%pnode) + divu*e%shape(1:e%pnode,e%igaus)*dvol
      
      !p p mat, without stabilization
      call ElementMatData%GetInverseVolumetricDeformationModulus(e%igaus,invK)
      GaussElInternalForces(e%ndime+1,1:e%pnode) = GaussElInternalForces(e%ndime+1,1:e%pnode) + invK*e%shape(1:e%pnode,e%igaus)*gppre(1)*dvol
   
      !Stabilization contribution 
      call ElementMatData%GetSecantCNorm(e%igaus,G)
      call e%elmlen
      tau = a%staco(1)*e%hleng(1)**2/(2*G)
      
      !Split OSS, velocity subscales
      call e%interpg(e%ndime,elMomentumResidual,gpMomentumResidual)
      !gpMomentumResidual = 0.0_rp
      
      call GetExternalForces(e,a,elext)
      
      DisplacementSubscales = tau*(gradP-elext(1:e%ndime)-gpMomentumResidual)
      !qp stabilization
      forall (inode = 1:e%pnode)
         GaussElInternalForces(e%ndime+1,inode) = GaussElInternalForces(e%ndime+1,inode) + dot_product(e%cartd(:,inode),DisplacementSubscales)*dvol
      end forall
   
      !Store the subscales for postprocessing
      if (a%UPStoreSubscales) then
         a%UPSubscales(ielem)%a(1:e%ndime,e%igaus) = DisplacementSubscales
      endif
      
 
   end subroutine
   

   
end module
