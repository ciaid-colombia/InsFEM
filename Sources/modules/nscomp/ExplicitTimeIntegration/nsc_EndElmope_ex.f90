module Mod_nsc_EndElmope_ex
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_Mesh
   use Mod_NSCompressibleElement
   use Mod_NSCompressibleExplicit
   use Mod_NSCompressibleExplicitElement
   use Mod_NSCompressibleSubroutines
   use Mod_ConvectiveElement  
   use Mod_Element
   use Mod_NscExacso    

   implicit none
   class(NSCompressibleExplicitProblem), pointer :: a => NULL()
   type(NscExacso) :: exacso     
   
   !Pointers
   type :: PPointer
      procedure(), NOPASS, pointer :: ExternalForces => NULL()
      procedure(), NOPASS, pointer :: ComputeConvectiveRes => NULL()
      procedure(), NOPASS, pointer :: ComputeUnsteadyRes => NULL()
      procedure(), NOPASS, pointer :: ComputeDiffusiveRes => NULL()
      procedure(), NOPASS, pointer :: ComputeForcesRes => NULL()
      procedure(), NOPASS, pointer :: ComputeDCUnsteadyRes => NULL()
      procedure(), NOPASS, pointer :: ComputeInPoint => NULL()
   end type
   type(PPointer) :: ProcPointer

   !Hooks
   type :: PHook
      procedure(), NOPASS, pointer :: Initializations => NULL()
      procedure(), NOPASS, pointer :: Gathers => NULL()
      procedure(), NOPASS, pointer :: OnIeltyChange => NULL()
      procedure(), NOPASS, pointer :: ElmatsToZero => NULL()
      procedure(), NOPASS, pointer :: PreGauss => NULL()
      procedure(), NOPASS, pointer :: InGauss => NULL()
      procedure(), NOPASS, pointer :: Interpolates => NULL()
      procedure(), NOPASS, pointer :: PhysicalProp => NULL()         
      procedure(), NOPASS, pointer :: ComputeVariables => NULL()
      procedure(), NOPASS, pointer :: ComputeTaus => NULL()
      procedure(), NOPASS, pointer :: InGaussElmats => NULL()
      procedure(), NOPASS, pointer :: InGaussElmatsAssembly => NULL()
      procedure(), NOPASS, pointer :: Assembly => NULL()
      procedure(), NOPASS, pointer :: Finalizations => NULL()
      procedure(), NOPASS, pointer :: PostLoop => NULL()         
    end type
   type(PHook) :: ProcHook
   
   class(FiniteElement) , pointer     :: e => NULL()
   integer(ip) :: ielem,nelem

   real(rp), pointer :: vmass(:) => NULL()
  
   real(rp), allocatable :: wrhg(:,:,:),elrhg(:,:,:)
   real(rp)              :: dvol

   integer(ip)           :: ielty0 = 0 !Previous element type

   real(rp), allocatable :: elden(:,:)
   real(rp), allocatable :: elmom(:,:,:)
   real(rp), allocatable :: elene(:,:)
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: elext(:)
   real(rp), allocatable :: elexd(:),elexm(:),elexe(:)
   real(rp), allocatable :: elrgm(:,:,:), elrge(:,:)
   real(rp), allocatable :: AGradV(:), hessV(:,:,:), laplV(:) 
   real(rp), allocatable :: gpden(:), gpmom(:,:), gpene(:)
   real(rp), allocatable :: gpomo(:), gpove(:)
   real(rp), allocatable :: novst(:,:),elvst(:,:)  
   real(rp), allocatable :: grden(:), grmom(:,:), grene(:)
   real(rp), allocatable :: elgep(:,:), elgmp(:,:,:)
   real(rp), allocatable :: geprj(:), gmprj(:,:)
   real(rp), allocatable :: orthe(:), orthm(:,:)
   real(rp), allocatable :: dcrem(:), dcree(:)
   real(rp), allocatable :: visten(:,:), cndten(:,:)
   real(rp), allocatable :: heden(:,:,:),hemom(:,:,:),heene(:,:,:)
   real(rp), allocatable :: elred(:), elrem(:), elree(:)
   real(rp), allocatable :: elcd(:), elcm(:), elce(:)
   real(rp), allocatable :: vgmom(:),eldvs(:)
   real(rp), allocatable :: dvist(:) 
   real(rp), allocatable :: wrepro(:,:),elres(:,:) 
   real(rp), allocatable :: elrep(:,:),reprj(:) 
   real(rp), allocatable :: gpR(:)
   real(rp), allocatable :: Subscales(:) 

   real(rp)    :: elexh
   real(rp)    :: chale(2),timom(3)
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: arvis,artco
   real(rp)    :: gppre,gptem,gpspd
   real(rp)    :: gpode,gpoen
   real(rp)    :: gpmno,divmom
   real(rp)    :: vgden,vgene
   real(rp)    :: acgamma,actcn,invgpd,invcvh
   real(rp)    :: aux,sqinvgpd,sqgpmn,sqgpvn,aux_t,aux_d
   real(rp)    :: dhflx
   integer(ip) :: igaus,idime,npoin,ipoin
   
   integer(ip) :: calculate_primitives, compute_statistics
   !ALE
   logical :: isALE
   real(rp), allocatable :: elmve(:,:),gpmve(:)
   real(rp), pointer :: meshve(:,:,:) => NULL()
   real(rp), allocatable :: ALEGradV(:), mvgmom(:)
   real(rp)    :: mvgden,mvgene
   
   character(6) :: itask
   

   !#$COMPOSEPROCS 50
#include "COMPOSEPROCS_POINTERS_50.i90"   
   
contains

#include "COMPOSEPROCS_SUBROUTINES_50.i90"
    
   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointers
      implicit none
      
      !External Procedures
      procedure() :: NULLSUB

      integer(ip) :: kfl_nonlinear, nelty
      
      !Reset Procedure Composition
      call ResetProcedureComposition

      !Pointers
      ProcPointer%ExternalForces  => ComputeForces
      ProcPointer%ComputeConvectiveRes => ComputeConvectiveRes
      ProcPointer%ComputeUnsteadyRes => ComputeUnsteadyRes
      ProcPointer%ComputeDiffusiveRes => NULLSUB
      ProcPointer%ComputeForcesRes => NULLSUB
      ProcPointer%ComputeDCUnsteadyRes => NULLSUB
      ProcPointer%ComputeInPoint => NULLSUB

      !Hooks
      ProcHook%Initializations => NULLSUB
      ProcHook%Gathers         => NULLSUB
      ProcHook%OnIeltyChange   => NULLSUB
      ProcHook%ElmatsToZero   => NULLSUB
      ProcHook%PreGauss        => NULLSUB
      ProcHook%InGauss         => NULLSUB
      ProcHook%Interpolates    => NULLSUB
      ProcHook%PhysicalProp    => NULLSUB    
      ProcHook%ComputeVariables => NULLSUB
      ProcHook%ComputeTaus     => NULLSUB 
      ProcHook%InGaussElmats   => NULLSUB
      ProcHook%InGaussElmatsAssembly   => NULLSUB
      ProcHook%Assembly => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      ProcHook%PostLoop        => NULLSUB    
      
      call ConcatenateProcedures(ProcHook%Initializations,AllocateBase)
      call ConcatenateProcedures(ProcHook%Gathers,GatherBase)
      call ConcatenateProcedures(ProcHook%InGauss,CalculateElmchl)
      call ConcatenateProcedures(ProcHook%Interpolates,InterpolateBase)
      call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeLinearVariables)
      call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeResidual)
      call ConcatenateProcedures(ProcHook%Finalizations,DeallocateBase)

      !-----------------------------------------------------------
      !ExactSol
      if (a%kfl_exacs/=0) then  
         call ConcatenateProcedures(ProcHook%Initializations,AllocateExactSol)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocateExactSol)
         ProcPointer%ExternalForces => ExactSolutionForces
         call ConcatenateProcedures(ProcPointer%ComputeForcesRes,ComputeExactSolForcesRes)
      endif    

      !-----------------------------------------------------------
      !Non-linear subscales
      if (a%kfl_nolsg == 1) then
         call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeNonlinearVariables)
      endif
   
      call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeAuxiliarVariables)

      !State law 
      if (a%lawde /= 0) then
         if (a%lawde == 1) then
            call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeIdealAuxVariables)
            call ConcatenateProcedures(ProcPointer%ComputeConvectiveRes,ComputeConvectiveRes_i)
         else if (a%lawde /= 1) then
            call runend('Nsc_elmope: Non-ideal state law not ready')
         endif
      endif

      call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeSoundSpeed)
      call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeBaseGradients)

      call ConcatenateProcedures(ProcHook%ComputeTaus,ComputeTaus)

      !-----------------------------------------------------------
      !Non-linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%Initializations,AllocateHighOrderDerivatives)
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
         call ConcatenateProcedures(ProcPointer%ComputeDiffusiveRes,ComputeDiffusiveRes)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocateHighOrderDerivatives)
         !State law
         if (a%lawde /= 0) then
            if (a%lawde == 1) then
               call ConcatenateProcedures(ProcPointer%ComputeDiffusiveRes,ComputeDiffusiveRes_i)
            else if (a%lawde /= 1) then
               call runend('Nsc_elmope: Non-ideal state law not ready')
            endif
         endif
      endif
      call a%Mesh%GetALE(isALE)
      if (isALE) then
         call a%Mesh%GetMeshVeloc(meshve)
         call ConcatenateProcedures(ProcHook%Initializations,AllocMeshVeloc)
         call ConcatenateProcedures(ProcHook%Gathers,GatherMeshVeloc)
         call ConcatenateProcedures(ProcHook%Interpolates,InterpolateMeshVeloc)
         call ConcatenateProcedures(ProcHook%ComputeVariables,ComputeALEGradients)
         call ConcatenateProcedures(ProcPointer%ComputeConvectiveRes,ALEConvectiveContribution)
         call ConcatenateProcedures(ProcHook%Finalizations,DeallocMeshVeloc)
      end if

      !---------------------------------------------------------------
      if (itask .eq. 'Endpro' ) then 
              call a%Mesh%GetNpoin(npoin)
              call ConcatenateProcedures(ProcHook%PostLoop,PointsLoop)
              !Residual projection
              if (a%kfl_repro == 1) then 
                 ProcPointer%ComputeUnsteadyRes => NULLSUB
                 call ConcatenateProcedures(ProcHook%Initializations,AllocateRep)
                 call ConcatenateProcedures(ProcHook%ElmatsToZero,ElmatsToZeroRep)
                 call ConcatenateProcedures(ProcHook%InGaussElmatsAssembly,GpresToElres)
                 call ConcatenateProcedures(ProcHook%Assembly,AssemblyResidual)
                 call ConcatenateProcedures(ProcPointer%ComputeInPoint,ResidualProjSol)
                 call ConcatenateProcedures(ProcHook%PostLoop,DeallocateRep)
              end if
              !Gradient Projection calculation only for output
              if (a%kfl_shock == 2) then
                 call a%Mesh%GetNpoin(npoin)
                 call ConcatenateProcedures(ProcHook%Initializations,AllocateGradProj)
                 call ConcatenateProcedures(ProcHook%ElmatsToZero,ZeroGradProj)
                 call ConcatenateProcedures(ProcHook%InGaussElmats,CalculateGradProj)
                 call ConcatenateProcedures(ProcHook%Assembly,AssemblyGradProj)
                 call ConcatenateProcedures(ProcPointer%ComputeInPoint,GradientProjSol)
                 call ConcatenateProcedures(ProcHook%PostLoop,DeallocateGradProj)
              endif
      end if

      !---------------------------------------------------------------
      if (itask .eq. 'Endste' ) then 
              call ConcatenateProcedures(ProcHook%InGaussElmats,TotalResidual)
              !Orthogonal subscales
              if (a%kfl_repro == 1) then 
                 ProcPointer%ComputeUnsteadyRes => NULLSUB
                 call ConcatenateProcedures(ProcHook%Initializations,AllocateOSS)
                 call ConcatenateProcedures(ProcHook%Finalizations,DeallocateOSS)
                 call ConcatenateProcedures(ProcHook%Gathers,GatherOSS)
                 call ConcatenateProcedures(ProcHook%Interpolates,InterpolateOSS)
                 call ConcatenateProcedures(ProcHook%InGaussElmats,ResidualOSS)
                 if (a%kfl_shock == 1) then
                    call ConcatenateProcedures(ProcPointer%ComputeDCUnsteadyRes,ComputeDCUnsteadyRes)
                 endif
              endif 

              if (a%kfl_trasg == 1 .and. a%kfl_tacsg == 0) then 
                 call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeStaticSubscales)
              endif 

              !-----------------------------------------------------------
              !ShockCapturing added diffusion output
              if (a%npp_stepi(9)>0) then
                 if (a%kfl_shock == 1) then
                    call ConcatenateProcedures(ProcHook%Initializations,AllocateResDC)
                    call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeDCResidual)
                    call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeResDiff)
                    call ConcatenateProcedures(ProcHook%Finalizations,DeallocateResDC)
                    if (a%kfl_repro == 1) then
                       call ConcatenateProcedures(ProcPointer%ComputeDCUnsteadyRes,ComputeDCUnsteadyRes)
                    endif
                 else if (a%kfl_shock == 2) then
                    call ConcatenateProcedures(ProcHook%Initializations,AllocateGradProjDC)
                    call ConcatenateProcedures(ProcHook%Gathers,GatherGradProjDC)
                    call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGradProjDC)
                    call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGradDiff)
                    call ConcatenateProcedures(ProcHook%Finalizations,DeallocateGradProjDC)
                 else if (a%kfl_shock > 2) then
                    call runend('Nsc_EndElmope: Other shock capturing methods not ready')
                 endif
              call ConcatenateProcedures(ProcHook%InGaussElmats,OutputSCDiffusion)
              endif

              !-----------------------------------------------------------
              !Subgrid added diffusion output
              if (a%npp_stepi(10)>0) then
                 call ConcatenateProcedures(ProcHook%ComputeTaus,OutputSGDiffusion)
              endif
              !-----------------------------------------------------------
              !Gauss point residual output
              if (a%npp_stepi(13)>0) then
                 call ConcatenateProcedures(ProcHook%InGaussElmats,OutputResidual)
              endif

      end if

      !---------------------------------------------------------------
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange

      
   end subroutine

   !----------------------------------------------------------
   subroutine AllocateBase
      
      call a%Memor%alloc(1,elred,'elred','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,elrem,'elrem','nsc_EndElmope_ex')
      call a%Memor%alloc(1,elree,'elree','nsc_EndElmope_ex')
      call a%Memor%alloc(1,elcd,'elcd','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,elcm,'elcm','nsc_EndElmope_ex')
      call a%Memor%alloc(1,elce,'elce','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','nsc_EndElmope_ex')
      call a%Memor%alloc(e%mnode,a%ncomp-1,elden,'elden','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%mnode,a%ncomp-1,elmom,'elmom','nsc_EndElmope_ex')
      call a%Memor%alloc(e%mnode,a%ncomp-1,elene,'elene','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,elext,'elext','nsc_EndElmope_ex')
      call a%Memor%alloc(a%ncomp-1,gpden,'gpden','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,a%ncomp-1,gpmom,'gpmom','nsc_EndElmope_ex')
      call a%Memor%alloc(a%ncomp-1,gpene,'gpene','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,grden,'grden','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,grmom,'grmom','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,grene,'grene','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,gpomo,'gpomo','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,gpove,'gpove','nsc_EndElmope_ex')
      call a%Memor%alloc(e%mnode,AGradV,'AGradV','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,elvst,'elvst','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,novst,'novst','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,visten,'visten','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,cndten,'cndten','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,vgmom,'vgmom','nsc_EndElmope_ex')
      call a%Memor%alloc(a%ndofn,gpR,'gpR','nsc_EndElmope_ex')
      call a%Memor%alloc(a%ndofn,Subscales,'Subscales','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine DeallocateBase
      
      call a%Memor%dealloc(1,elred,'elred','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,elrem,'elrem','nsc_EndElmope_ex')
      call a%Memor%dealloc(1,elree,'elree','nsc_EndElmope_ex')
      call a%Memor%dealloc(1,elcd,'elcd','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,elcm,'elcm','nsc_EndElmope_ex')
      call a%Memor%dealloc(1,elce,'elce','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%mnode,a%ncomp-1,elden,'elden','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%mnode,a%ncomp-1,elmom,'elmom','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%mnode,a%ncomp-1,elene,'elene','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,elext,'elext','nsc_EndElmope_ex')
      call a%Memor%dealloc(a%ncomp-1,gpden,'gpden','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,a%ncomp-1,gpmom,'gpmom','nsc_EndElmope_ex')
      call a%Memor%dealloc(a%ncomp-1,gpene,'gpene','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,grden,'grden','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,grmom,'grmom','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,grene,'grene','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,gpomo,'gpomo','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,gpove,'gpove','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%mnode,AGradV,'AGradV','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,elvst,'elvst','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,novst,'novst','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,visten,'visten','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,cndten,'cndten','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,vgmom,'vgmom','nsc_EndElmope_ex')
      call a%Memor%dealloc(a%ndofn,gpR,'gpR','nsc_EndElmope_ex')
      call a%Memor%dealloc(a%ndofn,Subscales,'Subscales','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine GatherBase
      
      call e%gather(1_ip,elden(:,1),a%densf(:,1))
      call e%gather(1_ip,elden(:,2),a%densf(:,3))
      call e%gather(e%ndime,elmom(:,:,1),a%momen(:,:,1))
      call e%gather(e%ndime,elmom(:,:,2),a%momen(:,:,3))
      call e%gather(1_ip,elene(:,1),a%energ(:,1))
      call e%gather(1_ip,elene(:,2),a%energ(:,3))

   end subroutine

   !----------------------------------------------------------
   subroutine CalculateElmchl
      
      call nsc_ComputeElementVelocity(e,elden(:,1),elmom(:,:,1),elvel)
      call elmchl(e,1_ip,elvel,chale)
      
   end subroutine
   
   !----------------------------------------------------------
   subroutine InterpolateBase
      
      call e%interpg(1_ip,elden(:,1),gpden(1))
      call e%interpg(1_ip,elden(:,2),gpden(2))
      call e%interpg(e%ndime,elmom(:,:,1),gpmom(:,1))
      call e%interpg(e%ndime,elmom(:,:,2),gpmom(:,2))
      call e%interpg(1_ip,elene(:,1),gpene(1))
      call e%interpg(1_ip,elene(:,2),gpene(2))
      
   end subroutine

   !----------------------------------------------------------
   subroutine ComputeLinearVariables
      implicit none
      
      gpode = gpden(1)
      gpomo(:) = gpmom(:,1)
      gpoen = gpene(1)

   end subroutine

   subroutine ComputeNonlinearVariables
      implicit none
      
      gpode = gpode + a%cosgs(ielem)%a(1,e%igaus)
      gpomo(:) = gpomo(:) + a%mosgs(ielem)%a(:,1,e%igaus)
      gpoen = gpoen + a%ensgs(ielem)%a(1,e%igaus)

   end subroutine

   subroutine ComputeAuxiliarVariables
      implicit none
      
      !Advection momentum norm 
      call vecnor(gpomo,e%ndime,gpmno,2)

      invcvh = 1.0_rp / accvh 
      invgpd = 1.0_rp / gpode
      gpove = gpomo / gpode
      sqinvgpd = invgpd * invgpd
      sqgpmn = gpmno * gpmno
      sqgpvn = sqgpmn * sqinvgpd

   end subroutine

   ! Primitive auxiliary variables
   subroutine ComputeIdealAuxVariables
      implicit none
     
      acgamma = accph * invcvh
      actcn = invgpd * invcvh
      aux = acgamma - 1.0_rp
      aux_t = gpoen - (sqgpmn*invgpd/2.0_rp) 
      aux_d = aux * sqgpvn / 2.0_rp
      gptem = invgpd * invcvh * aux_t 
      gppre = (acgamma - 1.0_rp) * aux_t
      
   end subroutine

   subroutine ComputeSoundSpeed
      implicit none
     
      call nsc_ComputeSoundSpeed(accph,accvh,gptem,gpspd)
      
   end subroutine

   subroutine ComputeBaseGradients
      
      !Compute mom·grad(V)
      call ComputeAGradV(e,gpove,AGradV)

      ! Compute element variables gradient 
      call e%gradient(1_ip,elden(:,1),grden)
      call e%gradient(e%ndime,elmom(:,:,1),grmom)
      call e%gradient(1_ip,elene(:,1),grene)

      ! Compute element momentum divergence
      call e%divergence(elmom(:,:,1),divmom)

      ! Compute (mom/rho)·grad Var
      call nsc_vgvar(e,elden(:,1),AGradV,vgden)
      call nsc_vgvec(e,elmom(:,:,1),AGradV,vgmom)
      call nsc_vgvar(e,elene(:,1),AGradV,vgene)

   end subroutine

   !----------------------------------------------------------
   !Computation of Tau
   subroutine ComputeTaus
      implicit none
      
      call nsc_ComputeTau(e,gpspd,acvis*invgpd,actco*invgpd/accph,gpmno*invgpd,a%staco,chale,timom)

   end subroutine

   !----------------------------------------------------------
   subroutine ComputeResidual
      
      elext=0.0_rp
      elexh=0.0_rp
      elvst=0.0_rp
      novst=0.0_rp
      visten = 0.0_rp
      cndten = 0.0_rp

      !Compute vector of external forces
      call ProcPointer%ExternalForces   
      
      !Residual matrices to zero
      elred=0.0_rp
      elrem=0.0_rp
      elree=0.0_rp
      elcd=0.0_rp
      elcm=0.0_rp
      elce=0.0_rp

      !Convective residual
      call ProcPointer%ComputeConvectiveRes

      elred = elcd
      elrem = elcm
      elree = elce

      !Unsteady residual
      call ProcPointer%ComputeUnsteadyRes
      
      !Reaction (transformed forces) term residual
      call nsc_res_mf(e,gpode,elext,elrem)
      call nsc_res_ef(e,gpomo,elext,elree)
      call nsc_res_es(e,gpode,elexh,elree)

      !Viscous residual
      call nsc_elmvst(e,gpode,gpomo,grden,vgden,grmom,divmom,novst)
      call nsc_MolecularDiff(e,acvis,actco,visten,cndten)
      call ProcPointer%ComputeDiffusiveRes

      !External forces residual
      call ProcPointer%ComputeForcesRes

   end subroutine

   !-------------------------------------------------------------------
   !Force term   
   subroutine ComputeForces
      implicit none    
      !Compute vector of momentum external forces
      call nsc_ComputeMomentumExternalForces(e,a%grnor,a%gravi,elext)  
      !Compute energy external heat
      call nsc_ComputeEnergyExternalHeat(e,a%srce,elexh)  
   end subroutine   

   !Exact solution forces 
   subroutine AllocateExactSol
      
      call a%Memor%alloc(1,elexd,'elexd','nsc_elmope_ex')
      call a%Memor%alloc(e%ndime,elexm,'elexm','nsc_elmope_ex')
      call a%Memor%alloc(1,elexe,'elexe','nsc_elmope_ex')
      
   end subroutine

   subroutine DeallocateExactSol

      call a%Memor%dealloc(1,elexd,'elexd','nsc_elmope_ex')
      call a%Memor%dealloc(e%ndime,elexm,'elexm','nsc_elmope_ex')
      call a%Memor%dealloc(1,elexe,'elexe','nsc_elmope_ex')
      
   end subroutine

   subroutine  ExactSolutionForces    
      use Mod_Mesh    
      implicit none      
      real(rp)    :: gpcod(e%ndime)
      !Interpolate
      call e%interpg(e%ndime,e%elcod,gpcod)

      elexd(1) = 0.0_rp
      elexm = 0.0_rp
      elexe(1) = 0.0_rp

      !Compute vector of external forces  
      call exacso%nsc_ComputeSolution(e%ndime,gpcod,a)
      call exacso%nsc_GetConservativeForce(e%ndime,elexd,elexm,elexe,a)        
    
   end subroutine

   subroutine ComputeExactSolForcesRes
      implicit none    
   
      call nsc_res_dm(e,elexd(1),elred)
      call nsc_res_mf(e,1.0_rp,elexm,elrem)
      call nsc_res_dm(e,elexe(1),elree)

   end subroutine   

   !-------------------------------------------------------------------
   !Unsteady residual 
   subroutine ComputeUnsteadyRes
      implicit none    
   
      !Unsteady residual
      call nsc_res_unst(e,1_ip,gpden,a%dtinv,elred)
      call nsc_res_unst(e,e%ndime,gpmom,a%dtinv,elrem)
      call nsc_res_unst(e,1_ip,gpene,a%dtinv,elree)

   end subroutine   

   !-------------------------------------------------------------------
   !Convective residual 
   subroutine ComputeConvectiveRes
      implicit none    

      !Mass residual
      call nsc_res_dm(e,divmom,elcd)
      !Momentum residual
      call nsc_res_md(e,gpove,vgden,elcm)
      call nsc_res_mmd(e,gpove,divmom,elcm)
      call nsc_res_mmg(e,vgmom,elcm)
      !Energy residual
      call nsc_res_edvar(e,gpoen*invgpd,vgden,elce)
      call nsc_res_edvar(e,gppre*invgpd,vgden,elce)
      call nsc_res_emvar(e,gpoen*invgpd,divmom,elce)
      call nsc_res_emvar(e,gppre*invgpd,divmom,elce)

   end subroutine   

   !-------------------------------------------------------------------
   !Ideal state law residual convective terms
   subroutine ComputeConvectiveRes_i
      implicit none
      
      !Momentum residual
      call nsc_res_mdp(e,aux_d,grden,elcm)
      call nsc_res_mmp(e,aux,gpove,grmom,elcm)
      call nsc_res_mep(e,aux,grene,elcm)
      !Energy residual
      call nsc_res_edpp(e,aux_d,vgden,elce)
      call nsc_res_empp(e,aux,gpove,vgmom,elce)
      call nsc_res_eep(e,acgamma,vgene,elce)

   end subroutine   

   !----------------------------------------------------------
   !NonLinear Elements   

   subroutine AllocateHighOrderDerivatives
      
      call a%Memor%alloc(e%ndime,e%ndime,1_ip,heden,'heden','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,e%ndime,hemom,'hemom','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,1_ip,heene,'heene','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,dvist,'dvist','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,hessV,'hessV','nsc_EndElmope_ex')
      call a%Memor%alloc(e%mnode,laplV,'laplV','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine DeallocateHighOrderDerivatives
      
      call a%Memor%dealloc(e%ndime,e%ndime,1_ip,heden,'heden','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,e%ndime,hemom,'hemom','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,1_ip,heene,'heene','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,dvist,'dvist','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,hessV,'hessV','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%mnode,laplV,'laplV','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine InGaussNonLinear
      implicit none
      
      call e%elmder
      call e%elmhes

   end subroutine

   !-------------------------------------------------------------------
   subroutine ComputeDiffusiveRes
      implicit none    

      !Hessian
      call nsc_hessian(e,1_ip,elden(:,1),heden)
      call nsc_hessian(e,e%ndime,elmom(:,:,1),hemom)
      call nsc_hessian(e,1_ip,elene(:,1),heene)

      !Momentum residual
      dvist= 0.0_rp
      call nsc_res_mdp(e,(2/3)*vgden*sqinvgpd,grden,dvist)
      call nsc_dvist_dh(e,gpove*invgpd/3,heden(:,:,1),dvist)
      call nsc_dvist_dd(e,gpove*sqinvgpd,grden,dvist)
      call nsc_dvist_dl(e,gpove*invgpd,heden(:,:,1),dvist)
      call nsc_dvist_dm1(e,sqinvgpd/3,grden,grmom,dvist)
      call nsc_dvist_dm2(e,2*sqinvgpd,grden,grmom,dvist)
      call nsc_res_mdp(e,-sqinvgpd*divmom/3,grden,dvist)
      call nsc_dvist_mh(e,hemom*invgpd/3,dvist)
      call nsc_dvist_ml(e,hemom*invgpd,dvist)
      elrem = elrem + matmul(visten,dvist)

      !Energy residual
      elvst = matmul(novst,transpose(visten))
      call nsc_res_diff_ed(e,gpove*sqinvgpd,elvst,grden,elree)
      call nsc_res_diff_em(e,elvst*invgpd,grmom,elree)
      call nsc_res_diff_ev(e,gpove,visten,dvist,elree)

   end subroutine   

   !-------------------------------------------------------------------
   !Ideal state law residual diffusive terms
   subroutine ComputeDiffusiveRes_i
      implicit none
      
      !Terms of the energy equation
      dhflx = 0.0_rp
      call nsc_dhflx_edd(e,actcn*(3*sqgpvn*invgpd-2*gpoen*sqinvgpd),grden,dhflx)
      call nsc_dhflx_edl(e,actcn*(gpoen*invgpd-sqgpvn),heden(:,:,1),dhflx)
      call nsc_dhflx_ede(e,2*invgpd*actcn,grden,grene,dhflx)
      call nsc_dhflx_edm(e,4*invgpd*actcn,gpove,grden,grmom,dhflx)
      call nsc_dhflx_emm(e,invgpd*actcn,grmom,dhflx)
      call nsc_dhflx_emh(e,actcn,gpove,hemom,dhflx)
      call nsc_dhflx_eeh(e,actcn,heene(:,:,1),dhflx)
      elree(1) = elree(1) + actco*dhflx

   end subroutine   

   !-------------------------------------------------------------------

   subroutine TotalResidual
      implicit none
     
      gpR(1) = elred(1) 
      gpR(2:e%ndime+1) = elrem
      gpR(a%ndofn) = elree(1)

   end subroutine   

   !Orthogonal Projection Subscales
   
   subroutine AllocateOSS

      call a%Memor%alloc(a%ndofn,e%mnode,elrep,'elrep','nsc_EndElmope_ex')
      call a%Memor%alloc(a%ndofn,reprj,'reprj','nsc_EndElmope_ex')

   end subroutine
   
   subroutine DeallocateOSS

      call a%Memor%dealloc(a%ndofn,e%mnode,elrep,'elrep','nsc_EndElmope_ex')
      call a%Memor%dealloc(a%ndofn,reprj,'reprj','nsc_EndElmope_ex')

   end subroutine

   subroutine GatherOSS
      
      call e%gather(a%ndofn,elrep,a%repro)

   end subroutine

   subroutine InterpolateOSS

      call e%interpg(a%ndofn,elrep,reprj)
      
   end subroutine


   subroutine ResidualOSS
      implicit none
     
      gpR(1) = elred(1) - reprj(1) 
      gpR(2:e%ndime+1) = elrem(:) - reprj(2:e%ndime+1)
      gpR(a%ndofn) = elree(1) - reprj(a%ndofn) 

   end subroutine   

   !-------------------------------------------------------------------

   subroutine ComputeStaticSubscales
      implicit none

      Subscales(1)            = timom(1)*gpR(1) 
      Subscales(2:e%ndime+1)  = timom(2)*gpR(2:e%ndime+1) 
      Subscales(a%ndofn)      = timom(3)*gpR(a%ndofn) 

      !The calculated residual is positive, must change sign
      a%cosgs(ielem)%a(1,e%igaus)   = -Subscales(1)
      a%mosgs(ielem)%a(:,1,e%igaus) = -Subscales(2:e%ndime+1)
      a%ensgs(ielem)%a(1,e%igaus)   = -Subscales(a%ndofn)

   end subroutine   

   !---------------------------------------------------------------
   !Loop by points
   subroutine PointsLoop
      
     do ipoin = 1,npoin

         call ProcPointer%ComputeInPoint

     enddo
      
   end subroutine

   !---------------------------------------------------------------
   !Residual Projection 
   
   subroutine AllocateRep
      implicit none
      call a%Memor%alloc(a%ndofn,e%mnode,elres,'elres','nsc_EndElmope_ex')
      call a%Memor%alloc(a%ndofn,npoin,wrepro,'wrepro','nsc_EndElmope_ex')
      
     wrepro(:,:) = 0.0_rp

   end subroutine
   
   subroutine DeallocateRep
      implicit none
      call a%Memor%dealloc(a%ndofn,e%mnode,elres,'elres','nsc_EndElmope_ex')
      call a%Memor%dealloc(a%ndofn,npoin,wrepro,'wrepro','nsc_EndElmope_ex')

   end subroutine
   
   subroutine ElmatsToZeroRep
      
      elres(:,:) = 0.0_rp

   end subroutine
   
   subroutine GpresToElres
      
      elres(1,1:e%pnode) = elres(1,1:e%pnode) + e%shape(1:e%pnode,e%igaus)*elred(1)*dvol
      do idime = 1,e%ndime
         elres(idime+1,1:e%pnode) = elres(idime+1,1:e%pnode) + e%shape(1:e%pnode,e%igaus)*elrem(idime)*dvol
      enddo
      elres(e%ndime+2,1:e%pnode) = elres(e%ndime+2,1:e%pnode) + e%shape(1:e%pnode,e%igaus)*elree(1)*dvol

   end subroutine
   
   subroutine AssemblyResidual
      
      wrepro(:,e%lnods(1:e%pnode)) = wrepro(:,e%lnods(1:e%pnode)) + elres(:,1:e%pnode)

   end subroutine

   subroutine ResidualProjSol
      
       a%repro(:,ipoin) = vmass(ipoin)*wrepro(:,ipoin)
      
   end subroutine
   
   !----------------------------------------------------------
   !Gradient Orthogonal Projection calculation 
   subroutine AllocateGradProj
      
      call a%Memor%alloc(e%ndime+1,e%ndime,npoin,wrhg,'wrhg','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime+1,e%ndime,e%mnode,elrhg,'elrhg','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elrgm,'elrgm','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%mnode,elrge,'elrge','nsc_EndElmope_ex')

      wrhg = 0.0_rp

   end subroutine

   subroutine DeallocateGradProj
      
      call a%Memor%dealloc(e%ndime+1,e%ndime,npoin,wrhg,'wrhg','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime+1,e%ndime,e%mnode,elrhg,'elrhg','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elrgm,'elrgm','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%mnode,elrge,'elrge','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine ZeroGradProj
      
      !ElmatsToZero
      elrhg=0.0_rp
      elrgm=0.0_rp
      elrge=0.0_rp

   end subroutine

   !----------------------------------------------------------
   subroutine CalculateGradProj
      
       call grad_proj(e,e%ndime,dvol,grmom,elrgm)
       call grad_proj(e,1_ip,dvol,grene,elrge)
      
   end subroutine

   !----------------------------------------------------------
   subroutine AssemblyGradProj
      
      elrhg(1:e%ndime,:,1:e%pnode) = elrgm(1:e%ndime,:,1:e%pnode) &
            + elrhg(1:e%ndime,:,1:e%pnode)
   
      elrhg(e%ndime+1,:,1:e%pnode) = elrge(:,1:e%pnode) &
            + elrhg(e%ndime+1,:,1:e%pnode)

      !GlobalAssembly

     wrhg(:,:,e%lnods(1:e%pnode)) = wrhg(:,:,e%lnods(1:e%pnode)) + elrhg(:,:,1:e%pnode)
      
   end subroutine

   !----------------------------------------------------------
   subroutine GradientProjSol
      
       a%grprj(1:e%ndime+1,:,ipoin) = vmass(ipoin)*wrhg(1:e%ndime+1,:,ipoin)
      
   end subroutine

   !----------------------------------------------------------
   !Shock capturing
   subroutine ComputeDCResidual

      dcrem = elrem
      dcree(1) = elree(1)

      call ProcPointer%ComputeDCUnsteadyRes

   end subroutine

   subroutine ComputeDCUnsteadyRes
      implicit none    

      call nsc_res_unst(e,e%ndime,gpmom,a%dtinv,dcrem)
      call nsc_res_unst(e,1_ip,gpene,a%dtinv,dcree(1))

   end subroutine

   subroutine ComputeResDiff

      call nsc_ArtificialVis(e,a%shock,chale,dcrem,grmom,arvis)
      call nsc_ArtificialCnd(e,a%shock,chale,dcree(1),grene,artco)

   end subroutine

   subroutine ComputeGradDiff

      orthm = grmom - gmprj
      orthe = grene - geprj
      
      call nsc_GradOrthVis(e,a%shock,chale,gpmno*invgpd,orthm,grmom,arvis)
      call nsc_GradOrthCnd(e,a%shock,chale,gpmno*invgpd,orthe,grene,artco)

   end subroutine

   !----------------------------------------------------------
   subroutine OutputSCDiffusion

      a%scdiffarray(1,ielem)%a(e%igaus) = arvis
      a%scdiffarray(2,ielem)%a(e%igaus) = artco

   end subroutine

   subroutine OutputSGDiffusion
      implicit none
      
      a%sgdiffarray(1,ielem)%a(e%igaus) = sqgpvn*timom(2)
      a%sgdiffarray(2,ielem)%a(e%igaus) = sqgpvn*timom(3)

   end subroutine

   subroutine OutputResidual
      implicit none
      
      a%coResGP(ielem)%a(e%igaus) = -elred(1)
      a%moResGP(ielem)%a(:,e%igaus) = -elrem
      a%enResGP(ielem)%a(e%igaus) = -elree(1)

   end subroutine

   !----------------------------------------------------------
   !Residual based Shock Capturing
   subroutine AllocateResDC
      
      call a%Memor%alloc(e%ndime,dcrem,'dcrem','nsc_EndElmope_ex')
      call a%Memor%alloc(1,dcree,'dcree','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine DeallocateResDC
      
      call a%Memor%dealloc(e%ndime,dcrem,'dcrem','nsc_EndElmope_ex')
      call a%Memor%dealloc(1,dcree,'dcree','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   !Gradient Orthogonal Projection Shock Capturing
   subroutine AllocateGradProjDC
      
      call a%Memor%alloc(e%ndime,e%ndime,orthm,'orthm','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,gmprj,'gmprj','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgmp,'elgmp','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,orthe,'orthe','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,geprj,'geprj','nsc_EndElmope_ex')
      call a%Memor%alloc(e%ndime,e%mnode,elgep,'elgep','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine DeallocateGradProjDC
      
      call a%Memor%dealloc(e%ndime,e%ndime,orthm,'orthm','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,gmprj,'gmprj','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgmp,'elgmp','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,orthe,'orthe','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,geprj,'geprj','nsc_EndElmope_ex')
      call a%Memor%dealloc(e%ndime,e%mnode,elgep,'elgep','nsc_EndElmope_ex')

   end subroutine

   !----------------------------------------------------------
   subroutine GatherGradProjDC
      
      do idime=1,e%ndime
         elgmp(idime,:,1:e%pnode) = a%grprj(idime,:,e%lnods(1:e%pnode))
      enddo
      elgep(:,1:e%pnode) = a%grprj(e%ndime+1,:,e%lnods(1:e%pnode))

   end subroutine

   !----------------------------------------------------------
   subroutine InterpolateGradProjDC
      
      do idime=1,e%ndime
         gmprj(idime,:) = matmul(elgmp(idime,:,:),e%shape(1:e%pnode,e%igaus))
      enddo
      geprj = matmul(elgep,e%shape(1:e%pnode,e%igaus))
      
   end subroutine

   !--------------------------------------------------------------------
   !Multiple type of elements
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
   !---------------------------------------------------------------
   !ALE
   subroutine AllocMeshVeloc
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elmve,'elmve','nsc_elmope_ex')
      call a%Memor%alloc(e%ndime,gpmve,'gpmve','nsc_elmope_ex')
      call a%Memor%alloc(e%mnode,ALEGradV,'ALEGradV','nsc_elmope_ex')
      call a%Memor%alloc(e%ndime,mvgmom,'mvgmom','nsc_elmope_ex')
   end subroutine
   
   subroutine DeallocMeshVeloc
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,elmve,'elmve','nsc_elmope_ex')
      call a%Memor%dealloc(e%ndime,gpmve,'gpmve','nsc_elmope_ex')
      call a%Memor%dealloc(e%mnode,ALEGradV,'ALEGradV','nsc_elmope_ex')
      call a%Memor%dealloc(e%ndime,mvgmom,'mvgmom','nsc_elmope_ex')
   end subroutine

   subroutine GatherMeshVeloc
      implicit none
      
      call e%gather(e%ndime,elmve,meshve(:,:,1))
   end subroutine
   
   subroutine InterpolateMeshVeloc
      implicit none
      
      call e%interpg(e%ndime,elmve,gpmve)
   end subroutine
   
   subroutine ComputeALEGradients
      implicit none
      
      !Compute mom·grad(V)
      call ComputeAGradV(e,gpmve,ALEGradV)
      ! Compute (mom/rho)·grad Var
      call nsc_vgvar(e,elden(:,1),ALEGradV,mvgden)
      call nsc_vgvec(e,elmom(:,:,1),ALEGradV,mvgmom)
      call nsc_vgvar(e,elene(:,1),ALEGradV,mvgene)
   end subroutine

   subroutine ALEConvectiveContribution
      implicit none
      
      !Mass residual
      call nsc_res_dm(e,-mvgden,elcd)
      !Momentum residual
      call nsc_res_mmg(e,-mvgmom,elcm)
      !Energy residual
      call nsc_res_eep(e,1.0_rp,-mvgene,elce)
   end subroutine

end module


subroutine nsc_EndElmope_ex(NSCompExplicitProblem,task)
   use Mod_nsc_EndElmope_ex
   use Mod_NSCompressibleExplicit
   implicit none
   class(NSCompressibleExplicitProblem), target :: NSCompExplicitProblem
   character(6) :: task
   
   a=>NSCompExplicitProblem
   itask = task

   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !We force closed rule for lumping
   call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'nsc_EndElmope_ex')

   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)

   !Hook
   call ProcHook%Initializations
   
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetVmass(vmass)

   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)    
      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Elmats to Zero
      call ProcHook%ElmatsToZero
      
      !Hook
      call ProcHook%Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      !Hook
      call ProcHook%PreGauss

      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Hook
         call ProcHook%Interpolates
         
         !Hook
         call ProcHook%PhysicalProp          

         !Operating variables      
         call ProcHook%ComputeVariables
         
         !Compute the stability parameters
         call ProcHook%ComputeTaus
            
         !InGaussElmats
         call ProcHook%InGaussElmats
         
         !InGaussElmats Assembly
         call ProcHook%InGaussElmatsAssembly
         
      enddo gauss_points
      
      !Assembly 
      call ProcHook%Assembly
      
   enddo elements
   
   !Hook
   call ProcHook%Finalizations
   
  
   !Operations to be done after the Elemental Loop
   call ProcHook%PostLoop
   
    !Element
   call a%Mesh%ElementDeAlloc(e,a%Memor,a%EndLoopQuadrature,'nsc_EndElmope_ex')

   
end subroutine
