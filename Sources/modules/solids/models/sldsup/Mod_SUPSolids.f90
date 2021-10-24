module Mod_SUPSolids
   use typre
   use Mod_UPSolids
   implicit none
   private
   public SUPSolidsProblem,zesup

   type,abstract, extends(UPSolidsProblem) ::  SUPSolidsProblem

      !norms of sigma and devstrain
      real(rp) :: &
              signorm = 0.0_rp ,&
              devnorm = 0.0_rp

      ! Physical problem
      real(rp), allocatable :: &
              sigma_cp(:,:)

      !Stabilization parameters for SUP Solids
      real(rp) :: tau_s = 0.2_rp             !Stab Coef for stresses

      logical :: kfl_printDevStrain   = .false.   !Used to know how to print stress and strain
      logical :: kfl_useSecantModulus = .false.   !Calculate tau_u with ||s||/||e_dev||

      type(r2p), allocatable :: &
         s_sgs(:)                            ! Sigma subgrid scales

      type(r2p), allocatable :: &
          residualS(:)               ! Stress residual at gp for postprocess

   contains

      !Basic procedures
      procedure :: SetExmod                => sldsup_SetExmod
      procedure :: SetNdofn                => sldsup_SetNdofn
      procedure :: SetNdofbc               => sldsup_SetNdofbc

      procedure :: SolidSpecificTurnon     => sldsup_turnon
      procedure :: TurnofExtend            => sldsup_turnof 
      procedure :: UPSolidSpecificTurnof   => sldsup_NULLSUB
      procedure :: SUPSolidSpecificTurnof  => sldsup_Specificturnof
      procedure :: SUPModelTurnon          => sldsup_turnonModel
      procedure :: SpecificIniunk          => sldsup_iniunk
      procedure :: ReaousExtend            => sldsup_reaous

      procedure :: SpecificEndite          => sldsup_endite
      procedure :: UpdateDynamicComponents => sldsup_updateDynamicComp

      !NULLSUBS
      procedure :: SpecificExaerr    => sldsup_NULLSUB
      procedure :: SpecificUpdbcs    => sldsup_NULLSUB
      procedure :: SpecificSolite    => sldsup_NULLSUBitask
      procedure :: SpecificCrankNicolsonEndste => sldsup_NULLSUB
      procedure :: Getste            => sldsup_NULLSUBdtinv
       
      procedure :: SpecificOuterr           => sldsup_outerr
      procedure :: MemallExtend             => sldsup_memall
      procedure :: MemallExtend2            => sldsup_memallSpecific
      procedure :: InnerResiduals           => sldsup_InnerResiduals
      procedure :: Cvgunk                   => sldsup_cvgunk
      procedure :: SolidSpecificReampi      => sldsup_reampi
      procedure :: SolidSpecificReanut      => sldsup_reanut
      procedure :: OutputExtend             => sldsup_output
      procedure :: SolidSpecificRestart     => sldsup_restar
      procedure :: PointTracking            => sldsup_outtpo
      procedure :: SolidSpecificRefine      => sldsup_refine
      procedure :: SolidSpecificRefineGauss => sldsup_refineGauss

      !Parent procedures
      procedure :: SpecificReaphy  => sld_reaphy
      procedure :: SpecificReanut  => sld_reanut
      procedure :: SpecificReaMPI  => sld_reaMPI
      procedure :: SpecificReabcs  => sld_Reabcs
      procedure :: SpecificBegste  => sld_begste
      procedure :: SpecificEndste  => sld_endste
      procedure :: SpecificRestart => sld_restar
      procedure :: SpecificRefine  => sld_refine
      procedure :: GetRefinementCriteria => sld_GetRefinementCriteria

      !Boundary operations 
      procedure :: EndBouope                => sldsup_EndBouope
      procedure :: Bouope                   => sld_bouop0
      procedure :: ModifyBouopeRHS          => sldsup_NullsubRhsIn
      procedure :: ModifyPointForces        => sldsup_NullsubPointForcesIn

      ! ROM
      procedure :: SetInitialConditions => sldsup_SetInitialConditions
      procedure :: GetInitialConditions => sldsup_GetInitialConditions
      procedure :: SpecificGetUnkno     => sldsup_GetUnkno

      !ALE
      procedure :: GetFoldedElements     => sldsup_GetFoldedElements
      procedure :: GetMatrixOrganization => sldsup_GetMatrixOrganization
      
   end type

   interface
     
      subroutine sldsup_SetExmod(a)
         import 
         implicit none         
         class(SUPSolidsProblem) :: a      
      end subroutine
      
      subroutine sldsup_SetNdofn(a)
         import 
         implicit none         
         class(SUPSolidsProblem) :: a      
      end subroutine
      
      subroutine sldsup_SetNdofbc(a)
         import 
         implicit none         
         class(SUPSolidsProblem) :: a      
      end subroutine
     
      subroutine sld_reaphy(a,itask)
         import 
         implicit none         
         integer(ip) :: itask
         class(SUPSolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reanut(a,itask)
         import 
         implicit none         
         integer(ip) :: itask
         class(SUPSolidsProblem) :: a      
      end subroutine
   
      subroutine sldsup_reaous(a,itask)
         import 
         implicit none         
         integer(ip) :: itask
         class(SUPSolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reaMPI(a)
         import 
         implicit none         
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine sldsup_outerr(a)
         import
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine sldsup_memall(a)
         import
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine sldsup_memallSpecific(a)
         import
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine
      
      subroutine sld_begste(a)
         import
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine
      
      subroutine sldsup_endite(a,itask)
         use typre
         import
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip) :: itask     
      end subroutine
      
      subroutine sldsup_cvgunk(a,itask)
         use typre
         import
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip), intent(in) :: itask   
      end subroutine
      
      subroutine sldsup_InnerResiduals(a,resld,srsld,prsld)
          use typre
          import
          implicit none
          class(SUPSolidsProblem) :: a
          real(rp)    :: resld
          real(rp), optional :: srsld,prsld
      end subroutine

      subroutine sldsup_reampi(a)
         use typre
         import 
         implicit none         
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine sldsup_reanut(a,itask)
         use typre
         import 
         implicit none         
         integer(ip) :: itask
         class(SUPSolidsProblem) :: a      
      end subroutine
   
      subroutine sldsup_output(a,itask)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip) :: itask    
      end subroutine

      subroutine sldsup_outtpo(a)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
      end subroutine
      
      subroutine sld_endste(a,itask)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip) :: itask    
      end subroutine
            
      subroutine sldsup_iniunk(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a
      end subroutine 

      subroutine Elmope(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine
      
      subroutine sld_bouop0(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a     
      end subroutine
      
      subroutine sldsup_turnon(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine sldsup_turnof(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine sldsup_Specificturnof(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a      
      end subroutine

      subroutine EndElmope(a,task)
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         character(6) :: task
      end subroutine
      
      subroutine sldsup_EndBouope(a)
         import 
         implicit none
         class(SUPSolidsProblem) :: a
      end subroutine

      subroutine sld_restar(a,itask)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip), intent(in) :: itask       
         !Restart calculations
      end subroutine

      subroutine sldsup_restar(a,itask)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip), intent(in) :: itask
         !Restart calculations
      end subroutine

      subroutine sld_refine(a,itask)
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sld_GetRefinementCriteria(a,markel)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

      subroutine sldsup_refine(a,itask)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sldsup_refineGauss(a,itask)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sld_Reabcs(a,itask,kflag)
         use typre
         import 
         implicit none
         class(SUPSolidsProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      end subroutine

      subroutine sldsup_updateDynamicComp(a)
          import 
          implicit none
          class(SUPSolidsProblem) :: a
      end subroutine

      subroutine sldsup_GetStress(sld,ndime,sz,stress_t)
         use typre 
         import 
         implicit none
         class(SUPSolidsProblem) :: sld
         integer(ip)          :: ndime,sz
         real(rp)             :: stress_t(ndime,ndime)

      end subroutine

  end interface

  contains

   subroutine sldsup_NULLSUB(a)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
   
   end subroutine

   subroutine sldsup_NULLSUBrestart(a,itask)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine

   subroutine sldsup_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      integer(ip) :: itask
      
   end subroutine
   
   subroutine sldsup_NULLSUBitaskc(a,itask)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      character(6) :: itask
      
   end subroutine

   subroutine sldsup_NULLSUBitaskintentin(a,itask)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine
   
    subroutine sldsup_NULLSUBitaskkflag(a,itask,kflag)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
      
   end subroutine

    subroutine sldsup_NULLSUBdtinv(a,dtinv)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      real(rp) :: dtinv
      
   end subroutine

   subroutine sldsup_NullsubPointForcesIn(a,ndof,pfsize,pf)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      integer(ip), intent(in) :: ndof,pfsize
      real(rp), intent(inout) :: pf(ndof,pfsize)
   
   end subroutine

   subroutine sldsup_NullsubRhsIn(a,ndof,mnode,rhs)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      integer(ip), intent(in) :: ndof,mnode
      real(rp), intent(inout) :: rhs(ndof,mnode)
   
   end subroutine

   subroutine sldsup_GetMatrixOrganization(a,u1,uf,s1,sf,p1,bcstart)
      use typre
      implicit none
      class(SUPSolidsProblem)    :: a
      integer(ip), intent(inout) :: u1,uf,s1,sf,p1,bcstart
      integer(ip)                :: nd,tn

      !Update unknowns
      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2

      !u1 = 1
      !uf = nd
      !s1 = uf+1
      !sf = uf+tn
      !p1 = sf+1
      !bcstart = 0

      s1 = 1
      sf = tn
      u1 = sf+1
      uf = sf+nd
      p1 = uf+1
      bcstart = tn

   end subroutine

   subroutine sldsup_SetInitialConditions(a,unkno)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      real(rp) :: unkno(:,:,:)
      
      integer(ip) :: npoin,ndime,ispos,ipoin,itwost,ndofn
      integer(ip) :: u1,uf,s1,sf,p1,bc
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%GetNdofn(ndofn)

      call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)
      
      do ipoin = 1,npoin
         a%disp(1:ndime,ipoin,3) = unkno(u1:uf,ipoin,1)
      enddo
      
      a%disp(:,:,1) = a%disp(:,:,3)
      
      a%disp(:,:,2) = a%disp(:,:,3)
      
      !Higher order integration schemes
      if (a%ncomp>3) then
         do itwost = 1,a%ncomp-3
            do ipoin = 1,npoin
               a%disp(1:ndime,ipoin,3+itwost) = unkno(u1:uf,ipoin,itwost+1)
            enddo
         enddo
      endif
   end subroutine
   
   subroutine sldsup_GetInitialConditions(a,unkno)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      real(rp) :: unkno(:,:,:)
      integer(ip) :: u1,uf,s1,sf,p1,bc
      
      integer(ip) :: npoin,ndime,ispos,ipoin,itwost,ndofn
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%GetNdofn(ndofn)

      call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)
      
      do ipoin = 1,npoin
         unkno(u1:uf,ipoin,1) = a%disp(1:ndime,ipoin,3) 
      enddo
      
      !Higher order integration schemes
      if (a%ncomp>3) then
         do itwost = 1,a%ncomp-3
            do ipoin = 1,npoin
               unkno(u1:uf,ipoin,itwost+1) = a%disp(1:ndime,ipoin,3+itwost)
            enddo
         enddo
      endif
   end subroutine

   subroutine sldsup_GetUnkno(a,unkno)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      real(rp) :: unkno(:,:)
            
   end subroutine

   subroutine sldsup_GetFoldedElements(a,kfl_folded)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a
      logical::kfl_folded

      kfl_folded = a%kfl_foldedElements

  end subroutine

  subroutine sldsup_turnonModel(a)
      use typre
      implicit none
      class(SUPSolidsProblem) :: a      

  end subroutine

end module Mod_SUPSolids
