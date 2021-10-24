module Mod_UPSolids
   use typre
   use Mod_Solids
   implicit none
   private
   public UPSolidsProblem,zesup

   real(rp), parameter   :: &
     zesup = epsilon(0.0_rp)        ! zero

   type,abstract, extends(SolidsProblem) ::  UPSolidsProblem

      ! Physical problem
      real(rp), allocatable :: &
              devstrain(:,:)

      real(rp), allocatable :: &
              j2(:)                 !J2 stresses

      type(r2p), allocatable :: &
              j2_g(:)                    !J2 stresses

      ! Physical problem
      real(rp), allocatable :: &
              press_cp(:)

      ! Numerical treatment
      integer(ip) ::&
        mtrit      = 99,&                    ! Maximum Number of iterations of tracking
        kfl_repro  = 99,&                    ! Stabilization based on residual projection
        kfl_trasg  = 99,&                    ! Tracking of subgrid scale
        kfl_nolsg  = 99,&                    ! Non-linearity of the subscales
        kfl_tacsg  = 99                      ! Time accuracy of subscales

      !Stabilization parameters for SUP Solids
      real(rp) :: tau_u = 1.0_rp             !Stab Coef for displacements
      real(rp) :: tau_p = 1.0_rp             !Stab Coef for pressure

      real(rp) :: kfl_PrbCharL = 1.0_rp      !Characteristic length for tau_u

      type(r3p), allocatable :: &
         u_sgs(:)                            ! Displacement subgrid scales
      type(r1p), allocatable :: &
         p_sgs(:)                            ! Pressure subgrid scales

      type(r2p), allocatable :: &
          residualU(:)                       ! Velocity residual at gp for postprocess
      type(r1p), allocatable :: &
          residualP(:)                       ! Pressure residual at gp for postprocess

      !Auxiliar working arrays
      real(rp), allocatable :: &
         repro(:,:)                          ! Residual Projection

      logical :: kfl_printJ2Stresses  = .false.   !Print J2 stresses duh
      logical :: kfl_printResiduals   = .false.   !Print residual duh

   contains

      !Basic procedures
      procedure :: SetExmod                 => sldup_SetExmod
      procedure :: SetNdofn                 => sldup_SetNdofn
      procedure :: SetNdofbc                => sldup_SetNdofbc

      procedure :: SolidSpecificTurnon      => sldup_turnon
      procedure :: SolidSpecificTurnof      => sldup_turnof 
      procedure :: UPSolidSpecificTurnof    => sldup_Specificturnof
      procedure :: TurnofExtend             => sldup_NULLSUB
      procedure :: UPModelTurnon            => sldup_turnonModel
      procedure :: SpecificIniunk           => sldup_iniunk
      procedure :: SolidSpecificReaous      => sldup_reaous
      procedure :: ReaousExtend             => sldup_NULLSUBitask

      procedure :: SpecificEndite           => sldup_endite
      procedure :: UpdateDynamicComponents  => sldup_updateDynamicComp

      procedure :: SolidSpecificMemall      => sldup_memall
      procedure :: MemallExtend             => sldup_memallSpecific
      procedure :: InnerResiduals           => sldup_InnerResiduals
      procedure :: Cvgunk                   => sldup_cvgunk
      procedure :: SolidSpecificReampi      => sldup_reampi
      procedure :: SolidSpecificReanut      => sldup_reanut
      procedure :: SolidSpecificOutput      => sldup_output
      procedure :: OutputExtend             => sldup_NULLSUBitask
      procedure :: SolidSpecificRestart     => sldup_restar
      procedure :: PointTracking            => sldup_outtpo
      procedure :: SolidSpecificRefine      => sldup_refine
      procedure :: SolidSpecificRefineGauss => sldup_refineGauss

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
      procedure :: EndBouope                => sldup_EndBouope
      procedure :: Bouope                   => sld_bouop0
      procedure :: ModifyBouopeRHS          => sldup_modifyBouopeRHS
      procedure :: ModifyPointForces        => sldup_modifyPointForces

      procedure :: GetMatrixOrganization => sldup_GetMatrixOrganization
      
   end type

   interface
     
      subroutine sldup_SetExmod(a)
         import 
         implicit none         
         class(UPSolidsProblem) :: a      
      end subroutine
      
      subroutine sldup_SetNdofn(a)
         import 
         implicit none         
         class(UPSolidsProblem) :: a      
      end subroutine
      
      subroutine sldup_SetNdofbc(a)
         import 
         implicit none         
         class(UPSolidsProblem) :: a      
      end subroutine
     
      subroutine sld_reaphy(a,itask)
         import 
         implicit none         
         integer(ip) :: itask
         class(UPSolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reanut(a,itask)
         import 
         implicit none         
         integer(ip) :: itask
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine sldup_reaous(a,itask)
         import 
         implicit none         
         integer(ip) :: itask
         class(UPSolidsProblem) :: a      
      end subroutine
   
      subroutine sld_reaMPI(a)
         import 
         implicit none         
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine sldup_memall(a)
         import
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine sldup_memallSpecific(a)
         import
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine
      
      subroutine sld_begste(a)
         import
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine
      
      subroutine sldup_endite(a,itask)
         use typre
         import
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip) :: itask     
      end subroutine
      
      subroutine sldup_cvgunk(a,itask)
         use typre
         import
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip), intent(in) :: itask   
      end subroutine
      
      subroutine sldup_InnerResiduals(a,resld,srsld,prsld)
          use typre
          import
          implicit none
          class(UPSolidsProblem) :: a
          real(rp)    :: resld
          real(rp), optional :: srsld,prsld
      end subroutine

      subroutine sldup_reampi(a)
         use typre
         import 
         implicit none         
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine sldup_reanut(a,itask)
         use typre
         import 
         implicit none         
         integer(ip) :: itask
         class(UPSolidsProblem) :: a      
      end subroutine
   
      subroutine sldup_output(a,itask)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip) :: itask    
      end subroutine

      subroutine sldup_outtpo(a)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
      end subroutine
      
      subroutine sld_endste(a,itask)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip) :: itask    
      end subroutine
            
      subroutine sldup_iniunk(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a
      end subroutine 

      subroutine Elmope(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine
      
      subroutine sld_bouop0(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a     
      end subroutine
      
      subroutine sldup_turnon(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine sldup_turnof(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine sldup_Specificturnof(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a      
      end subroutine

      subroutine EndElmope(a,task)
         import 
         implicit none
         class(UPSolidsProblem) :: a
         character(6) :: task
      end subroutine
      
      subroutine sldup_EndBouope(a)
         import 
         implicit none
         class(UPSolidsProblem) :: a
      end subroutine

      subroutine sld_restar(a,itask)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip), intent(in) :: itask       
         !Restart calculations
      end subroutine

      subroutine sldup_restar(a,itask)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip), intent(in) :: itask
         !Restart calculations
      end subroutine

      subroutine sld_refine(a,itask)
         import 
         implicit none
         class(UPSolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sld_GetRefinementCriteria(a,markel)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip) :: markel(*)
      end subroutine

      subroutine sldup_refine(a,itask)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sldup_refineGauss(a,itask)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         character(6) :: itask
      end subroutine

      subroutine sld_Reabcs(a,itask,kflag)
         use typre
         import 
         implicit none
         class(UPSolidsProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      end subroutine

      subroutine sldup_updateDynamicComp(a)
          import 
          implicit none
          class(UPSolidsProblem) :: a
      end subroutine

      subroutine sldsup_GetStress(sld,ndime,sz,stress_t)
         use typre 
         import 
         implicit none
         class(UPSolidsProblem) :: sld
         integer(ip)          :: ndime,sz
         real(rp)             :: stress_t(ndime,ndime)

      end subroutine

  end interface

  contains

   subroutine sldup_NULLSUB(a)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
   
   end subroutine

   subroutine sldup_NULLSUBrestart(a,itask)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine

   subroutine sldup_NULLSUBitaskc(a,itask)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      character(6) :: itask
      
   end subroutine

   subroutine sldup_NULLSUBitaskintentin(a,itask)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      integer(ip), intent(in) :: itask
      
   end subroutine
   
    subroutine sldup_NULLSUBitaskkflag(a,itask,kflag)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
      
   end subroutine

    subroutine sldup_NULLSUBdtinv(a,dtinv)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      real(rp) :: dtinv
      
   end subroutine

   subroutine sldup_modifyPointForces(a,ndof,pfsize,pf)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      integer(ip), intent(in) :: ndof,pfsize
      real(rp), intent(inout) :: pf(ndof,pfsize)

      pf = pf/(2.0_rp*a%mu)
   
   end subroutine

   subroutine sldup_modifyBouopeRHS(a,ndof,mnode,rhs)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      integer(ip), intent(in) :: ndof,mnode
      real(rp), intent(inout) :: rhs(ndof,mnode)

      rhs = rhs/(2.0_rp*a%mu)
   
   end subroutine

   subroutine sldup_GetMatrixOrganization(a,u1,uf,s1,sf,p1,bcstart)
      use typre
      implicit none
      class(UPSolidsProblem)    :: a
      integer(ip), intent(inout) :: u1,uf,s1,sf,p1,bcstart
      integer(ip)                :: nd

      !Update unknowns
      call a%Mesh%GetNdime(nd)

      s1 = 0
      sf = 0
      u1 = 1
      uf = nd
      p1 = uf+1
      bcstart = 0

   end subroutine

   subroutine sldup_SetInitialConditions(a,unkno)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
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
   
   subroutine sldup_GetInitialConditions(a,unkno)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
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

   subroutine sldup_GetUnkno(a,unkno)
      use typre
      implicit none
      class(UPSolidsProblem) :: a
      real(rp) :: unkno(:,:)
            
   end subroutine

  subroutine sldup_turnonModel(a)
      use typre
      implicit none
      class(UPSolidsProblem) :: a      

  end subroutine

  subroutine sldup_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(UPSolidsProblem) :: a      
      integer(ip) :: itask

  end subroutine
   

end module Mod_UPSolids
