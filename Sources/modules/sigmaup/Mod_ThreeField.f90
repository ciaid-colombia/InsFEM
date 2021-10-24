module Mod_ThreeField
   use typre
   use Mod_Mesh
   use Mod_Timer
   use Mod_NavierStokes, only : NavierStokesProblem
   use Mod_ParallelSystemInterface
   private
   public ThreeFieldNSProblem,ThreeFieldNSProblem_Const
   
   type, extends(NavierStokesProblem) :: ThreeFieldNSProblem
   
      ! Physical problem
      real(rp), allocatable :: &
            sigma(:,:,:), &                  !Stress field
            reproGrad(:,:), &                !Residual Projection for Discont Capturing
            reproSGrad(:,:),&                !Gradient Projection
            exsigma(:,:,:),&
            exveloc(:,:,:),&
            express(:,:), &
            sigmaold(:,:), &
            psiReal(:,:),&
            sigma_cp(:,:)
            
            
      ! Oss Split      
      real(rp), allocatable :: &
            reproSDiv(:,:), &     
            reproUGradU(:,:), & 
            reproGradP(:,:), &
            reproLapla(:,:), &
            reproDivU(:,:), &
            reproGradU(:,:), &
            reproUGradS(:,:), &
            reproSGradU(:,:), &
            reproExpS(:,:)
            
            
      real(rp) :: &
            ersh1t(2), &
            erslit(2), &
            !Temperature model
            ReferenceTemp, &
            c1_WLF, &
            c2_WLF, &
            nu0_LCR, &
            alpha_Arrhenius
      
      !Statistics  
      real(rp) :: &      
            tasigmin,&                             ! Minimum tau
            tasigmax,&                             ! Maximum tau
            tasigmea                               ! Mean tau
            
      real(rp)              ::      resis,incremental   !Residual of outer iteration
      integer(ip)  :: ntens, &
            kfl_bc_number, &
            kfl_splitOSSMomentum,&     !Flag to turn on or turn off the cross-stabilizing terms for momentum equation
            kfl_splitOSSConstitutive,& !Flag to turn on or turn off the cross-stabilizing terms for constitutive equation
            kfl_reproBoundzero, &      !Flag to turn on or turn off the Residual projection to 0 in boundaries
            Giesekus_model,&           !Giesekus model in  use
            PTT_model,&                !PTT model in use
            kfl_cotem_WLF,&            !WLF_model, temperature coupling
            kfl_cotem_Arrhenius,&      !Arrhenius model, temperature coupling   
            kfl_cotem_Boussinesq,&     !Boussinesq model, temperature coupling
            LogFormulation,&           !Log conformation reformulation
            kfl_LogFormulation,&                  !Flag to indicates the existence of proporcionality between conformation tensor and identity tensor
            kfl_linearConstitutiveTerms,&        !Flag to choose linearization to non-linear terms. 
            kfl_linearConvectiveTerm,& !Flag to choose linearization to the convective term. 
            kfl_ntausmooth,&
            kfl_reproTemporalTermZero,&!Flag to eliminate the temporal term in stabilization if the orthogonal projection is chosen.
            kfl_linconvec              !Flag to choose linearization to the convective term 
            
      type(r3p), allocatable ::vesgs3(:) !Velocity subgrid scales needed in dyn Split-Oss 
      type(r3p), allocatable ::sisgs(:), sisgs2(:), sisgs3(:)  !Stress subgrid scales
            
      type(r1p), allocatable :: lambdarray(:)    !Gauss point lambda (relaxation parameter)
      type(r1p), allocatable :: alpha3array(:)   !Third stabilization parameter
      type(r2p), allocatable :: Term(:)
      real(rp), allocatable  :: taudet(:) 
      type(r1p), allocatable :: tau_mom(:)
      type(r1p), allocatable :: tau_sig(:)

      
      
contains

      procedure :: SetExmod             => sup_SetExmod
      procedure :: SetNdofn             => sup_SetNdofn
      procedure :: SetNdofbc            => sup_SetNdofbc      
      procedure :: SpecificIniunk       => sup_iniunk
      procedure :: SpecificBegite       => sup_begite
      procedure :: SpecificEndite       => sup_endite
      procedure :: Cvgunk               => sup_cvgunk
      procedure :: Memall               => sup_memall
      procedure :: Bouope               => sup_bouop0      
      procedure :: SpecificBegste       => sup_begste  
      procedure :: SpecificEndste       => sup_endste
      procedure :: InnerResiduals       => sup_InnerResiduals
      procedure :: Elmope               => supm_elmope
      procedure :: EnditeElmope         => supm_EnditeElmope
      procedure :: SpecificTurnof       => sup_turnof
      procedure :: SpecificReaphy       => sup_reaphy 
      procedure :: Ifconf               => sup_ifconf  
      procedure :: SpecificExaerr       => sup_exaerr
      procedure :: NstincSpecificReampi => sup_reampi
      procedure :: SpecificGetste       => sup_getste
      procedure :: SpecificReabcs       => sup_reabcs
      procedure :: SpecificReanut       => sup_reanut
      procedure :: SpecificUpdbcs       => sup_updbcs
      procedure :: IncrementalLambda    => sup_incremental
      procedure :: GetSigmaArray
      !Boundary operations 
      procedure :: EndBouope                 => sup_EndBouope 
      procedure :: SpecificReadOnNodes       => sup_ReadOnNodes
      procedure :: PointTracking             => sup_outtpo
      
      !Statistics
      procedure :: InitStats                 => supm_InitStats
      procedure :: InGaussStats_sup          => supm_InGaussStats_sup
      procedure :: FinalizeStats             => supm_FinalizeStats

      procedure :: NstincSpecificOutput      => sup_output
      procedure :: NstincSpecificRestart     => sup_restar
      procedure :: GaussRestart              => sup_gaussRestart

   end type

  interface ThreeFieldNSProblem_Const
      procedure constructor
  end interface ThreeFieldNSProblem_Const
   
   interface

      subroutine sup_outtpo(a)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine
      
      subroutine sup_SetExmod(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine 

      subroutine sup_SetNdofn(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine 
      
      subroutine sup_SetNdofbc(a)
         import
         implicit none         
         class(ThreeFieldNSProblem) :: a      
      end subroutine      
         
      subroutine sup_iniunk(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine 

      subroutine sup_begite(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine 
      
      subroutine sup_endite(a,itask)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip)                :: itask
      end subroutine  

      subroutine sup_cvgunk(a,itask)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip), intent(in)     :: itask
      end subroutine  
      
      subroutine sup_memall(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine 
      
      subroutine sup_getste(a,dtinv)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         real(rp) :: dtinv      
      end subroutine      

      subroutine sup_output(a,itask)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip)                :: itask
      end subroutine 

      subroutine sup_begste(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine 

      subroutine sup_exaerr(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine       
      
      subroutine sup_endste(a,itask)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip)                :: itask
      end subroutine 

      subroutine sup_restar(a,itask)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine      

      subroutine sup_gaussRestart(a,itask)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         integer(ip), intent(in) :: itask
      end subroutine      
      
      subroutine sup_InnerResiduals(a,rinsi,rprnsi) !,rsrnsi)
         use typre
         import
         implicit none
         class(ThreeFieldNSProblem) :: a
         real(rp) :: rinsi,rprnsi !,rsrnsi   
      end subroutine  

      subroutine supm_elmope(a)
         import
         implicit none
         class(ThreeFieldNSProblem) :: a      
      end subroutine
     
     subroutine supm_EnditeElmope(a,task)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a   
        character(6) :: task
     end subroutine   

     subroutine sup_ReadOnNodes(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a      
     end subroutine    
      
     subroutine sup_turnof(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
     end subroutine  
      
     subroutine sup_updbcs(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
     end subroutine 
     
     subroutine sup_reabcs(a,itask,kflag)
        use typre
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
        integer(ip) :: itask
        integer(ip), optional :: kflag
     end subroutine 
     
     subroutine sup_reanut(a,itask)
        use typre
        import
        implicit none         
        integer(ip) :: itask
        class(ThreeFieldNSProblem) :: a      
     end subroutine      
     
     subroutine sup_reaphy(a,itask)
        use typre
        import
        implicit none
        integer(ip) :: itask
        class(ThreeFieldNSProblem) :: a
     end subroutine  
     
     subroutine sup_ifconf(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
     
     end subroutine
     
      subroutine sup_reampi(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
     
     end subroutine
           
     subroutine sup_bouop0(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
     
     end subroutine 
     
     subroutine sup_EndBouope(NSProblem,task)
        import
        implicit none
        class(ThreeFieldNSProblem),target :: NSProblem
        character(6) :: task
     end subroutine 
     
     subroutine supm_InitStats(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a      
     end subroutine
     
     subroutine supm_InGaussStats_sup(a,acden,acvis,gpvno,chale,timom,tisig)
        use typre
        import
        implicit none
        class(ThreeFieldNSProblem) :: a
        real(rp) :: acden,acvis,gpvno,chale(2),timom,tisig
     end subroutine
     
     subroutine supm_FinalizeStats(a)
        import
        implicit none
        class(ThreeFieldNSProblem) :: a      
     end subroutine
      
   end interface
 
contains 

    function constructor()
        class(ThreeFieldNSProblem), pointer :: constructor

        allocate(constructor)

    end function constructor
    
    subroutine sigma_SetCtime(a,ctime)
      use typre
      implicit none
      class(ThreeFieldNSProblem) :: a
      real(rp) :: ctime

      a%ctime = ctime
   end subroutine

    
   subroutine GetSigmaArray(a,sigma)
      use typre
      implicit none
      class(ThreeFieldNSProblem), target :: a
      real(rp), pointer :: sigma(:,:) 
      if(a%LogFormulation==0)  sigma => a%sigma(:,:,1)
      if(a%LogFormulation==1)  sigma => a%sigmaold(:,:)
      
   end subroutine
   
   
    subroutine sup_incremental(a,imat,lambdaIn, lambda)
      use typre
      implicit none
      class(ThreeFieldNSProblem) :: a
      integer(ip), intent(in) :: imat
      real(rp), intent(in)    :: lambdaIn
      real(rp), intent(inout) :: lambda 
      
          
      
      if(a%kfl_timei==0) then
         if(a%itera<=a%incremental)then
            lambda= (lambdaIn/a%incremental)*a%itera
         else
            lambda = lambdaIn
         end if
            
       else if(a%kfl_timei==1)then
         if(a%istep<=a%incremental)then
            lambda= (lambdaIn/a%incremental)*a%istep
         else 
            lambda = lambdaIn
         end if         
      end if       

   end subroutine 
   
 

 
end module
