module Mod_SUPFractionalStep
   use typre
   use Mod_Mesh
   use Mod_Memor   
   use Mod_ThreeField
   use Mod_ParallelSystemInterface
   private
   public SUPFractionalStepProblem,SUPFractionalStepProblem_Const   
   type, extends(ThreeFieldNSProblem) :: SUPFractionalStepProblem
 
      !Linear System
      class(ParallelSystemInterface), pointer   :: LinearSystemC => NULL()
      class(ParallelSystemInterface), pointer   :: LinearSystemS => NULL()
      class(ParallelSystemInterface), pointer   :: LinearSystemST => NULL()
      real(rp), allocatable                     :: unknoC(:,:),unknoS(:,:),unknoST(:,:)
      integer(ip) :: &
         kfl_goiteS, &          !Keep iterating constituive equation
         iteraS, &              !Internal iteration counter constituive equation 
         kfl_goiteY, &          !Keep iterating Yosida momentum equation        
         iteraY                 !Internal iteration counter Yosida momentum
      real(rp) :: kdicCap
      
      integer(ip),  allocatable ::&
      kfl_fixpr(:)                    ! Nodal fixity for the pressure (NSF)
      real(rp), allocatable :: &
      bvpress(:)
 
contains   


      procedure :: SetExmod           => supf_SetExmod
      procedure :: LinearSystemMemall => supf_LinearSystemMemall
      procedure :: LinearSystemTurnof => supf_LinearSystemTurnof
      procedure :: SpecificCrankNicolsonEndste => NULLSUB
      procedure :: Solite               => supf_solite
      procedure :: Endite               => supf_endite
      procedure :: SpecificBegite       => supf_begite     
      procedure :: SpecificTurnof       => supf_turnof       
      procedure ::  SpecificReabcs      => supf_Reabcs
      procedure :: SpecificBegste => supf_begste
      procedure :: SpecificOuterr         => supf_outerr
      
    end type
   
    interface
    
      subroutine supf_solite(a)
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
   
      end subroutine  
      
      subroutine supf_endite(a,itask)
         use typre
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
         integer(ip)                    :: itask
      end subroutine  
      
      subroutine supf_begite(a)
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a      
      end subroutine     
       
      
      subroutine supf_begste(a)
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
      end subroutine  
      
   end interface
   

  interface SUPFractionalStepProblem_Const
      procedure constructor
  end interface SUPFractionalStepProblem_Const

  contains

      function constructor()
          class(SUPFractionalStepProblem), pointer :: constructor

          allocate(constructor)

      end function constructor


   subroutine supf_LinearSystemMemall(a)
      use typre   
      implicit none
      class(SUPFractionalStepProblem) :: a
      
      integer(ip) :: npoin,ndime,ntens
      character(150) :: exmod,exfile
      character(150) :: auxstring      
   
      exfile = 'nsi'
      call a%Mesh%GetNdime(ndime)      
      ntens=(ndime-1)*(ndime-1)+2
      
      !Velocity linear system
      exmod = 'nsf'      
      auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exfile))//'.sol '
      call a%ParallelLibrary%CreateSystem(a%LinearSystem, a%Memor) 
      call a%LinearSystem%SetFlush(a%kfl_flush)
      call a%LinearSystem%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Mesh%InitializeSystem(ndime,0,a%LinearSystem,a%Memor,auxstring,exmod,a%lun_solve)
   
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(ndime,npoin,a%unkno,'unkno',trim(exmod)//'_LinearSystemMemall')

      !Pressure Linear System
      exmod = 'nsf'
      auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exfile))//'.sol '
      call a%ParallelLibrary%CreateSystem(a%LinearSystemC, a%Memor) 
      call a%LinearSystemC%SetFlush(a%kfl_flush)
      call a%LinearSystemC%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Mesh%InitializeSystem(1,0,a%LinearSystemC,a%Memor,auxstring,exmod,a%lun_solve)   
      call a%Memor%alloc(1,npoin,a%unknoC,'unknoC',trim(exmod)//'_LinearSystemMemall')
      
      !Stress Linear System
      exmod = 'nsf'
      auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exfile))//'.sol '
      call a%ParallelLibrary%CreateSystem(a%LinearSystemS, a%Memor) 
      call a%LinearSystemS%SetFlush(a%kfl_flush)
      call a%LinearSystemS%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Mesh%InitializeSystem(ntens,0,a%LinearSystemS,a%Memor,auxstring,exmod,a%lun_solve)
      call a%Memor%alloc(ntens,npoin,a%unknoS,'unknoS',trim(exmod)//'_LinearSystemMemall')     
      
      if(a%kfl_inist==1)then
         !Stokes Aditional Problem
         exmod = 'nsf'
         auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exfile))//'.sol '
         call a%ParallelLibrary%CreateSystem(a%LinearSystemST, a%Memor) 
         call a%LinearSystemST%SetFlush(a%kfl_flush)
         call a%LinearSystemST%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
         call a%Mesh%InitializeSystem(ndime+1,0,a%LinearSystemST,a%Memor,auxstring,exmod,a%lun_solve)
         call a%Memor%alloc(ndime+1,npoin,a%unknoST,'unknoST',trim(exmod)//'_LinearSystemMemall')      
      end if
      
   end subroutine
   
   subroutine supf_LinearSystemTurnof(a)
      use typre
      use Mod_PhysicalProblem
      implicit none
      class(SUPFractionalStepProblem) :: a
      
         integer(ip) :: npoin,ndime,ntens
         
         
         call a%Mesh%GetNpoin(npoin)
         call a%Mesh%GetNdime(ndime)         
         ntens=(ndime-1)*(ndime-1)+2
         
         !Velocity Linear System
         call a%Memor%dealloc(ndime,npoin,a%unkno,'unkno','supf_LinearSystemTurnof')
         call a%LinearSystem%Deallocate
         call a%ParallelLibrary%DeallocateSystem(a%LinearSystem, a%Memor) 
         !call a%Mesh%DeallocateSystem(a%LinearSystem,a%Memor)
         
         !Pressure Linear System
         call a%Memor%dealloc(1,npoin,a%unknoC,'unknoC','supf_LinearSystemTurnof')
         call a%LinearSystemC%Deallocate
         call a%ParallelLibrary%DeallocateSystem(a%LinearSystemC, a%Memor) 
         !call a%Mesh%DeallocateSystem(a%LinearSystemC,a%Memor)  
         
         !Stress Linear System
         call a%Memor%dealloc(ntens,npoin,a%unknoS,'unknoS','supf_LinearSystemTurnof')
         call a%LinearSystemS%Deallocate
         call a%ParallelLibrary%DeallocateSystem(a%LinearSystemS, a%Memor) 
         !call a%Mesh%DeallocateSystem(a%LinearSystemS,a%Memor)  
         
         if(a%kfl_inist==1)then
            !Stokes aditional problem Linear System
            call a%Memor%dealloc(ndime+1,npoin,a%unknoST,'unknoST','supf_LinearSystemTurnof')
            call a%LinearSystemST%Deallocate
            call a%ParallelLibrary%DeallocateSystem(a%LinearSystemST, a%Memor) 
            !call a%Mesh%DeallocateSystem(a%LinearSystemST,a%Memor)           
         end if
         
   end subroutine
   
   subroutine supf_SetExmod(a)
      use typre
      implicit none
      
      class(SUPFractionalStepProblem) :: a
      
      a%exmod = 'nsi'
      a%namod = 'SUPFractional'
      
   end subroutine


   subroutine supf_outerr(a)
      use typre
      implicit none
      class(SUPFractionalStepProblem) :: a
      integer(ip) :: itask
      
      interface
         subroutine nsi_outerr(a)
            use typre
            use Mod_NavierStokes
            implicit none
            class(NavierStokesProblem) :: a
         end subroutine   
      end interface
      
      
      call nsi_outerr(a)
      if(a%kfl_timei==0)then
         call runend('FRACTIONAL STEP CAN NOT BE USED IN NON TRANSIENT CASES')
      end if

   end subroutine
   
      
   subroutine supf_Reabcs(a,itask,kflag)
      implicit none
      class(SUPFractionalStepProblem) :: a
      integer(ip) :: itask
      integer(ip), optional :: kflag
      
      integer(ip) :: npoin
      interface
         subroutine nsi_Reabcs(a,itask,kflag)
            use typre
            use Mod_NavierStokes
            implicit none
            class(NavierStokesProblem) :: a
            integer(ip) :: itask
            integer(ip), optional :: kflag
         end subroutine  
         
         subroutine supf_pbcpre(a)
            use typre
            import SUPFractionalStepProblem
            implicit none
            class(SUPFractionalStepProblem) :: a
         end subroutine  
      end interface
      
      
      call nsi_Reabcs(a,itask,kflag)
      if (itask == 100) then
         call a%Mesh%GetNpoin(npoin)
         call a%Memor%alloc(npoin,a%kfl_fixpr,'kfl_fixpr','nsf_reabcsMPI')
         call a%Memor%alloc(npoin,a%bvpress,'bvpress','nsf_reabcsMPI')
         call supf_pbcpre(a)
      endif
      
   end subroutine
   
   
   subroutine supf_turnof(a)
      implicit none
      class (SUPFractionalStepProblem) :: a
      
      integer(ip) :: npoin
      interface
         subroutine sup_turnof(a)
            use typre
            use Mod_ThreeField
            implicit none
            class(ThreeFieldNSProblem) :: a
         end subroutine   
      end interface
      
      call sup_turnof(a)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(npoin,a%kfl_fixpr,'kfl_fixpr','nsf_reabcsMPI')
      call a%Memor%dealloc(npoin,a%bvpress,'bvpress','nsf_reabcsMPI')
      
   end subroutine
   
   subroutine NULLSUB(a)
      use typre
      implicit none
      class(SUPFractionalStepProblem) :: a
   end subroutine


end module
