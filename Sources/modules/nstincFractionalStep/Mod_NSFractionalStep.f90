module Mod_NSFractionalStep
   use typre
   use Mod_Mesh
   use Mod_Memor
   
   use Mod_NavierStokes
   use Mod_ParallelSystemInterface
   private
   public NSFractionalStepProblem,NSFractionalStepProblem_Const
   
   type, extends(NavierStokesProblem) :: NSFractionalStepProblem
      
      !Linear System
      class(ParallelSystemInterface), pointer   :: LinearSystemP => NULL()
      real(rp), allocatable                     :: unknoP(:,:)
      integer(ip),  allocatable :: kfl_fixpr(:)       ! Nodal fixity for the pressure (NSF)
      real(rp), allocatable :: bvpress(:)
      
      integer(ip) :: kfl_restar = 0_ip                ! I need to know if this was a restart problem
      !to decide wether to initialize hydrostatic pressure field or not
      integer(ip) :: FractionalExternalIterations = 1
      
contains   
   
      procedure :: SetExmod                           => nsf_SetExmod
      procedure :: LinearSystemMemall                 => nsf_LinearSystemMemall
      procedure :: LinearSystemTurnof                 => nsf_LinearSystemTurnof
      procedure :: SpecificReanut                     => nsf_reanut
      procedure :: SpecificReaMPI                     => nsf_reampi
      procedure :: SpecificCrankNicolsonEndste        => NULLSUB
      !procedure :: SpecificReabcsMPI                  => nsf_reabcsMPI
      procedure :: Solite                             => nsf_solite
      procedure :: Endite                             => nsf_endite
      procedure :: SpecificTurnof                     => nsf_turnof
      procedure :: InnerResiduals                     => nsf_InnerResiduals
      procedure ::  SpecificReabcs                    => nsf_Reabcs
      
      !For Initializing the hydrostatic pressure field
      !Need to note down if this is a initial restart problem
      procedure :: SpecificRestart                    => nsi_restar
      procedure :: SpecificBegste                     => nsf_begste
      procedure :: SpecificOuterr                     => nsf_outerr
      procedure :: SpecificRefine                     => nsf_refine
   end type
   
   interface
      subroutine nsf_solite(a)
         import NSFractionalStepProblem
         implicit none
         class(NSFractionalStepProblem) :: a
   
      end subroutine  
      
      subroutine nsf_endite(a,itask)
         use typre
         import NSFractionalStepProblem
         implicit none
         class(NSFractionalStepProblem) :: a
         integer(ip)                    :: itask
      end subroutine  
      
      
      subroutine nsf_InnerResiduals(a,rinsi,rprnsi)
         use typre
         import NSFractionalStepProblem
         implicit none
         class(NSFractionalStepProblem) :: a
         real(rp) :: rinsi,rprnsi
      end subroutine  
      
      subroutine nsi_restar(a,itask)
         use typre
         import NSFractionalStepProblem
         implicit none
         class(NSFractionalStepProblem) :: a
         integer(ip), intent(in)                   :: itask
      end subroutine  
      
      subroutine nsf_begste(a)
         import NSFractionalStepProblem
         implicit none
         class(NSFractionalStepProblem) :: a
      end subroutine  
      
      subroutine nsf_refine(a,itask)
         use typre
         import NSFractionalStepProblem
         implicit none
         class(NSFractionalStepProblem) :: a
         character(6) :: itask
      end subroutine
      
   end interface
   
  interface NSFractionalStepProblem_Const
      procedure constructor
  end interface NSFractionalStepProblem_Const

  contains

    function constructor()
        class(NSFractionalStepProblem), pointer :: constructor

        allocate(constructor)

    end function constructor

   subroutine nsf_LinearSystemMemall(a)
      implicit none
      class(NSFractionalStepProblem) :: a
      
      integer(ip) :: npoin,ndime
      character(150) :: exmod,exfile
      character(150) :: auxstring
      
      
   
      exfile = 'nsi'
      
      !Velocity linear system
      exmod = 'nsf'
      call a%Mesh%GetNdime(ndime)
      auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exfile))//'.sol '
      
      !Create ParallelSystem
      call a%ParallelLibrary%CreateSystem(a%LinearSystem, a%Memor) 
      call a%LinearSystem%SetFlush(a%kfl_flush)
      call a%LinearSystem%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Mesh%InitializeSystem(ndime,0,a%LinearSystem,a%Memor,auxstring,exmod,a%lun_solve)
   
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(ndime,npoin,a%unkno,'unkno_nsf',trim(exmod)//'_LinearSystemMemall')

      !Pressure Linear System
      exmod = 'nsp'
      auxstring = trim(a%InputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exfile))//'.sol '
      call a%ParallelLibrary%CreateSystem(a%LinearSystemP, a%Memor) 
      call a%LinearSystemP%SetFlush(a%kfl_flush)
      call a%LinearSystemP%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
      call a%Mesh%InitializeSystem(1,0,a%LinearSystemP,a%Memor,auxstring,exmod,a%lun_solve)
   
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(1,npoin,a%unknoP,'unknoP',trim(exmod)//'_LinearSystemMemall')
      
      
   end subroutine
   
   subroutine nsf_LinearSystemTurnof(a)
      use typre
      use Mod_PhysicalProblem
      implicit none
      class(NSFractionalStepProblem) :: a
   
      integer(ip) :: npoin,ndime
      
      !Velocity Linear System
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Memor%dealloc(ndime,size(a%unkno,2),a%unkno,'unkno_nsf','nsf_LinearSystemTurnof')
      call a%LinearSystem%Deallocate
      call a%ParallelLibrary%DeallocateSystem(a%LinearSystem, a%Memor) 
      !call a%Mesh%DeallocateSystem(a%LinearSystem,a%Memor)
      
      !Pressure Linear System
      call a%Memor%dealloc(1,size(a%unknoP,2),a%unknoP,'unknoP','nsf_LinearSystemTurnof')
      call a%LinearSystemP%Deallocate
      call a%ParallelLibrary%DeallocateSystem(a%LinearSystemP, a%Memor) 
      !call a%Mesh%DeallocateSystem(a%LinearSystemP,a%Memor)
      
   
end subroutine
   
   subroutine nsf_SetExmod(a)
      use typre
      implicit none
      
      class(NSFractionalStepProblem) :: a
      
      a%exmod = 'nsi'
      a%namod = 'NSTINC_FractionalStep'
      
   end subroutine
   
   subroutine nsf_reanut(a,itask)
      use typre
      use MPI
      implicit none
      class(NSFractionalStepProblem) :: a
      integer(ip) :: itask
      interface
         subroutine nsi_reanut(a,itask)
            use typre
            use Mod_NavierStokes
            implicit none
            class(NavierStokesProblem) :: a
            integer(ip) :: itask
         end subroutine   
      end interface
      
      !For Listener
      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
      integer(ip) :: ierr
      
      call nsi_reanut(a,itask)
      if (itask == 0) then
         a%kfl_repro_SkipFE = 1    !Default is skip Finite Element part of the residual
         a%FractionalExternalIterations = 1
         
      elseif (itask == 1) then
         call a%Listener%getarrs(words,param,nnpar,nnwor)
         if(words(1)=='EXTIT') then
            a%FractionalExternalIterations = param(1)
         endif
         
      elseif (itask == 100) then
         if (a%kfl_repro == 0) 	a%kfl_repro = 1			
      endif
   end subroutine
   
   subroutine nsf_reampi(a)
      use MPI
      implicit none
      class(NSFractionalStepProblem) :: a
      interface
         subroutine nsi_reampi(a)
            use typre
            use Mod_NavierStokes
            implicit none
            class(NavierStokesProblem) :: a
         end subroutine   
      end interface
      
      integer(ip) :: ierr
      
      call nsi_reampi(a)
      CALL MPI_BCAST(a%FractionalExternalIterations, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)  
   end subroutine
   
   
   subroutine nsf_outerr(a)
      use typre
      implicit none
      class(NSFractionalStepProblem) :: a
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
      if(a%kfl_inist==1)then
         call runend('INISTOKES NOT READY IN FRACTIONAL STEP')
      end if
   end subroutine
      
   subroutine nsf_Reabcs(a,itask,kflag)
      implicit none
      class(NSFractionalStepProblem) :: a
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
         
         subroutine nsf_pbcpre(a)
            use typre
            import NSFractionalStepProblem
            implicit none
            class(NSFractionalStepProblem) :: a
         end subroutine  
      end interface
      
      
      call nsi_Reabcs(a,itask,kflag)
      if (itask == 100) then
         call a%Mesh%GetNpoin(npoin)
         call a%Memor%alloc(npoin,a%kfl_fixpr,'kfl_fixpr','nsf_reabcsMPI')
         call a%Memor%alloc(npoin,a%bvpress,'bvpress','nsf_reabcsMPI')
         call nsf_pbcpre(a)
      endif
      
   end subroutine
   
   
   subroutine nsf_turnof(a)
      implicit none
      class (NSFractionalStepProblem) :: a
      
      integer(ip) :: npoin
      interface
         subroutine nsi_turnof(a)
            use typre
            use Mod_NavierStokes
            implicit none
            class(NavierStokesProblem) :: a
         end subroutine   
      end interface
      
      call nsi_turnof(a)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(size(a%kfl_fixpr),a%kfl_fixpr,'kfl_fixpr','nsf_reabcsMPI')
      call a%Memor%dealloc(size(a%bvpress),a%bvpress,'bvpress','nsf_reabcsMPI')
      
   end subroutine
   
   subroutine NULLSUB(a)
      use typre
      implicit none
      class(NSFractionalStepProblem) :: a
   end subroutine



end module
