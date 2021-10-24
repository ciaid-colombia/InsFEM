module Mod_NSCompressibleExplicit
   use typre
   use Mod_Listen
   use Mod_NSCompressible
   implicit none
   private
   public NSCompressibleExplicitProblem, NSCompressibleExplicitProblem_Const
   
   type, extends(NSCompressibleProblem) :: NSCompressibleExplicitProblem
      
      real(rp), allocatable :: &
        femti(:,:)                             !Finite element scale stage evaluation

! Subscales indices for the explicit time integration scheme:
!position 1: stage increment
!position 2: stage contribution
!position 3: previous time step 

      type(r2p), allocatable :: &
            sgsti(:)                            !Stage evaluation for the subscales explicit
                                                !time integration scheme   
      character(5) :: &
         kfl_tsche_ex                           !Temporal Scheme: 1st order time derivative in eqns

contains

      procedure :: Cvgunk            => nsc_cvgunk_ex
      procedure :: Elmope            => nsc_NULLSUB
      procedure :: Bouope            => nsc_NULLSUB

      procedure :: Solite          => nsc_ExplicitSolite 
      procedure :: EndElmope      => nsc_EndElmope_ex 
      procedure :: DynamicSGSUpdate  => nsc_SGSTimeUpdate

      procedure :: SpecificNSCompReaphy  => nsc_NULLSUBitask
      procedure :: SpecificNSCompReanut  => nsc_reanut_ex
      procedure :: SpecificNSCompReaMPI  => nsc_reaMPI_ex
      procedure :: SpecificNSCompReabcs  => nsc_reabcs_ex
      procedure :: SpecificNSCompMemall  => nsc_memall_ex
      procedure :: SpecificNSCompTurnof  => nsc_turnof_ex
      procedure :: SpecificNSCompBegste  => nsc_NULLSUB
      procedure :: SpecificNSCompBegite  => nsc_NULLSUB
      procedure :: SpecificNSCompEndite  => nsc_NULLSUBitask
      procedure :: SpecificNSCompEndste  => nsc_endste_ex
      procedure :: SpecificNSCompExaerr  => nsc_exaerr_ex
      procedure :: SpecificNSCompRefine  => nsc_NULLSUBcharitask
      procedure :: SpecificNSCompRefinementCriteria  => nsc_NULLSUBmarkel
      procedure :: SpecificNSCompRestar    => nsc_restar_ex

   end type

   interface
     
      subroutine nsc_cvgunk_ex(a,itask)
         use typre
         import NSCompressibleExplicitProblem
         implicit none
         class(NSCompressibleExplicitProblem) :: a      
         integer(ip), intent(in) :: itask
      
      
      end subroutine
      
      subroutine nsc_EndElmope_ex(NSCompExplicitProblem,task)
         import NSCompressibleExplicitProblem
         implicit none
         class(NSCompressibleExplicitProblem), target :: NSCompExplicitProblem
         character(6) :: task
      end subroutine

   end interface
    
  interface NSCompressibleExplicitProblem_Const
      procedure constructor
  end interface NSCompressibleExplicitProblem_Const
   
contains

function constructor()
    class(NSCompressibleExplicitProblem), pointer :: constructor

    allocate(constructor)

end function constructor

      subroutine nsc_reanut_ex(a,itask)
         use typre
         implicit none
   
         class(NSCompressibleExplicitProblem) :: a      
         integer(ip) :: itask
         
         !For Listener
         real(rp), pointer     :: param(:) => NULL()
         character(5), pointer :: words(:) => NULL()
         integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
         
         if (itask == 0) then
            a%kfl_tsche_ex = 'RK4'                             ! Runge kutta 4 explicit scheme
            
         elseif (itask == 1) then
            call a%Listener%getarrs(words,param,nnpar,nnwor)
            if(words(1)=='EXPLI') then
            a%kfl_tsche_ex = words(2)
            end if

         elseif (itask == 100) then
            
         endif

      end subroutine
      
      subroutine nsc_reaMPI_ex(a)
         use MPI
         implicit none
         class(NSCompressibleExplicitProblem) :: a      

         integer(ip) :: ierr
         
         CALL MPI_BCAST(a%kfl_tsche_ex,len(a%kfl_tsche_ex),MPI_CHAR, a%MPIroot, a%MPIcomm, ierr)  
      end subroutine

      subroutine nsc_reabcs_ex(a,itask)
         use typre
         use MPI
         use Mod_NSCompressible
         use Mod_Mesh
         use Mod_NscExacso
         implicit none
   
         class(NSCompressibleExplicitProblem) :: a      
         integer(ip) :: itask

         type(NscExacso) :: exacso     
         !Exact Values
         real(rp) :: exden, extem
         real(rp), allocatable   :: exdeg(:),exvel(:),exveg(:,:),exteg(:)
         !Nodal coordinates
         real(rp), pointer       :: coord(:) => NULL()   
         real(rp), pointer :: exnor(:,:) => NULL()
         integer(ip) :: ndime,ipoin,npoin,ibopo
         
         if (itask == 0) then
            
         elseif (itask == 1) then

         elseif (itask == 100) then
            
            !Exact solution boundary conditions
            if(a%kfl_exacs/=0) then
               ! Allocate exact components
               call a%Mesh%GetNdime(ndime)
               call a%Mesh%Getnpoin(npoin)
               call a%Memor%alloc(ndime,exdeg,'exdeg','nsc_reabcs_ex')   
               call a%Memor%alloc(ndime,exvel,'exvel','nsc_reabcs_ex')   
               call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsc_reabcs_ex')     
               call a%Memor%alloc(ndime,exteg,'exteg','nsc_reabcs_ex')   
               a%kfl_fixbo = -1
               a%kfl_fixno = -1
               a%bvess = 0
               
               do ipoin = 1,npoin
                  call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
                  if (ibopo /= 0) then
                     call a%Mesh%GetPointCoord(ipoin,coord)
                     call exacso%nsc_ComputeSolution(ndime,coord,a)
                     call exacso%nsc_GetDensity(ndime,exden,exdeg)           
                     call exacso%nsc_GetVelocity(ndime,exvel,exveg)           
                     call exacso%nsc_GetTemperature(ndime,extem,exteg)           
                     a%kfl_fixno(:,ipoin) = 1_ip
                     a%bvess(1,ipoin,1) = exden
                     a%bvess(2,ipoin,1) = exvel(1)
                     a%bvess(3,ipoin,1) = exvel(2)
                     a%bvess(4,ipoin,1) = extem
                  endif
               enddo   
               ! Deallocate exact components
               call a%Memor%dealloc(ndime,exdeg,'exdeg','nsc_reabcs_ex')   
               call a%Memor%dealloc(ndime,exvel,'exvel','nsc_reabcs_ex')   
               call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsc_reabcs_ex')     
               call a%Memor%dealloc(ndime,exteg,'exteg','nsc_reabcs_ex')   
      
            endif
         endif

      end subroutine

      subroutine nsc_memall_ex(a)
         use typre
         use Mod_Mesh
         use Mod_Memor
         implicit none
         class(NSCompressibleExplicitProblem) :: a      

         integer(ip) :: ncomp,ndime,npoin,nelem,ielem,pnode,&
                  pgaus,sgsti_coun
                  
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         call a%Mesh%GetNelem(nelem)
         
         a%ncomp = 3
         ncomp = a%ncomp
         
         !Unknowns
         call a%Memor%alloc(npoin,ncomp,a%densf,'densf','nsc_memall_ex')
         call a%Memor%alloc(ndime,npoin,ncomp,a%momen,'momen','nsc_memall_ex')
         call a%Memor%alloc(npoin,ncomp,a%energ,'energ','nsc_memall_ex')
         !Explicit Time dependent finite element scales
         call a%Memor%alloc(a%ndofn,npoin,a%femti,'femti','nsc_memall_ex')

         !Primitives
         call a%Memor%alloc(ndime,npoin,ncomp,a%veloc,'veloc','nsc_memall_ex')
         call a%Memor%alloc(npoin,ncomp,a%press,'press','nsc_memall_ex')
         call a%Memor%alloc(npoin,ncomp,a%tempe,'tempe','nsc_memall_ex')

         !Explicit Time dependent subgrid scales
         if(a%kfl_tacsg ==1) then
            sgsti_coun = 0
            
            call a%Memor%alloc(nelem,a%sgsti,'sgsti','nsc_memall_ex')
      
            do ielem=1,nelem
               call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
               allocate(a%sgsti(ielem)%a(a%ndofn,pgaus))
               a%sgsti(ielem)%a = 0.0_rp
               sgsti_coun = sgsti_coun + a%ndofn*pgaus
            end do
      
            call a%Memor%allocObj(0,'sgsti%a','nsc_memall_ex',sgsti_coun*rp)
      
         end if

      end subroutine

      subroutine nsc_turnof_ex(a)
         use typre
         use Mod_Mesh
         use Mod_Memor
         implicit none
         class(NSCompressibleExplicitProblem) :: a      

         integer(ip) :: ncomp,ndime,npoin,nelem,pgaus,ielem,pnode,sgsti_coun

         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         call a%Mesh%GetNelem(nelem)
         
         ncomp = a%ncomp

         !Unknowns
         call a%Memor%dealloc(npoin,ncomp,a%densf,'densf','nsc_turnof_ex')
         call a%Memor%dealloc(ndime,npoin,ncomp,a%momen,'momen','nsc_turnof_ex')
         call a%Memor%dealloc(npoin,ncomp,a%energ,'energ','nsc_turnof_ex')

         !Primitives
         call a%Memor%dealloc(ndime,npoin,ncomp,a%veloc,'veloc','nsc_turnof_ex')
         call a%Memor%dealloc(npoin,ncomp,a%press,'press','nsc_turnof_ex')
         call a%Memor%dealloc(npoin,ncomp,a%tempe,'tempe','nsc_turnof_ex')


         !Explicit Time dependent finite element scales
         call a%Memor%dealloc(a%ndofn,npoin,a%femti,'femti','nsc_turnof_ex')

         !Explicit Time dependent subgrid scales
         if(a%kfl_tacsg ==1) then
            sgsti_coun = 0
         
            do ielem=1,nelem
               call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
               deallocate(a%sgsti(ielem)%a)
               sgsti_coun = sgsti_coun + a%ndofn*pgaus
            end do
         
            call a%Memor%deallocObj(0,'sgsti%a','nsc_turnof_ex',sgsti_coun*rp)
            call a%Memor%dealloc(nelem,a%sgsti,'sgsti','nsc_turnof_ex')
         
         end if

      end subroutine

      subroutine nsc_ExplicitSolite(a)
         implicit none
         class(NSCompressibleExplicitProblem) :: a      
      
         interface
            subroutine nsc_solite_ex(a)
               import NSCompressibleExplicitProblem
               implicit none
               class(NSCompressibleExplicitProblem) :: a
            end subroutine
            subroutine nsc_elmdir_ex(a)
               import NSCompressibleExplicitProblem
               implicit none
               class(NSCompressibleExplicitProblem) :: a
            end subroutine
            subroutine nsc_GhostCommunicate(a)
               import NSCompressibleProblem
               implicit none
               class(NSCompressibleProblem) :: a
            end subroutine
         end interface

         integer(ip) :: ndime

         !General operations for explicit solite
      
         !Update inner a%iteration counter and write headings in the solver file.
         a%itera = a%itera + 1
         if (a%MPIrank == a%MPIroot) then
            if(a%itera==1) write(a%lun_solve,100) a%istep
            write(a%lun_solve,101) a%itera
         endif
      
         call a%Timer%BuildLinearSystem%Tic
      
         call nsc_solite_ex(a)
      
         call a%Timer%BuildLinearSystem%Toc
         
         !No Algebraic system to solve
         call a%Timer%SolveLinearSystem%Tic
         
         call a%Timer%SolveLinearSystem%Toc
      
         !Update unknowns
         call a%Mesh%GetNdime(ndime)
   
         a%densf(:,1)   = a%unkno(1,:)
         a%momen(:,:,1) = a%unkno(2:ndime+1,:)
         a%energ(:,1)   = a%unkno(a%ndofn,:)
   
         call nsc_elmdir_ex(a)
   
         !Ghostcommunicate
         call nsc_GhostCommunicate(a)
   
         !HangingNodes for adaptivity

         !Formats. 
         100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
         101 format('------------------------------------------------------------', &
               /,'   INNER ITERATION NUMBER: ',i5)

      end subroutine

      subroutine nsc_SGSTimeUpdate(a)
         use typre
         use Mod_Mesh
         implicit none
         class(NSCompressibleExplicitProblem) :: a
      
         integer(ip)                :: ielem,nelem,ndime
      
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNelem(nelem)
      
         do ielem = 1,nelem
            a%cosgs(ielem)%a(1,:) = a%cosgs(ielem)%a(2,:)
            a%mosgs(ielem)%a(:,1,:)=a%mosgs(ielem)%a(:,2,:)
            a%ensgs(ielem)%a(1,:) = a%ensgs(ielem)%a(2,:)
            a%cosgs(ielem)%a(3,:) = a%cosgs(ielem)%a(2,:)
            a%mosgs(ielem)%a(:,3,:)=a%mosgs(ielem)%a(:,2,:)
            a%ensgs(ielem)%a(3,:) = a%ensgs(ielem)%a(2,:)
         end do
      
      end subroutine

      subroutine nsc_endste_ex(a,itask)
         use typre
         use def_parame
         implicit none

         class(NSCompressibleExplicitProblem) :: a      
         integer(ip) :: itask

         if (itask == 1) then
              
            !Residual and gradient projection 
            !Maybe redundant, but necessary in output cases
            if(a%kfl_repro == 1 .or. a%kfl_shock == 2)then
               call a%EndElmope('Endpro')
            end if
            
            call a%EndElmope('Endste')
      
            call a%CalculatePrimitivesLogic

         elseif (itask == 2) then

            if(a%kfl_timei==1) then
               a%densf(:,3) = a%densf(:,1)
               a%momen(:,:,3) = a%momen(:,:,1)
               a%energ(:,3) = a%energ(:,1)
            end if
            
         endif
         
      end subroutine

      subroutine nsc_exaerr_ex(a)
         use typre
         use def_parame
         implicit none

         class(NSCompressibleExplicitProblem) :: a      
   
         call a%CalculatePrimitives(1)

      end subroutine

      subroutine nsc_restar_ex(a,itask)
         use typre
         implicit none
         class(NSCompressibleExplicitProblem) :: a      
         integer(ip) :: itask

          
         if (itask == 1) then
               
            call a%CalculateConservatives(a%ncomp)

         elseif (itask == 2) then
            
            call a%CalculatePrimitives(a%ncomp)

         endif
   
      end subroutine

      subroutine nsc_NULLSUB(a)
         use typre
         implicit none
         class(NSCompressibleExplicitProblem) :: a
      
      end subroutine
   
      subroutine nsc_NULLSUBitask(a,itask)
         use typre
         implicit none
         class(NSCompressibleExplicitProblem) :: a
         integer(ip) :: itask
         
      end subroutine

      subroutine nsc_NULLSUBcharitask(a,itask)
         use typre
         implicit none
         class(NSCompressibleExplicitProblem) :: a
         character(6) :: itask
         
      end subroutine
   
      subroutine nsc_NULLSUBmarkel(a,markel)
         use typre
         implicit none
         class(NSCompressibleExplicitProblem) :: a
         integer(ip) :: markel(*)
         
      end subroutine

end module Mod_NSCompressibleExplicit

