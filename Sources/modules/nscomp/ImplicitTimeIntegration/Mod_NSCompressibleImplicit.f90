module Mod_NSCompressibleImplicit
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Mesh
   use Mod_NSCompressible
   implicit none
   private
   public NSCompressibleImplicitProblem, NSCompressibleImplicitProblem_Const
   
   type, extends(NSCompressibleProblem) :: NSCompressibleImplicitProblem
      
   real(rp), allocatable :: &
        dgrprj(:,:)                         !Density Gradient finite element projection
  
   real(rp) :: subrelax                           !Subrelaxation for nonlinearity

   real(rp) :: referencevelocity                   !Reference velocity for scaled norm

contains

      procedure :: Cvgunk            => nsc_cvgunk_im
      procedure :: Elmope            => nsc_elmope_im
      procedure :: Bouope            => nsc_bouope_im

      procedure :: EndElmope      => nsc_EndElmope_im 

      procedure :: SpecificNSCompReaphy    => nsc_NULLSUBitask
      procedure :: SpecificNSCompReanut    => nsc_reanut_im
      procedure :: SpecificNSCompReaMPI    => nsc_reampi_im
      procedure :: SpecificNSCompReabcs    => nsc_reabcs_im
      procedure :: SpecificNSCompMemall    => nsc_memall_im
      procedure :: SpecificNSCompTurnof    => nsc_turnof_im
      procedure :: SpecificNSCompBegste    => nsc_NULLSUB
      procedure :: SpecificNSCompBegite    => nsc_NULLSUB
      procedure :: SpecificNSCompEndite    => nsc_endite_im
      procedure :: SpecificNSCompEndste    => nsc_endste_im
      procedure :: SpecificNSCompExaerr    => nsc_exaerr_im
      procedure :: SpecificNSCompRefine    => nsc_refine_im
      procedure :: SpecificNSCompRefinementCriteria  => nsc_GetRefinementCriteria_im
      procedure :: SpecificSubscalesRefCriteria  => nsc_SubscalesRefCriteria
      procedure :: SpecificExactSolRefCriteria  => nsc_ExactSolRefCriteria
      procedure :: SpecificNSCompRestar    => nsc_restar_im

   end type

   interface
     
      subroutine nsc_cvgunk_im(a,itask)
         use typre
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem) :: a      
         integer(ip), intent(in) :: itask
      
      
      end subroutine

      subroutine nsc_elmope_im(a)
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem) :: a      
      
      end subroutine

      subroutine nsc_bouope_im(a)
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem) :: a      
      
      end subroutine

      subroutine nsc_EndElmope_im(NSCompImplicitProblem,task)
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem), target :: NSCompImplicitProblem
         character(6) :: task
      end subroutine

      subroutine nsc_endite_im(a,itask)
         use typre
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem) :: a      
         integer(ip) :: itask
      
      
      end subroutine

      subroutine nsc_SubscalesRefCriteria(a,error,TotalEstimatedError)
         use typre
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem) :: a
         real(rp) :: error(:), TotalEstimatedError
      end subroutine

      subroutine nsc_ExactSolRefCriteria(a,error)
         use typre
         import NSCompressibleImplicitProblem
         implicit none
         class(NSCompressibleImplicitProblem) :: a
         real(rp) :: error(:)
      end subroutine

   end interface
    
  interface NSCompressibleImplicitProblem_Const
      procedure constructor
  end interface NSCompressibleImplicitProblem_Const
   
contains

function constructor()
    class(NSCompressibleImplicitProblem), pointer :: constructor

    allocate(constructor)

end function constructor

      subroutine nsc_reanut_im(a,itask)
         use typre
         use MPI
         implicit none
   
         class(NSCompressibleImplicitProblem) :: a      
         integer(ip) :: itask
         
         !For Listener
         real(rp), pointer     :: param(:) => NULL()
         character(5), pointer :: words(:) => NULL()
         integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
         
         if (itask == 0) then
            a%subrelax = 1                                       ! Relaxation parameter
            a%referencevelocity = 1.0_rp                         ! Reference velocity
            
         elseif (itask == 1) then
            call a%Listener%getarrs(words,param,nnpar,nnwor)
            if(words(1)=='SUBRE') then
               if(words(2) == 'ON ') a%subrelax = param(2)
            elseif(words(1)=='REFVE') then
               if(words(2) == 'ON ') a%referencevelocity = param(2)
            end if

         elseif (itask == 100) then
            
         endif

      end subroutine

      subroutine nsc_reaMPI_im(a)
         use MPI
         implicit none
         class(NSCompressibleImplicitProblem) :: a      

         integer(ip) :: ierr
         
         CALL MPI_BCAST(a%subrelax, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%referencevelocity, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      end subroutine

      subroutine nsc_reabcs_im(a,itask)
         use typre
         use MPI
         use Mod_NscExacso
         implicit none
   
         class(NSCompressibleImplicitProblem) :: a      
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
               call a%Memor%alloc(ndime,exdeg,'exdeg','nsc_reabcs_im')   
               call a%Memor%alloc(ndime,exvel,'exvel','nsc_reabcs_im')   
               call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsc_reabcs_im')     
               call a%Memor%alloc(ndime,exteg,'exteg','nsc_reabcs_im')   
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
               call a%Memor%dealloc(ndime,exdeg,'exdeg','nsc_reabcs_im')   
               call a%Memor%dealloc(ndime,exvel,'exvel','nsc_reabcs_im')   
               call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsc_reabcs_im')     
               call a%Memor%dealloc(ndime,exteg,'exteg','nsc_reabcs_im')   
      
            endif
         endif

      end subroutine

      subroutine nsc_memall_im(a)
         use typre
         use Mod_Mesh
         use Mod_Memor
         use Mod_TimeIntegrator
         implicit none
         class(NSCompressibleImplicitProblem) :: a      

         type(TimeIntegratorDt1) :: Integrator
         integer(ip)             :: nsteps

         integer(ip) ::  ncomp,ndime,npoin

         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         
         !Set the Number of components necessary for the arrays
         call Integrator%Init(a%kfl_tsche_1st_datafile)
         call Integrator%GetNumberOfTimeSteps(nsteps)
         a%ncomp = 1 + nsteps
         ncomp = a%ncomp

         !Unknowns
         call a%Memor%alloc(npoin,ncomp,a%densf,'densf','nsc_memall_im')
         call a%Memor%alloc(ndime,npoin,ncomp,a%momen,'momen','nsc_memall_im')
         call a%Memor%alloc(npoin,ncomp,a%energ,'energ','nsc_memall_im')

         !Primitives
         call a%Memor%alloc(ndime,npoin,ncomp,a%veloc,'veloc','nsc_memall_im')
         call a%Memor%alloc(npoin,ncomp,a%press,'press','nsc_memall_im')
         call a%Memor%alloc(npoin,ncomp,a%tempe,'tempe','nsc_memall_im')

        !Density Gradient Orthogonal projection
        if(a%kfl_shock==2 .or. a%ErrorEstimatorTypeOfSubscales == 1) then
           call a%Memor%alloc(ndime,npoin,a%dgrprj,'dgrprj','nsc_memall_im')
        end if

      end subroutine

      subroutine nsc_turnof_im(a)
         use typre
         use Mod_Mesh
         use Mod_Memor
         use Mod_TimeIntegrator
         implicit none
         class(NSCompressibleImplicitProblem) :: a      

         integer(ip) :: ncomp,ndime,npoin

         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         
         ncomp = a%ncomp

         !Unknowns
         call a%Memor%dealloc(npoin,ncomp,a%densf,'densf','nsc_turnof_im')
         call a%Memor%dealloc(ndime,npoin,ncomp,a%momen,'momen','nsc_turnof_im')
         call a%Memor%dealloc(npoin,ncomp,a%energ,'energ','nsc_turnof_im')

         !Primitives
         call a%Memor%dealloc(ndime,npoin,ncomp,a%veloc,'veloc','nsc_turnof_im')
         call a%Memor%dealloc(npoin,ncomp,a%press,'press','nsc_turnof_im')
         call a%Memor%dealloc(npoin,ncomp,a%tempe,'tempe','nsc_turnof_im')

        !Density Gradient Orthogonal projection
        if(a%kfl_shock==2 .or. a%ErrorEstimatorTypeOfSubscales == 1) then
           call a%Memor%dealloc(ndime,npoin,a%dgrprj,'dgrprj','nsc_turnof_im')
        end if
      end subroutine

      subroutine nsc_endste_im(a,itask)
         use typre
         use def_parame
         implicit none

         class(NSCompressibleImplicitProblem) :: a      
         integer(ip) :: itask,icomp

         if (itask == 1) then
               
           call a%EndElmope('Endste')
      
           call a%CalculatePrimitivesLogic

         elseif (itask == 2) then
            
            !Higher order components
            if (a%ncomp > 3) then
               do icomp = a%ncomp, 4,-1
                  a%densf(:,icomp) = a%densf(:,icomp-1)
                  a%momen(:,:,icomp) = a%momen(:,:,icomp-1)
                  a%energ(:,icomp) = a%energ(:,icomp-1)
               enddo
            endif
      
            if(a%kfl_timei==1) then
               a%densf(:,3) = a%densf(:,1)
               a%momen(:,:,3) = a%momen(:,:,1)
               a%energ(:,3) = a%energ(:,1)
            end if

         endif
         
      end subroutine

      subroutine nsc_exaerr_im(a)
         use typre
         use def_parame
         implicit none

         class(NSCompressibleImplicitProblem) :: a      
   
         call a%CalculatePrimitives(1)

      end subroutine
      
      subroutine nsc_refine_im(a,itask)
         use typre
         use Mod_phpRefineArrays
         implicit none
         class(NSCompressibleImplicitProblem) :: a      
         character(6) :: itask

         integer(ip) :: oldnpoin,newnpoin,icomp,ndime
         real(rp), allocatable :: auxdensf(:,:), auxmomen(:,:,:), auxenerg(:,:)
         real(rp), allocatable :: auxpress(:,:), auxveloc(:,:,:), auxtempe(:,:)
         real(rp), allocatable :: auxdgrprj(:,:)
          
         !This subroutine transforms the tempe arrays from one mesh to the other 
         !Adaptive Mesh Refinement
         !Old Dimensions
         oldnpoin = size(a%densf,1)
      
         !We need to modify all the arrays
         !We assume that the mesh has already been updated
         call a%Mesh%GetNpoin(newnpoin)
         call a%Mesh%GetNdime(ndime)
         
         !Conservatives
         call a%Memor%alloc(newnpoin,a%ncomp,auxdensf,'densf','nsc_Refine_im')
         call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxmomen,'momen','nsc_Refine_im')
         call a%Memor%alloc(newnpoin,a%ncomp,auxenerg,'energ','nsc_Refine_im')
         do icomp = 1,a%ncomp
            if (itask == 'Refine') then
               call a%Refiner%UpdateVariable(1_ip,a%densf(:,icomp),auxdensf(:,icomp))
               call a%Refiner%UpdateVariable(ndime,a%momen(:,:,icomp),auxmomen(:,:,icomp))
               call a%Refiner%UpdateVariable(1_ip,a%energ(:,icomp),auxenerg(:,icomp))
            elseif (itask == 'Rebala') then
               call a%Refiner%RebalanceVariable(1_ip,a%densf(:,icomp),auxdensf(:,icomp))
               call a%Refiner%RebalanceVariable(ndime,a%momen(:,:,icomp),auxmomen(:,:,icomp))
               call a%Refiner%RebalanceVariable(1_ip,a%energ(:,icomp),auxenerg(:,icomp))
            endif   
         enddo
         call move_alloc(auxdensf,a%densf)
         call move_alloc(auxmomen,a%momen)
         call move_alloc(auxenerg,a%energ)
         call a%Memor%deallocObj(0,'densf','nsc_Refine_im',rp*oldnpoin*a%ncomp)
         call a%Memor%deallocObj(0,'momen','nsc_Refine_im',rp*ndime*oldnpoin*a%ncomp)
         call a%Memor%deallocObj(0,'energ','nsc_Refine_im',rp*oldnpoin*a%ncomp)
      
         !Primitives
         
         call a%Memor%alloc(newnpoin,a%ncomp,auxpress,'press','nsc_Refine_im')
         call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxveloc,'veloc','nsc_Refine_im')
         call a%Memor%alloc(newnpoin,a%ncomp,auxtempe,'tempe','nsc_Refine_im')
         do icomp = 1,a%ncomp
            if (itask == 'Refine') then
               call a%Refiner%UpdateVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
               call a%Refiner%UpdateVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
               call a%Refiner%UpdateVariable(1_ip,a%tempe(:,icomp),auxtempe(:,icomp))
            elseif (itask == 'Rebala') then
               call a%Refiner%RebalanceVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
               call a%Refiner%RebalanceVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
               call a%Refiner%RebalanceVariable(1_ip,a%tempe(:,icomp),auxtempe(:,icomp))
            endif   
         enddo
         call move_alloc(auxpress,a%press)
         call move_alloc(auxveloc,a%veloc)
         call move_alloc(auxtempe,a%tempe)
         call a%Memor%deallocObj(0,'press','nsc_Refine_im',rp*oldnpoin*a%ncomp)
         call a%Memor%deallocObj(0,'veloc','nsc_Refine_im',rp*ndime*oldnpoin*a%ncomp)
         call a%Memor%deallocObj(0,'tempe','nsc_Refine_im',rp*oldnpoin*a%ncomp)
         
         !Residual at the gauss points
         if ((a%kfl_shock == 1) .or. (a%npp_stepi(13) /= 0)) then
            call runend('Postprocess residual: Adaptive refinement not ready to transfer info on gauss points')
         endif
   
         !Gradient Orthogonal projection
         if(a%kfl_shock==2 .or. a%ErrorEstimatorTypeOfSubscales == 1) then
            call a%Memor%alloc(ndime,newnpoin,auxdgrprj,'dgrprj','nsc_Refine_im')
            if (itask == 'Refine') then
               call a%Refiner%UpdateVariable(ndime,a%dgrprj(:,:),auxdgrprj(:,:))
            elseif (itask == 'Rebala') then
               call a%Refiner%RebalanceVariable(ndime,a%dgrprj(:,:),auxdgrprj(:,:))
            endif   
            call move_alloc(auxdgrprj,a%dgrprj)
            call a%Memor%deallocObj(0,'dgrprj','nsc_Refine_im',rp*ndime*oldnpoin)
         end if

      end subroutine

      subroutine nsc_GetRefinementCriteria_im(a,markel)
         use typre
         use Mod_ZZErrorEstimator
      
         implicit none
         class(NSCompressibleImplicitProblem) :: a
         integer(ip) :: markel(*)
         integer(ip) :: ndime,nelem
         real(rp)    :: TotalEstimatedError
         real(rp), allocatable :: error(:)
         
         call a%Mesh%GetNelem(nelem)
         call a%Memor%alloc(nelem,error,'error','nsc_GetRefinementCriteria_im')
         call a%Mesh%GetNdime(ndime)
         if (a%RefinerErrorEstimator == 'ZZ   ') then 
            call ZZErrorEstimator(a%Mesh,a%Memor,ndime,a%momen(:,:,1),error)
         elseif (a%RefinerErrorEstimator == 'GRADI') then
            call GradientErrorEstimator(a%Mesh,a%Memor,ndime,a%momen(:,:,1),error)
         elseif (a%RefinerErrorEstimator == 'SUBSC') then
            call a%SpecificSubscalesRefCriteria(error,TotalEstimatedError)
         endif   
      
         call a%FilePostpr%postgp(error,'Error',a%istep,a%ctime,a%Mesh)
         
         call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)
         
         call a%Memor%dealloc(nelem,error,'error','nsc_GetRefinementCriteria_im')
            
      end subroutine

      subroutine nsc_restar_im(a,itask)
         use typre
         implicit none
         class(NSCompressibleImplicitProblem) :: a      
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
         class(NSCompressibleImplicitProblem) :: a
      
      end subroutine

      subroutine nsc_NULLSUBitask(a,itask)
         use typre
         implicit none
         class(NSCompressibleImplicitProblem) :: a
         integer(ip) :: itask
         
      end subroutine

end module Mod_NSCompressibleImplicit

