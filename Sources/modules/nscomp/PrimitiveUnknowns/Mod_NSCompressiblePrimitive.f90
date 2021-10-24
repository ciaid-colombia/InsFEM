module Mod_NSCompressiblePrimitive
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Mesh
   use Mod_NSCompressible
   implicit none
   private
   public NSCompressiblePrimitiveProblem, NSCompressiblePrimitiveProblem_Const

   type, extends(NSCompressibleProblem) :: NSCompressiblePrimitiveProblem

      !Compressible/Incompressible limit for velocity in tau
      real(rp) :: compinc

      !Subrelaxation for nonlinearity
      real(rp) :: subrelax
      
      !Conservative Fluxes 
      real(rp), allocatable :: ConservativeFluxes(:,:)

   contains

      procedure :: Getste            => nsc_pr_getste
      procedure :: Cvgunk            => nsc_pr_cvgunk
      procedure :: Elmope            => nsc_pr_elmope
      procedure :: Bouope            => nsc_pr_bouope
      procedure :: SpecificIniunk    => nsc_pr_iniunk
      procedure :: InnerResiduals    => nsc_pr_InnerResiduals
      procedure :: EndElmope         => nsc_pr_EndElmope
      procedure :: EndBouope         => nsc_pr_EndBouope

      procedure :: SpecificNSCompReaphy    => nsc_pr_reaphy
      procedure :: SpecificNSCompReanut    => nsc_pr_reanut
      procedure :: SpecificNSCompReaMPI    => nsc_pr_reampi
      procedure :: SpecificNSCompReabcs    => nsc_pr_reabcs
      procedure :: SpecificNSCompMemall    => nsc_pr_memall
      procedure :: SpecificNSCompTurnof    => nsc_pr_turnof
      procedure :: SpecificNSCompBegste    => nsc_pr_begste
      procedure :: SpecificNSCompBegite    => nsc_pr_begite
      procedure :: SpecificNSCompEndite    => nsc_pr_endite
      procedure :: SpecificNSCompEndste    => nsc_pr_endste
      procedure :: SpecificNSCompExaerr    => nsc_NULLSUB
      procedure :: SpecificNSCompRefine    => nsc_pr_refine
      procedure :: SpecificNSCompRefinementCriteria  => nsc_pr_GetRefinementCriteria
      procedure :: SpecificSubscalesRefCriteria  => nsc_pr_SubscalesRefCriteria
      procedure :: SpecificExactSolRefCriteria  => nsc_pr_ExactSolRefCriteria
      procedure :: SpecificNSCompRestar    => nsc_pr_restar
      
      ! ROM
      procedure :: SpecificGetUnkno     => nsc_GetUnkno

      procedure :: ComputeNscompNstincVelocity

      ! Conservation Restriction
      procedure :: ComputeConservativeFluxes => nsc_pr_ComputeConservativeFluxes

   end type

   interface

      subroutine nsc_pr_getste(a,dtinv)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      
         real(rp) :: dtinv


      end subroutine

      subroutine nsc_pr_cvgunk(a,itask)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      
         integer(ip), intent(in) :: itask


      end subroutine

      subroutine nsc_pr_iniunk(a)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      


      end subroutine

      subroutine nsc_pr_elmope(a)
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      

      end subroutine

      subroutine nsc_pr_bouope(a)
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      

      end subroutine

      subroutine nsc_pr_EndElmope(NSCompPrimitiveProblem,task)
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem), target :: NSCompPrimitiveProblem
         character(6) :: task
      end subroutine

      subroutine nsc_pr_endite(a,itask)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      
         integer(ip) :: itask


      end subroutine

      subroutine nsc_pr_endBouope(a)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      


      end subroutine

      subroutine nsc_pr_ComputeConservativeFluxes(a,task)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a      
         character(6) :: task


      end subroutine

      subroutine nsc_pr_SubscalesRefCriteria(a,error,TotalEstimatedError)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a
         real(rp) :: error(:), TotalEstimatedError

      end subroutine

      subroutine nsc_pr_ExactSolRefCriteria(a,error)
         use typre
         import NSCompressiblePrimitiveProblem
         implicit none
         class(NSCompressiblePrimitiveProblem) :: a
         real(rp) :: error(:)

      end subroutine

   end interface

   interface NSCompressiblePrimitiveProblem_Const
      procedure constructor
   end interface NSCompressiblePrimitiveProblem_Const

contains

   function constructor()
      class(NSCompressiblePrimitiveProblem), pointer :: constructor

      allocate(constructor)

   end function constructor

   subroutine nsc_pr_reaphy(a,itask)
      use typre
      implicit none

      class(NSCompressiblePrimitiveProblem) :: a      
      integer(ip) :: itask

      !For Listener
      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()

      if (itask == 0) then

      elseif (itask == 1) then

         call a%Listener%getarrs(words,param,nnpar,nnwor)
         if(words(1) == 'RESTR')then
            if(words(2)=='ON') then
               a%kfl_restrictions = 1_ip
            endif
         endif 

      elseif (itask == 2) then

      elseif (itask == 100) then

      endif

   end subroutine

   subroutine nsc_pr_reanut(a,itask)
      use typre
      implicit none

      class(NSCompressiblePrimitiveProblem) :: a      
      integer(ip) :: itask

      !For Listener
      real(rp), pointer     :: param(:) => NULL()
      character(5), pointer :: words(:) => NULL()
      integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()

      if (itask == 0) then
         a%compinc = 0.3       ! Physical compressible/incompressible limit
         a%subrelax = 1        ! Relaxation parameter

      elseif (itask == 1) then
         call a%Listener%getarrs(words,param,nnpar,nnwor)
         if(words(1)=='COMPI') then
            a%compinc = param(1)
         else if(words(1)=='SUBRE') then
            if(words(2) == 'ON ') a%subrelax = param(2)
         end if

      elseif (itask == 100) then

      endif

   end subroutine

   subroutine nsc_pr_reaMPI(a)
      use MPI
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a      

      integer(ip) :: ierr

      CALL MPI_BCAST(a%compinc, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(a%subrelax, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)

   end subroutine

   subroutine nsc_pr_reabcs(a,itask)
      use typre
      use Mod_NscExacso
      implicit none

      class(NSCompressiblePrimitiveProblem) :: a      
      integer(ip) :: itask

      type(NscExacso) :: exacso     

      !Exact Values
      real(rp) :: expre, extem
      real(rp), allocatable   :: exprg(:),exvel(:),exveg(:,:),exteg(:)
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
            call a%Memor%alloc(ndime,exprg,'exprg','nsc_pr_reabcs')   
            call a%Memor%alloc(ndime,exvel,'exvel','nsc_pr_reabcs')   
            call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsc_pr_reabcs')     
            call a%Memor%alloc(ndime,exteg,'exteg','nsc_pr_reabcs')   
            a%kfl_fixbo = -1
            a%kfl_fixno = -1
            a%bvess = 0

            do ipoin = 1,npoin
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
               if (ibopo /= 0) then
                  call a%Mesh%GetPointCoord(ipoin,coord)
                  call exacso%nsc_ComputeSolution(ndime,coord,a)
                  call exacso%nsc_GetPressure(ndime,expre,exprg)           
                  call exacso%nsc_GetVelocity(ndime,exvel,exveg)           
                  call exacso%nsc_GetTemperature(ndime,extem,exteg)           

                  !a%kfl_fixno(1,ipoin) = 1_ip
                  a%kfl_fixno(2,ipoin) = 1_ip
                  a%kfl_fixno(3,ipoin) = 1_ip
                  if (ndime ==3) a%kfl_fixno(4,ipoin) = 1_ip
                  a%kfl_fixno(a%ndofn,ipoin) = 1_ip
                  a%bvess(1,ipoin,1) = expre
                  a%bvess(2:ndime+1,ipoin,1) = exvel(:)
                  a%bvess(4,ipoin,1) = extem
               endif
            enddo   
            ! Deallocate exact components
            call a%Memor%dealloc(ndime,exprg,'exprg','nsc_pr_reabcs')   
            call a%Memor%dealloc(ndime,exvel,'exvel','nsc_pr_reabcs')   
            call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsc_pr_reabcs')     
            call a%Memor%dealloc(ndime,exteg,'exteg','nsc_pr_reabcs')   

         endif
      endif

   end subroutine

   subroutine nsc_pr_memall(a)
      use typre
      use Mod_Mesh
      use Mod_Memor
      use Mod_TimeIntegrator
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a      

      type(TimeIntegratorDt1) :: Integrator
      integer(ip)             :: nsteps

      integer(ip) :: ncomp,ndime,npoin

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)

      !Set the Number of components necessary for the arrays
      call Integrator%Init(a%kfl_tsche_1st_datafile)
      call Integrator%GetNumberOfTimeSteps(nsteps)
      a%ncomp = 1 + nsteps
      ncomp = a%ncomp

      !Unknowns
      call a%Memor%alloc(npoin,ncomp,a%press,'press','nsc_pr_memall')
      call a%Memor%alloc(ndime,npoin,ncomp,a%veloc,'veloc','nsc_pr_memall')
      call a%Memor%alloc(npoin,ncomp,a%tempe,'tempe','nsc_pr_memall')

      !Conservatives
      call a%Memor%alloc(npoin,ncomp,a%densf,'densf','nsc_pr_memall')
      call a%Memor%alloc(ndime,npoin,ncomp,a%momen,'momen','nsc_pr_memall')
      call a%Memor%alloc(npoin,ncomp,a%energ,'energ','nsc_pr_memall')

      !Interpolation with Restrictions   
      if (a%kfl_restrictions==1) then
         call a%Memor%alloc(5,2,a%ConservativeFluxes,'ConservativeFluxes','nsc_pr_memall')
      end if

   end subroutine

   subroutine nsc_pr_turnof(a)
      use typre
      use Mod_Mesh
      use Mod_Memor
      use Mod_TimeIntegrator
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a      

      integer(ip) :: ncomp,ndime,npoin

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)

      ncomp = a%ncomp

      !Unknowns
      call a%Memor%dealloc(npoin,ncomp,a%press,'press','nsc_pr_turnof')
      call a%Memor%dealloc(ndime,npoin,ncomp,a%veloc,'veloc','nsc_pr_turnof')
      call a%Memor%dealloc(npoin,ncomp,a%tempe,'tempe','nsc_pr_turnof')

      !Conservatives
      call a%Memor%dealloc(npoin,ncomp,a%densf,'densf','nsc_pr_turnof')
      call a%Memor%dealloc(ndime,npoin,ncomp,a%momen,'momen','nsc_pr_turnof')
      call a%Memor%dealloc(npoin,ncomp,a%energ,'energ','nsc_pr_turnof')

      if (a%kfl_restrictions==1) then
         call a%Memor%dealloc(5,2,a%ConservativeFluxes,'ConservativeFluxes','nsc_pr_memall')
      end if
   end subroutine

   subroutine nsc_pr_begste(a)
      use typre
      use def_parame
      implicit none

      class(NSCompressiblePrimitiveProblem) :: a      

      if(a%kfl_timei==1) then
         a%press(:,2) =   a%press(:,3)
         a%veloc(:,:,2) = a%veloc(:,:,3)
         a%tempe(:,2) =   a%tempe(:,3)
         a%press(:,1) =   a%press(:,2)
         a%veloc(:,:,1) = a%veloc(:,:,2)
         a%tempe(:,1) =   a%tempe(:,2)
      end if

      !Conservation Restrictions
      if(a%kfl_restrictions == 1)then
         call a%ComputeConservativeFluxes('BegSte')
      end if
   end subroutine

   subroutine nsc_pr_begite(a)
      use typre
      use def_parame
      implicit none

      class(NSCompressiblePrimitiveProblem) :: a      

      !Assign var(n,i,0) <-- var(n,i-1,*), initial guess for inner iterations

      if(a%kfl_timei==1) then
         a%press(:,1) =   a%press(:,2)
         a%veloc(:,:,1) = a%veloc(:,:,2)
         a%tempe(:,1) =   a%tempe(:,2)
      end if

      !Conservation Restrictions
      if(a%kfl_restrictions == 1)then
         call a%ComputeConservativeFluxes('BegIte')
      end if
   end subroutine

   subroutine nsc_pr_endste(a,itask)
      use typre
      use def_parame
      implicit none

      class(NSCompressiblePrimitiveProblem) :: a      
      integer(ip) :: itask,icomp,nelem,ielem

      if (itask == 1) then

         call a%EndElmope('Endste')

      elseif (itask == 2) then

         !Higher order components
         if (a%ncomp > 3) then
            do icomp = a%ncomp, 4,-1
               a%press(:,icomp) = a%press(:,icomp-1)
               a%veloc(:,:,icomp) = a%veloc(:,:,icomp-1)
               a%tempe(:,icomp) = a%tempe(:,icomp-1)
            enddo
         endif

         if(a%kfl_timei==1) then
            a%press(:,3) = a%press(:,1)
            a%veloc(:,:,3) = a%veloc(:,:,1)
            a%tempe(:,3) = a%tempe(:,1)
         end if

         ! Update the subgrid scales
         if(a%kfl_tacsg/=0) then
            call a%Mesh%GetNelem(nelem)
            if (a%kfl_tacsg>0) then
               !Crank Nicolson schemes
               if (a%kfl_tsche_1st_current == 'CN   ') then   
                  do ielem=1,nelem
                     a%cosgs(ielem)%a(2,:) = 2.0_rp*a%cosgs(ielem)%a(1,:) - a%cosgs(ielem)%a(2,:)
                     a%mosgs(ielem)%a(:,2,:)=2.0_rp*a%mosgs(ielem)%a(:,1,:)-a%mosgs(ielem)%a(:,2,:)
                     a%ensgs(ielem)%a(2,:) = 2.0_rp*a%ensgs(ielem)%a(1,:) -a%ensgs(ielem)%a(2,:)
                  end do
               else if(a%kfl_tacsg>0) then
                  do ielem=1,nelem
                     a%cosgs(ielem)%a(2,:) = a%cosgs(ielem)%a(1,:) 
                     a%mosgs(ielem)%a(:,2,:)=a%mosgs(ielem)%a(:,1,:)
                     a%ensgs(ielem)%a(2,:) = a%ensgs(ielem)%a(1,:) 
                  end do
               end if
            endif
         end if

         !Calculate conservative variables for postprocess
         if(a%kfl_restrictions /= 1) then
            if((a%npp_stepi(4)==1).or.(a%npp_stepi(5)==1).or.(a%npp_stepi(6)==1)) then
               call a%CalculateConservatives(1)
            endif
         endif

      endif

   end subroutine

   subroutine nsc_pr_InnerResiduals(a,drnsc,mrnsc,ernsc)
      use typre
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a
      real(rp) :: drnsc,mrnsc,ernsc
      real(rp) :: zensi = epsilon(0.0_rp)        ! zero

      integer(ip) :: ndime, npoinLocal

      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPIHeterogeneous(2_ip,a%ndofn,1,npoinLocal,1,1,1,a%unkno,a%press(1,1),drnsc,zensi,a%MPIcomm,a%MPIroot)     
      call vecresMPIHeterogeneous(2_ip,a%ndofn,ndime,npoinLocal,2,1,ndime,a%unkno,a%veloc(1,1,1),mrnsc,zensi,a%MPIcomm,a%MPIroot)
      call vecresMPIHeterogeneous(2_ip,a%ndofn,1,npoinLocal,a%ndofn,1,1,a%unkno,a%tempe(1,1),ernsc,zensi,a%MPIcomm,a%MPIroot)     

   end subroutine

   subroutine nsc_pr_refine(a,itask)
      use typre
      use Mod_phpRefineArrays
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a      
      character(6) :: itask

      integer(ip) :: oldnpoin,newnpoin,icomp,ndime
      real(rp), allocatable :: auxdensf(:,:), auxmomen(:,:,:), auxenerg(:,:)
      real(rp), allocatable :: auxpress(:,:), auxveloc(:,:,:), auxtempe(:,:)

      !This subroutine transforms the tempe arrays from one mesh to the other 
      !Adaptive Mesh Refinement
      !Old Dimensions
      oldnpoin = size(a%press,1)

      !We need to modify all the arrays
      !We assume that the mesh has already been updated
      call a%Mesh%GetNpoin(newnpoin)
      call a%Mesh%GetNdime(ndime)

      !Conservatives
      call a%Memor%alloc(newnpoin,a%ncomp,auxdensf,'densf','nsc_pr_Refine')
      call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxmomen,'momen','nsc_pr_Refine')
      call a%Memor%alloc(newnpoin,a%ncomp,auxenerg,'energ','nsc_pr_Refine')
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
      call a%Memor%deallocObj(0,'densf','nsc_pr_Refine',rp*oldnpoin*a%ncomp)
      call a%Memor%deallocObj(0,'momen','nsc_pr_Refine',rp*ndime*oldnpoin*a%ncomp)
      call a%Memor%deallocObj(0,'energ','nsc_pr_Refine',rp*oldnpoin*a%ncomp)

      !Primitives

      call a%Memor%alloc(newnpoin,a%ncomp,auxpress,'press','nsc_pr_Refine')
      call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxveloc,'veloc','nsc_pr_Refine')
      call a%Memor%alloc(newnpoin,a%ncomp,auxtempe,'tempe','nsc_pr_Refine')
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
      call a%Memor%deallocObj(0,'press','nsc_pr_Refine',rp*oldnpoin*a%ncomp)
      call a%Memor%deallocObj(0,'veloc','nsc_pr_Refine',rp*ndime*oldnpoin*a%ncomp)
      call a%Memor%deallocObj(0,'tempe','nsc_pr_Refine',rp*oldnpoin*a%ncomp)

      !Residual at the gauss points
      if ((a%kfl_shock == 1) .or. (a%npp_stepi(13) /= 0)) then
         call runend('Postprocess residual: Adaptive refinement not ready to transfer info on gauss points')
      endif

   end subroutine

   subroutine nsc_pr_GetRefinementCriteria(a,markel)
      use typre
      use Mod_ZZErrorEstimator

      implicit none
      class(NSCompressiblePrimitiveProblem) :: a
      integer(ip) :: markel(*)
      integer(ip) :: ndime,nelem
      real(rp)    :: TotalEstimatedError
      real(rp), allocatable :: error(:)

      call a%Mesh%GetNelem(nelem)
      call a%Memor%alloc(nelem,error,'error','nsc_pr_GetRefinementCriteria')
      call a%Mesh%GetNdime(ndime)
      if (a%RefinerErrorEstimator == 'ZZ   ') then 
         call ZZErrorEstimator(a%Mesh,a%Memor,ndime,a%veloc(:,:,1),error)
      elseif (a%RefinerErrorEstimator == 'GRADI') then
         call GradientErrorEstimator(a%Mesh,a%Memor,ndime,a%veloc(:,:,1),error)
      elseif (a%RefinerErrorEstimator == 'SUBSC') then
         call a%SpecificSubscalesRefCriteria(error,TotalEstimatedError)
      endif   

      call a%FilePostpr%postgp(error,'Error',a%istep,a%ctime,a%Mesh)

      call ApplyErrorCriteria(a%MPIcomm,nelem,a%RefinerErrorCriteria,a%RefinerErrorLimits,error,markel)

      call a%Memor%dealloc(nelem,error,'error','nsc_pr_GetRefinementCriteria')

   end subroutine

   subroutine nsc_pr_restar(a,itask)
      use typre
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a      
      integer(ip) :: itask


      if (itask == 1) then


      elseif (itask == 2) then


      endif

   end subroutine

   !Compute characteristic velocity for incompressible limit
   subroutine ComputeNscompNstincVelocity(a,vel,cspd,gpvst)
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a      
      real(rp),    intent(in)    :: vel,cspd
      real(rp),    intent(out)   :: gpvst
      real(rp)                   :: mach,xerf

      !Compressible/Incompressible characteristic velocity                                 
      mach = vel/cspd
      xerf = 2.0_rp-2.0_rp*(a%compinc-mach)/a%compinc
      gpvst = vel + erf(xerf)*cspd

   end subroutine

   subroutine nsc_NULLSUB(a)
      use typre
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a

   end subroutine

   subroutine nsc_NULLSUBitask(a,itask)
      use typre
      implicit none
      class(NSCompressiblePrimitiveProblem) :: a
      integer(ip) :: itask
   end subroutine

   subroutine nsc_GetUnkno(a,unkno)
      use typre
      implicit none
      class(NSCompressiblePrimitiveProblem)  :: a
      real(rp) :: unkno(:,:)
            
   end subroutine

end module Mod_NSCompressiblePrimitive

