subroutine php_reabcs(a)
!------------------------------------------------------------------------
!    This routine reads the boundary conditions for a general physical problem.
!
!    For conditions on nodes, gbvess(idime,ipoin,1) contains the
!    unkno value.
!    The different codes for kfl_fixno(ipoin) are:
!    = 1 ... Dirichlet
!    = 0 ... Free or initial
!
!    For conditions on boundaries, a%gbvnat(iboun) is allocated ONLY if
!    a condition is imposed on it. Overmore, its length depends on the
!    condition type. The different codes for kfl_fixbo(iboun) are:
!    = 1 ... Dirichlet ............ u
!    = 2 ... Neumann condition .... sig.n=-p n

!    At the end of the subroutine Dirichlet conditions on boundaries are
!    transfered to conditions on nodes.
!------------------------------------------------------------------------
   use MPI
   use typre
   use Mod_PhysicalProblem
   use Mod_ToCSR
   use Mod_OneToAllBuffer
   use Mod_OneToAllLoop
   use Mod_OneToAllReadOnBoundaries
   use Mod_OneToAllReadOnElements
   use Mod_OneToAllReadOnNodes   
   use Mod_OneToAllReadSource
   implicit none
   class(PhysicalProblem) :: a
   
   character(150) :: outstr
   integer(ip)    :: npoin,npoinLocal,nboun,ndime
   
   integer(ip), allocatable :: knodb(:),aux_knodb(:)
   
   integer(ip), allocatable :: iBounFlag(:), iProcList(:)
   
   integer(ip) :: kfl_Reading, kfl_doreceiver
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer            :: ierr
   integer            :: irequest01(a%MPIsize), irequest02(a%MPIsize), irequest03(a%MPIsize)
   integer            :: status(MPI_STATUS_SIZE)
   
   
   integer(ip) :: idummy
   integer(ip) :: irank,sendFlag,SendFlag_Buffer(a%MPIsize)
   
   integer(ip) :: ipoin,inodb,mnodb,pnodb,ifunc,ifunp,nfunp
   
   integer(ip) :: kfl_perio
   
   type(OneToAllReadOnBoundaries) :: OTAReadOnBoundaries
   type(OneToAllReadOnElements) :: OTAReadOnElements
   type(OneToAllReadOnNodes)      :: OTAReadOnNodes   
   type(OneToAllReadSource)      :: OTAReadSource   
   
   
   interface
      subroutine php_InitializeConstantBoundaryConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_InitializeInitialConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_InitializeOnBoundariesConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_InitializeOnElementsConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_ReadConstantBoundaryConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_ReadInitialConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_ReadOnBoundaries(a,knodb)
         use typre
         import
         implicit none
         class(PhysicalProblem) :: a
         integer(ip)            :: knodb(*)
      
      end subroutine
      
      subroutine php_ReadOnElements(a,lnode)
         use typre
         import
         implicit none
         class(PhysicalProblem) :: a
         integer(ip)            :: lnode(*)
      
      end subroutine
      
      subroutine php_ReadOnFunctions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
   
      subroutine php_ReadOnNodes(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_ScatterConstantBoundaryConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_ScatterInitialConditions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_ScatterOnFunctions(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_FinalizeOnBoundaries(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
      
      subroutine php_FinalizeOnElements(a)
         import
         implicit none
         class(PhysicalProblem) :: a
      
      end subroutine
   
   end interface
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   
   !Output
   outstr = adjustl(trim(a%exmod))//'_reabcs'
   
   !Dimensions
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetMnodb(mnodb)
   
   !-------------------------------------------------------------------------------
   !Initializations
   call php_InitializeConstantBoundaryConditions(a)
   call php_InitializeInitialConditions(a)
   call php_InitializeOnBoundariesConditions(a)
   call php_InitializeOnElementsConditions(a)
   
   !Specific Initializations
   call a%SpecificReabcs(0)
   
   !Initialize listener in root
   if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
      call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
      
      !Reach the boundary conditions section.
      call a%Listener%listen(outstr)
      do while(a%Listener%words(1)/='BOUND')
         call a%Listener%listen(outstr)
      end do
      
      !-------------------------------------------------------------------------------
      !CONSTANT/VARIABLE BOUNDARY CONDITIONS
      call php_ReadConstantBoundaryConditions(a)
   endif
   
   !-------------------------------------------------------------------------------
   !CONSTANT/VARIABLE BOUNDARY CONDITIONS
   if (a%kfl_ReadType == 0) call php_ScatterConstantBoundaryConditions(a)
   
   !-------------------------------------------------------------------------------
   !MEMORY ALLOCATION
   call a%Memor%alloc(a%ndofbc,npoin,a%kfl_fixno,'kfl_fixno',outstr)
   a%kfl_fixno = -1
   call a%Memor%alloc(nboun,a%kfl_fixbo,'kfl_fixbo',outstr)
   a%kfl_fixbo = -1
   call a%Memor%alloc(nboun,a%bvnat,'bvnat',outstr)
   if (a%kfl_conbc == 1) then
      call a%Memor%alloc(a%ndofbc,npoin,1,a%bvess,'bvess',outstr)
   else
      call a%Memor%alloc(npoin,a%kfl_funno,'kfl_funno',outstr)
      call a%Memor%alloc(nboun,a%kfl_funbo,'kfl_funbo',outstr)
      call a%Memor%alloc(10,2 ,a%kfl_funty,'kfl_funty',outstr)
      call a%Memor%alloc(10   ,a%funpa,'funpa',outstr)
      call a%Memor%alloc(a%ndofbc,npoin,2,a%bvess,'bvess',outstr)
   end if
   
   !Specific Header
   call a%SpecificReabcs(1)
   
   !Default is read and scatter boundary conditions
   if (a%kfl_BoundaryConditionsReadStrategy == 0) then
      !Initialize irequests
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         irequest01 = MPI_REQUEST_NULL
         irequest02 = MPI_REQUEST_NULL
         irequest03 = MPI_REQUEST_NULL
      endif

      !Start Reading
      kfl_Reading = 0
      do while (kfl_Reading /= -1)
         !Default for all is receive
         kfl_doreceiver = 1
         
         !First level
         if (kfl_Reading == 0) then
               
            !To do by root
            if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         
               !Listen
               call a%Listener%listen(outstr)
               
               sendFlag = 0
               if (a%Listener%words(1) == 'ENDBO') then
                  !Tell everybody we are done
                  sendFlag = -1
               elseif(a%Listener%words(1) == 'ONNOD') then
                  !Tell everybody we are on nodes
                  sendFlag = 1
               elseif(a%Listener%words(1) == 'INITI') then
               
                  !-------------------------------------------------------------------------------
                  !INITIAL CONDITIONS
                  
                  !Initial conditions have the data on the same line
                  call php_ReadInitialConditions(a)
                  !Continue reading
                  sendFlag = 0    
                  
               elseif(a%Listener%words(1) == 'ONBOU') then
                  !Tell everybody we are on boundary conditions
                  sendFlag = 3  
               elseif(a%Listener%words(1) == 'FUNCT') then
                  
                  !-------------------------------------------------------------------------------
                  !FUNCTIONS CONDITIONS
                  
                  !Read the function information, but do not send anything (yet)
                  call php_ReadOnFunctions(a)
                  !Continue reading
                  sendFlag = 6    
               elseif(a%Listener%words(1) == 'SOURC') then
                  !Tell everybody we are on boundary conditions
                  sendFlag = 4  
               elseif(a%Listener%words(1) == 'ONELE') then
                  !Tell everybody we are on boundary conditions
                  sendFlag = 5  
               else
                  !Module Specific Boundary conditions
                  call a%SpecificReabcs(2)
               endif
               kfl_Reading = sendFlag
               if (a%kfl_ReadType == 0) then
                  do irank = 0,a%MPIsize-1
                     !make sure all is sent before sending again
                     call MPI_WAIT(irequest01(irank+1), status, ierr)
                     sendFlag_buffer(irank+1) = sendFlag
                     call MPI_ISEND(sendFlag_buffer(irank+1),1, MPI_INTEGER4, irank, mtag1, a%MPIcomm,irequest01(irank+1), ierr) 
                  enddo
               endif
            endif
            
            !To do by everybody
            if (kfl_doreceiver == 1 .and. a%kfl_ReadType == 0) then
               !Receive information on what we are doing
               call MPI_RECV(kfl_Reading,1_ip, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,status, ierr) 
            endif
            
            
         !-------------------------------------------------------------------------------
         !ON NODES BOUNDARY CONDITIONS
         elseif (kfl_Reading == 1) then
            call a%OTAReadOnNodesDo
            
            kfl_Reading = 0  
         !-------------------------------------------------------------------------------
         !ON BOUNDARIES BOUNDARY CONDITIONS
         elseif (kfl_Reading == 3) then
            call a%OTAReadOnBoundariesDo
            
            kfl_Reading = 0      
         !-------------------------------------------------------------------------------
         !SOURCE TERM
         elseif (kfl_Reading == 4) then
            call a%OTAReadSourceDo
            
            kfl_Reading = 0  
         !ON ELEMENTS BOUNDARY CONDITIONS
         elseif (kfl_Reading == 5) then
            call a%OTAReadOnElementsDo
            kfl_Reading = 0      
         !SPECIAL BOUNDARY CONDITIONS
         !-------------------------------------------------------------------------------
         !ON FUNCTIONS BOUNDARY CONDITIONS
         elseif (kfl_Reading == 6) then
            if (a%kfl_ReadType == 0) then
               !On Functions Conditions
               call php_ScatterOnFunctions(a)
            endif
            
            kfl_Reading = 0      
         else
            !To do by root
            if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
               !Listen
               call a%Listener%listen(outstr)
               !Module Specific Boundary conditions
               call a%SpecificReabcs(3,kfl_Reading)
            endif
            
            !To do by everybody
            if (kfl_doreceiver == 1 .and. a%kfl_ReadType == 0) then
               call a%SpecificReabcs(4,kfl_Reading)
            endif
         endif
      enddo
      !Wait if necessary 
      !make sure all is sent and received
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         do irank = 0,a%MPIsize-1
            call MPI_WAIT(irequest01(irank+1), status, ierr)
            call MPI_WAIT(irequest02(irank+1), status, ierr)
            call MPI_WAIT(irequest03(irank+1), status, ierr)
         enddo
      endif
   
   !Manufactured boundary conditions
   elseif (a%kfl_BoundaryConditionsReadStrategy == 1) then
      call a%SpecificManufacturedBoundaryConditions
      
   endif
   
   !----------------------------------------------
   !SCATTERING
   if (a%kfl_ReadType == 0) then
      !Initial Conditions
      call php_ScatterInitialConditions(a)
   endif
   

   !----------------------------------------------------------------------------
   !ON BOUNDARIES
   call php_FinalizeOnBoundaries(a)
   
   !----------------------------------------------------------------------------
   !ON BOUNDARIES
   call php_FinalizeOnElements(a)
   
   !specific reabcs before Bounod

   call a%SpecificReabcs(5)

   call a%Bounod

   
   !FINALIZATIONS
   !-------------------------------------------------------------------------
   !COMMUNICATE NODAL GHOST DATA WITH NEIGHBOURS
   call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofbc,a%kfl_fixno)
   call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofbc,a%bvess(:,:,1))
   if (a%kfl_conbc /=1) then
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,a%kfl_funno)
      call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofbc,a%bvess(:,:,2))
   endif
   
   !This is not necessary, but convenient for tests, which had different values for kfl_fixno depending on type of communicator
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%MasterToSlave(a%ndofbc,a%kfl_fixno)
   
   !-------------------------------------------------------------------------
   !SPECIFIC FINALIZE
   call a%SpecificReabcs(100)
   
   !Post Process boundary conditions

   if (a%kfl_printBC) then
      if (a%ndofbc == 1) then
          call a%FilePostpr%postpr(a%bvess(1,:,1),adjustl(trim(a%exmod))//'_BC',0_ip,0.0_rp,a%Mesh)
          call a%FilePostpr%postpr(a%kfl_fixno(1,:),adjustl(trim(a%exmod))//'_FIXNO',0_ip,0.0_rp,a%Mesh)
      elseif (a%ndofbc == ndime) then
          call a%FilePostpr%postpr(a%bvess(:,:,1),adjustl(trim(a%exmod))//'_BC',0_ip,0.0_rp,a%Mesh)
          call a%FilePostpr%postpr(a%kfl_fixno,adjustl(trim(a%exmod))//'_FIXNO',0_ip,0.0_rp,a%Mesh)
      endif
   endif
   
end subroutine

subroutine php_OTAReadOnNodesDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllReadOnNodes
   class(PhysicalProblem) :: a
   type(OneToAllReadOnNodes) ::   OTAReadOnNodes   
   
   call OTAReadOnNodes%SetType(a%kfl_ReadType)
   call OTAReadOnNodes%InitializeOnNodes(a)
   call OTAReadOnNodes%Loop
   call OTAReadOnNodes%FinalizeOnNodes
end subroutine

subroutine php_OTAReadOnBoundariesDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllReadOnBoundaries
   class(PhysicalProblem) :: a
   type(OneToAllReadOnBoundaries) ::   OTAReadOnBoundaries   
   
   call OTAReadOnBoundaries%SetType(a%kfl_ReadType)
   call OTAReadOnBoundaries%Initialize(a)
   call OTAReadOnBoundaries%Loop
   call OTAReadOnBoundaries%Finalize
end subroutine

subroutine php_OTAReadOnElementsDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllReadOnElements
   class(PhysicalProblem) :: a
   type(OneToAllReadOnElements) ::   OTAReadOnElements   
   
   call OTAReadOnElements%SetType(a%kfl_ReadType)
   call OTAReadOnElements%Initialize(a)
   call OTAReadOnElements%Loop
   call OTAReadOnElements%Finalize
end subroutine

subroutine php_OTAReadSourceDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllReadSource
   class(PhysicalProblem) :: a
   type(OneToAllReadSource) ::   OTAReadSource   
   
   call OTAReadSource%SetType(a%kfl_ReadType)
   call OTAReadSource%InitializeSource(a)
   call OTAReadSource%Loop
   call OTAReadSource%FinalizeSource
end subroutine
