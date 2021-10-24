subroutine  plcd_Reabcs(a,itask,kflag)
   ! NAME
   !    plcd_reabcs
   ! DESCRIPTION
   !    This routine reads the specific boundary conditions for PLCD.
   !
   !    For conditions on nodes, bvess(idime,ipoin,1) must contain the
   !    velocity and temperature values multiplied by density.
   !    The different codes for kfl_fixno_nsi(ipoin) are:
   !    = 1 ... Dirichlet
   !    = 0 ... Free or initial
   !   For conditions on boundaries, bvnat(iboun) is allocated ONLY if
   !    a condition is imposed on it. Overmore, its length depends on the
   !    condition type. The different codes for kfl_fixbo(iboun) are:
   !    = 1 ... Dirichlet ............ u
   !    = 2 ... Traction ............. t
   !
   
   use typre
   use MPI
   use Mod_PLCD
   use Mod_plcd_Stages
   implicit none
  
   class(PLCDProblem), target :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag
   
   integer(ip) :: iboun,ipoin
   integer(ip) :: nboun,npoin,ndime
   
   character(150) :: outstr
   
   !Periodic boundary conditions
   integer(ip) :: kfl_perio,istage
   
   !MPI
   integer :: ierr
   
   !Stages pointer
   type(Stage), pointer :: s => NULL()
   
   !Initializations
   if (itask == 0) then
   
      a%kfl_confi = 0                                    ! Flow is not confined
      a%nodpr     = 0                                    ! Node where to impose pressure
   !Header   
   elseif(itask == 1) then
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         if (a%Listener%exists('FIXPR')) then
            a%kfl_confi =  1
            a%nodpr=a%Listener%getint('FIXPR',1,'#Node where to impose pressure')
            
            !If Periodic boundary conditions, then set the point to my master
            !(if I am a slave)
            call a%Mesh%GetPerio(kfl_perio)
            if (kfl_perio == 1) call a%Mesh%Initial2MasterInitial(a%nodpr,a%nodpr)
            
            call a%Mesh%Initial2Global(a%nodpr,a%nodpr)
         endif
      endif
      
      !Send information to all
      if (a%kfl_ReadType == 0) then
         CALL MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%nodpr, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      endif
      
      if (a%kfl_confi == 1) then
         if (a%nodpr == 0) then 
            a%nodpr = 1  !-1 means not chosen, 0 means chosen but not in my process
         endif
         call a%Mesh%Global2Local(a%nodpr,a%nodpr)
      endif   
      
    
    
    
    !Finalization
   elseif(itask == 100) then
   
      !Ghost communicate data for boundary conditions
      do istage = 1,a%NumberOfStages
         s => a%Stages(istage)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofbc,s%kfl_fixno)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofbc,s%bvess)
      enddo
   endif

end subroutine



subroutine plcd_ReadOnNodes(a)
   use typre
   use Mod_Mesh
   use Mod_MPIObject
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_PLCD
   use Mod_plcd_Stages
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   class(PLCDProblem), target :: a
   
   integer(ip)::ndime,idime,ipoin,ifun
   type(Stage), pointer :: s => NULL()

   !Stage change
   if (a%Listener%param(1) == -1) then
      a%CurrentReadingStage =>  a%Stages(int(a%Listener%param(2)))
   else
      s => a%CurrentReadingStage
   
      !Read point number
      ipoin                   = int(a%Listener%param(1)) 
      !To local numbering
      call a%Mesh%Global2Local(ipoin,ipoin)  
      !Read nodal Fixity conditions
      s%kfl_fixno(1,ipoin)   = int(a%Listener%param(2))
      s%bvess(1:a%ndofbc,ipoin)=a%Listener%param(3:2+a%ndofbc)
      call codfix(a%ndofbc,s%kfl_fixno(1,ipoin))
      
      !Boundary conditions are updated here
      if (a%kfl_StochasticTopologyOptimization == 1) then
         !Functions valued -1 are reserved for plcd_stochastic topology optimization
         ifun = int(a%Listener%param(2+a%ndofbc+1))
         call plcd_TDSTO_modifybcs(a,ipoin,s%bvess(1:a%ndofbc,ipoin),ifun)
         
      endif
   endif
end subroutine

subroutine plcd_ReadSource(a)
   use typre
   use Mod_Mesh
   use Mod_MPIObject
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_PLCD
   use Mod_plcd_Stages
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   class(PLCDProblem),target :: a
   
   integer(ip) :: idime,ndime,ipoin,ifun
   type(Stage), pointer :: s => NULL()
   
   !Stage change
   if (a%Listener%param(1) == -1) then
      a%CurrentReadingStage =>  a%Stages(int(a%Listener%param(2)))
   else
      s => a%CurrentReadingStage
      ipoin                  = int(a%Listener%param(1)) 
      call a%Mesh%Global2Local(ipoin,ipoin)

      s%NodalForces(:,ipoin) = a%Listener%param(2:1+size(s%NodalForces,1))
      
      
      

      !Boundary conditions are updated here
      if (a%kfl_StochasticTopologyOptimization == 1) then
         !Functions valued -1 are reserved for plcd_stochastic topology optimization
         ifun = int(a%Listener%param(2+size(s%NodalForces,1)))
         call plcd_TDSTO_modifybcs(a,ipoin,s%NodalForces(:,ipoin),ifun)
         
      endif
      
   endif
   
   
end subroutine



subroutine plcd_ReadOnBoundaries(a)
   use typre
   use Mod_Mesh
   use Mod_MPIObject
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_PLCD
   use Mod_plcd_Stages
   !use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   class(PLCDProblem),target :: a
   
   integer(ip), allocatable :: knodb(:)
   
   integer(ip) :: pnodb
   integer(ip) :: iboun,idime,inodb,ipsta,ndofi
   
   type(Stage), pointer :: s => NULL()
   
   !Stage change
   if (a%Listener%param(1) == -1) then
      a%CurrentReadingStage =>  a%Stages(int(a%Listener%param(2)))
   else
      s => a%CurrentReadingStage
      
      !Read boundary nodes
      pnodb=int(a%Listener%param(2))
      allocate(knodb(pnodb))
      !call a%Memor%alloc(pnodb,knodb,'knodb','plcd_reabcs')
      knodb(1:pnodb)=int(a%Listener%param(3:2+pnodb))
      
      !To local numbering
      call a%Mesh%Global2Local(pnodb,knodb(1:pnodb),knodb(1:pnodb))

      !Find which is the boundary
      call a%Mesh%Finbou(pnodb,knodb(1:pnodb),iboun)    
      if(iboun==0) call runend('php_ReadOnBoundaries: Boundary not found')
      
      ipsta=3+pnodb
      s%kfl_fixbo(iboun) = int(a%Listener%param(ipsta))
      if ((a%kfl_fixbo(iboun)>=0).AND.(a%kfl_fixbo(iboun) /= s%kfl_fixbo(iboun) )) call runend('php_ReadOnBoundaries: tried to re-assign a different condition to boundary')
      
      ! Dirichlet
      if(s%kfl_fixbo(iboun)==1) then
         if (associated(s%bvnat(iboun)%a)) then
            s%bvnat_coun = s%bvnat_coun - size(s%bvnat(iboun)%a)
            deallocate(s%bvnat(iboun)%a)
         endif
      
         allocate(s%bvnat(iboun)%a(a%ndofbc))
         ndofi=0
            do idime=1,a%ndofbc
               ipsta=ipsta+1
               ndofi=ndofi+1
               s%bvnat(iboun)%a(ndofi)=a%Listener%param(ipsta)
            end do
            
         s%kfl_funbo(iboun) = int(a%Listener%param(ipsta+1))
         
      !Neumann Boundary condition
      else if(s%kfl_fixbo(iboun)==2) then
         if (associated(s%bvnat(iboun)%a)) then
            s%bvnat_coun = s%bvnat_coun - size(s%bvnat(iboun)%a)
            deallocate(s%bvnat(iboun)%a)
         endif
   
         allocate(s%bvnat(iboun)%a(a%ndofbc))
         ndofi=0
            do idime=1,a%ndofbc
               ipsta=ipsta+1
               ndofi=ndofi+1
               s%bvnat(iboun)%a(ndofi)=a%Listener%param(ipsta)
            end do
      end if
      
      if (associated(s%bvnat(iboun)%a)) then
      s%bvnat_coun = s%bvnat_coun + size(s%bvnat(iboun)%a)
      endif
   endif
         
end subroutine

subroutine plcd_ReadOnElements(a)
   use typre
   use Mod_Listen
   use Mod_PLCD
   implicit none
   class(PLCDProblem) :: a

   integer(ip) :: ielem,imaterial,pnode,pgaus,ElementInfoInipos,angleSize,ndime
   
   
   
   !a%AuxRead_ElementMaterialList(a%gielem) = int(a%Listener%param(a%gipsta+1))
   
   ielem = a%gielem
   imaterial = int(a%Listener%param(a%gipsta+1))
   call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
   if (imaterial <= 0 .or. imaterial > a%NumberOfMaterials) call runend('Wrong material number')
   call a%Materials(imaterial)%p%CreateElementData(pgaus,a%ElementMaterialsData(ielem)%p)
   
   ElementInfoInipos = a%gipsta+2
   
   if (a%kfl_ElementRotators == 1) then
      call a%Mesh%GetNdime(ndime)
      if (ndime == 2) then
         angleSize = 1
      else
         anglesize = 3
      endif
      
      
      call a%ElementMaterialsData(ielem)%p%AllocateRotator(ndime,a%Listener%param(ElementInfoInipos:ElementInfoInipos+anglesize-1),a%Memor)
      
      if (ndime == 2) then 
         ElementInfoInipos = ElementInfoInipos+1
      else
         ElementInfoInipos = ElementInfoInipos+3
      endif
   endif
   
   call a%ElementMaterialsData(ielem)%p%SpecificReadElementData(a%Listener%param(ElementInfoInipos))
   

end subroutine


subroutine plcd_OTAReadOnNodesDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllPLCDReadOnNodes
   class(PhysicalProblem) :: a
   type(OneToAllPLCDReadOnNodes) ::   OTAReadOnNodes   
   
   call OTAReadOnNodes%SetType(a%kfl_ReadType)
   call OTAReadOnNodes%InitializeOnNodes(a)
   call OTAReadOnNodes%Loop
   call OTAReadOnNodes%FinalizeOnNodes
end subroutine

subroutine plcd_OTAReadOnBoundariesDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllPLCDReadOnBoundaries
   class(PhysicalProblem) :: a
   type(OneToAllPLCDReadOnBoundaries) ::   OTAReadOnBoundaries   
   
   call OTAReadOnBoundaries%SetType(a%kfl_ReadType)
   call OTAReadOnBoundaries%Initialize(a)
   call OTAReadOnBoundaries%Loop
   call OTAReadOnBoundaries%Finalize
end subroutine

subroutine plcd_OTAReadSourceDo(a)
   use Mod_PhysicalProblem
   use Mod_OneToAllPLCDReadSource
   class(PhysicalProblem) :: a
   type(OneToAllPLCDReadSource) ::   OTAReadSource   
   
   call OTAReadSource%SetType(a%kfl_ReadType)
   call OTAReadSource%InitializeSource(a)
   call OTAReadSource%Loop
   call OTAReadSource%FinalizeSource
end subroutine


subroutine plcd_ReadOnFunctions(a,ifunc)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem), target :: a
   integer(ip) :: ifunc
   
   if (a%Listener%words(2) == 'STOCH') then
      a%kfl_funty(ifunc,1)=-1
      a%kfl_funty(ifunc,2)=2
      !funpa(1) :: Module Variance
      !funpa(2) :: Direction angle variance
   elseif (a%Listener%words(2) == 'MULTI') then   
      a%kfl_funty(ifunc,1)=-2
      a%kfl_funty(ifunc,2)=2
      !funpa(1) :: Module Variance
      !funpa(2) :: Direction angle variance
      
   elseif (a%Listener%words(2) == 'STO3Z') then   
      a%kfl_funty(ifunc,1)=-3
      a%kfl_funty(ifunc,2)=2
      !funpa(1) :: nothing
      !funpa(2) :: angle variance
      
   endif   
end subroutine

subroutine plcd_ScatterOnFunctions(a)
   use typre
   use Mod_PLCD
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   class(PLCDProblem) :: a
   
   !Initializations for stochastic optimization   
   !Needs to be done here, after functions for boundary conditions are read
   !but before nodes boundary conditions are read :(
   if (a%kfl_StochasticTopologyOptimization == 1) then
      call plcd_TDSTO_turnon(a)
   endif
   
   
end subroutine
   