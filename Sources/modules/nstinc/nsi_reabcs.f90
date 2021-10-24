subroutine  nsi_Reabcs(a,itask,kflag)
   ! NAME
   !    nsi_reabcs
   ! DESCRIPTION
   !    This routine reads the specific boundary conditions for NSTINC.
   !
   !    For conditions on nodes, bvess(idime,ipoin,1) contains the
   !    velocity value.
   !    The different codes for kfl_fixno_nsi(ipoin) are:
   !    = 1 ... Dirichlet
   !    = 0 ... Free or initial
   !    = 3 ... Disconnect elements .. (dvol=0)
   !    = 5 ... Only Enforce Dirichlet boundary conditions if compressed fluid (p > 0). For Free surface problems
   !    = 6 ....Elastic boundary
   !
   !    For conditions on boundaries, bvnat(iboun) is allocated ONLY if
   !    a condition is imposed on it. Overmore, its length depends on the
   !    condition type. The different codes for kfl_fixbo(iboun) are:
   !    = 1 ... Dirichlet ............ u
   !    = 2 ... Pressure imposed  .... sig.n=-p n
   !    = 3 ... Wall law ............. sig.n.t=-rho*U_*^2*u/|u|
   !    = 4 ... Symmetry/Slip wall ... sig.n.t=0, u.n=0
   !    = 5 ... Dynamic pressure ..... sig.n=-1/2*rho*u^2
   !    = 6 ... Open flow ............ Assemble 2*mu*Sym(grad(u).n
   !    = 7 ... No slip wall ......... u=0
   !    At the end of the subroutine conditions on boundaries are
   !    transfered to conditions on nodes.
   !------------------------------------------------------------------------
   use MPI
   use typre
   use Mod_NavierStokes
   implicit none
  
   class(NavierStokesProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag
   
   integer(ip) :: iboun,ipoin,ndime
   integer(ip) :: nboun,npoin
   
   character(150) :: outstr
   
   !Periodic boundary conditions
   integer(ip) :: kfl_perio
   
   !MPI
   integer :: ierr
   
   integer(ip) :: iaux, idime
   
   real(rp) :: delta
   
   interface 
      subroutine nsi_BCNoSlipToNodes(a,iboun)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: iboun
      end subroutine
   end interface

   !Initializations
   if (itask == 0) then
      
      !Output
      outstr = adjustl(trim(a%exmod))//'_SPECIFIC_REABCS'
      
      a%kfl_inist = 0                                    ! First step is stokes
      a%kfl_confi = 0                                    ! Flow is not confined
      a%nodpr     = 0                                    ! Node where to impose pressure
      a%kfl_local = 0                                    ! No local system of reference
      a%kfl_ExchangeLocationWallLaw = 0                  ! Default is use a delta for the wall law, not the exchange location
      a%kfl_FixWallLawDeltaForAll = 0                    ! Default is each boundary can have its own delta
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNboun(nboun)
      call a%Mesh%GetNdime(ndime)
      !Arrays for conditions on boundaries
      call a%Memor%alloc(npoin,a%kfl_fixrs,'kfl_fixrs','nsi_memall')
      call a%Memor%alloc(nboun,a%kfl_bours,'kfl_bours', 'nsi_memall')
      call a%Memor%alloc(ndime + 1,npoin,a%PointwiseSource,'PointwiseSource','nsi_reabcs')
      
   !Header   
   elseif(itask == 1) then
   
      !Root Reads data
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         if(a%Listener%exists('FIXPR')) then
            a%kfl_confi =  1
            a%nodpr=a%Listener%getint('FIXPR',1,'#Node where to impose pressure')
            
            !If Periodic boundary conditions, then set the point to my master
            !(if I am a slave)
            call a%Mesh%GetPerio(kfl_perio)
            if (kfl_perio == 1) call a%Mesh%Initial2MasterInitial(a%nodpr,a%nodpr)
            
            call a%Mesh%Initial2Global(a%nodpr,a%nodpr)
         else if(a%Listener%exists('DONTF')) then
            a%kfl_confi = -1
         end if
         if(a%Listener%exists('STOKE')) then
            a%kfl_inist=1
         end if
         if(a%Listener%exists('EXCHA')) then
            !Exhange Location Wall law
            a%kfl_ExchangeLocationWallLaw = 1
         endif
         if(a%Listener%exists('FIXWA')) then
            !Fix Wall law delta value for all elements
            a%kfl_FixWallLawDeltaForAll = 1
            a%WallLawDeltaValue=a%Listener%getrea('FIXWA',0.001_rp,'#Wall law delta value')
         endif
         
      endif
      
      !Send information to all
      if (a%kfl_ReadType == 0) then
         CALL MPI_BCAST(a%kfl_inist, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%nodpr, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%kfl_ExchangeLocationWallLaw, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%kfl_FixWallLawDeltaForAll, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%WallLawDeltaValue, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      endif
      
      if (a%nodpr == 0) then 
         a%nodpr = -1  !-1 means not chosen, 0 means chosen but not in my process
      else   
         call a%Mesh%Global2Local(a%nodpr,a%nodpr)
      endif   
    
    !Finalization
   elseif(itask == 100) then
      
      !-------------------------------------------------------------------------
      !COMMUNICATE NODAL GHOST DATA WITH NEIGHBOURS
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%kfl_fixrs)
      
      !Final operations
      call a%Mesh%Getnboun(nboun)
      call a%Mesh%Getnpoin(npoin)
      
      !Wall law, change to slip or no slip if negative or too small
      do iboun=1,nboun
         if(a%kfl_fixbo(iboun)==3) then
            if (a%kfl_FixWallLawDeltaForAll == 0) then
               delta = a%bvnat(iboun)%a(1)
            elseif (a%kfl_FixWallLawDeltaForAll == 1) then
               delta = a%WallLawDeltaValue
            endif
            !Negative value for delta, switch to slip
            if (delta < -zensi) then
               a%kfl_fixbo(iboun)=4
            !Too small value for delta, switch to noslip   
            elseif( delta < zensi) then
               a%kfl_fixbo(iboun)=7
               call nsi_BCNoSlipToNodes(a,iboun)
            endif
         endif
      end do

     
      !Check if there is a local prescription or suction Dirichlet Boundary conditions
      a%kfl_ElasticBoundary = 0
      call a%Mesh%GetNdime(ndime)
      nodes: do ipoin=1,npoin
         if(a%kfl_fixrs(ipoin)/=0) then
            a%kfl_local=1
         end if
         do idime = 1,ndime
            if (a%kfl_fixno(idime,ipoin) == 5) then
               a%kfl_SuctionDirichletBC = 1
            elseif (a%kfl_fixno(idime,ipoin) == 6) then
               a%kfl_ElasticBoundary = 1
            endif
         enddo
      end do nodes
      
      !Reduce the info so that it is available to all processors
      call  MPI_AllREDUCE( a%kfl_local, iaux, 1, MPI_INTEGER4, MPI_MAX,a%MPIcomm, ierr )
      a%kfl_local = iaux
      
      call  MPI_AllREDUCE( a%kfl_SuctionDirichletBC, iaux, 1, MPI_INTEGER4, MPI_MAX,a%MPIcomm, ierr )
      a%kfl_SuctionDirichletBC = iaux
      
      call  MPI_AllREDUCE( a%kfl_ElasticBoundary, iaux, 1, MPI_INTEGER4, MPI_MAX,a%MPIcomm, ierr )
      a%kfl_ElasticBoundary = iaux
         
      if (a%kfl_ElasticBoundary == 1) then
         call a%Memor%alloc(ndime,npoin,a%EB_Displacement,'EB_Displacement','nsi_memall')
      endif
   endif
end subroutine

subroutine nsi_ReadOnNodes(a)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   logical                    :: isALE
   integer(ip)                :: ipoin,npoin
   
   !Reads the reference system
   a%kfl_fixrs(a%gipoin) = int(a%Listener%param(4+a%ndofbc))

   call a%Mesh%GetALE(isALE)   
   if (isALE) a%bvess(:,a%gipoin,2)=a%bvess(:,a%gipoin,1)

   if(a%kfl_fixno(1,a%gipoin)==3) then
      a%kfl_fixno(:,a%gipoin)=-3
   elseif(a%kfl_fixno(1,a%gipoin)==8) then
      a%kfl_fixno(:,a%gipoin)=-8
   end if
end subroutine

subroutine nsi_ReadSource(a)
   use typre
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip)::ndime,idime,localNode

   call a%Mesh%GetNdime(ndime)

   do idime=0,ndime
       a%PointwiseSource(idime+1,a%gipoin) = a%Listener%param(idime+1)
   enddo
end subroutine

subroutine nsi_ReadOnBoundaries(a)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip)                :: nbody,iboun,nboun,pnodb,ipsta
   
   !Dirichlet
   if(a%kfl_fixbo(a%giboun) == 1) then
      if(a%kfl_conbc==1) then
         a%kfl_bours(a%giboun) = int(a%Listener%param(a%gipsta+1))
      else
         a%kfl_bours(a%giboun) = int(a%Listener%param(a%gipsta+2))
      endif
      
   !Neumann  
   elseif (a%kfl_fixbo(a%giboun) == 2) then
      
   !Wall law  
   elseif (a%kfl_fixbo(a%giboun) == 3) then
      allocate(a%bvnat(a%giboun)%a(1))
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
   
   !Slip wall/symmetry  
   elseif (a%kfl_fixbo(a%giboun) == 4) then
               
   !Dynamic pressure: sig.n=-(-1/2*rho*u^2) n
   elseif (a%kfl_fixbo(a%giboun) == 5) then
   
   !Open boundary: assemble 2*mu*Sym(grad(u).n
   elseif (a%kfl_fixbo(a%giboun) == 6) then
   
   !No slip wall 
   elseif (a%kfl_fixbo(a%giboun) == 7) then
   
   !Atmospheric stress + pressure imposed
   else if(a%kfl_fixbo(a%giboun) == 8) then
      allocate(a%bvnat(a%giboun)%a(1))
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
     
   !Prescribed traction: sig.n = tract_0 (tract_0 given)
   else if(a%kfl_fixbo(a%giboun) == 9) then
      allocate(a%bvnat(a%giboun)%a(3))
      a%bvnat(a%giboun)%a(1:3) = a%Listener%param(a%gipsta+1:a%gipsta+3)
   endif
end subroutine

subroutine nsi_ReadOnFunctions(a,ifunc)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip) :: ifunc
   
   if (a%Listener%words(2) == 'SPANW') then
      a%kfl_funty(ifunc,1)=-1
      a%kfl_funty(ifunc,2)=3
   elseif (a%Listener%words(2) == 'BLASI') then
      a%kfl_funty(ifunc,1) = -2
      a%kfl_funty(ifunc,2) = 3
   elseif (a%Listener%words(2) == 'TURBU') then
      a%kfl_funty(ifunc,1) = -3
      a%kfl_funty(ifunc,2) = 0
   endif   
end subroutine


