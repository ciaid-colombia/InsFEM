subroutine  lmn_Reabcs(a,itask,kflag)
   ! NAME
   !    lmn_reabcs
   ! DESCRIPTION
   !    This routine reads the specific boundary conditions for LMACHN.
   !
   !    For conditions on nodes, bvess(idime,ipoin,1) must contain the
   !    velocity and temperature values multiplied by density.
   !    The different codes for kfl_fixno_nsi(ipoin) are:
   !    = 1 ... Dirichlet
   !    = 0 ... Free or initial
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
   use Mod_LowMach
   implicit none
  
   class(LowMachProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag
   integer(ip) :: iboun,ipoin,nboun,npoin
   character(150) :: outstr
   integer(ip) :: kfl_perio !Periodic boundary conditions
   integer :: ierr
   

   !Initializations
   if (itask == 0) then
      
      !Output
      outstr = adjustl(trim(a%exmod))//'_SPECIFIC_REABCS'
      
      a%kfl_confi = 0                                    ! Flow is not confined
      a%nodpr     = 0                                    ! Node where to impose pressure
      a%delta     = 0.0_rp                               ! Distance to the wall
      
      !Get Mesh data
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNboun(nboun)

      if(a%kfl_sourc==2) then
         call a%Mesh%GetNpoin(npoin)
         call a%Memor%alloc(npoin,a%PointwiseSource,'PointwiseSource','lmn_reabcs')
      end if
   
   !Header   
   elseif(itask == 1) then
      
      !Root Reads data
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         if(a%Listener%exists('FIXPR')) then
            a%nodpr=a%Listener%getint('FIXPR',1,'#Node where to impose pressure')
            
            !If Periodic boundary conditions, then set the point to my master
            !(if I am a slave)
            call a%Mesh%GetPerio(kfl_perio)
            if (kfl_perio == 1) call a%Mesh%Initial2MasterInitial(a%nodpr,a%nodpr)
            
            call a%Mesh%Initial2Global(a%nodpr,a%nodpr)
         else if(a%Listener%exists('DONTF')) then
            a%kfl_confi = -1
         end if
         if(a%Listener%exists('WALLD')) then
            a%delta=a%Listener%getrea('WALLD',0.0_rp,'#Distance to the wall')
         end if
      endif
      
      !Send information to all
      if (a%kfl_ReadType == 0) then
         CALL MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%nodpr, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%delta, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr) 
      endif
      
      if (a%nodpr == 0) then 
         a%nodpr = -1  !-1 means not chosen, 0 means chosen but not in my process
      else   
         call a%Mesh%Global2Local(a%nodpr,a%nodpr)
      endif   
 
 
    !Finalization
   elseif(itask == 100) then
      
      !Final operations
      call a%Mesh%Getnboun(nboun)
      call a%Mesh%Getnpoin(npoin)
      
      !Put slip wall condition to wall law boundaries if delta is negative
      if(a%delta>zelmn) then
      else if(a%delta<-zelmn) then
      if (a%MPIrank == a%MPIroot) write(a%lun_outph,*)&
                     'WALL DISTANCE IS NEGATIVE: WALL LAW CONDITION IS REPLACED BY A SLIP CONDITION'
         do iboun=1,nboun
            if(a%kfl_fixbo(iboun)==3) a%kfl_fixbo(iboun)=4
         end do
      else if( abs(a%delta) < zelmn) then
      if (a%MPIrank == a%MPIroot) write(a%lun_outph,*)&
                     'WALL DISTANCE IS VERY SMALL: WALL LAW CONDITION IS REPLACED BY A NON-SLIP CONDITION'
         do iboun=1,nboun
            if(a%kfl_fixbo(iboun)==3) a%kfl_fixbo(iboun)=7
         end do
      end if
 
   endif

end subroutine

subroutine lmn_ReadOnNodes(a)
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_Listen
   use Mod_PhysicalProblem
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   integer(ip) :: ndime
   
   call a%Mesh%GetNdime(ndime)
   if(a%kfl_sourc==2) a%PointwiseSource(a%gipoin) = a%Listener%param(4+ndime)

end subroutine


subroutine lmn_ReadOnBoundaries(a)
   use typre
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   
   !Neumann  
   if (a%kfl_fixbo(a%giboun) == 2) then
      
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
   
   ! Atmospheric stress + pressure imposed
   else if(a%kfl_fixbo(a%giboun) == 8) then
      allocate(a%bvnat(a%giboun)%a(1))
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)

   !Prescribed traction: sig.n = tract_0 (tract_0 given)
   else if(a%kfl_fixbo(a%giboun) == 9) then
      allocate(a%bvnat(a%giboun)%a(3))
      a%bvnat(a%giboun)%a(1:3) = a%Listener%param(a%gipsta+1:a%gipsta+3)
   endif

end subroutine


