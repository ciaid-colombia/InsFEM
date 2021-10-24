subroutine  nsc_Reabcs(a,itask)
   ! NAME
   !    nsc_reabcs
   ! DESCRIPTION
   !    This routine reads the specific boundary conditions for COMPRE.
   !
   !    For conditions on nodes, bvess(idime,ipoin,1) must contain the 
   !    density, velocity and temperature values for the mass, momentum
   !    and Total energy equations.
   !    The different codes for kfl_fixno(idofbc,ipoin) are:
   !    = 0 ... Free(Neumann null) or initial
   !    = 1 ... Dirichlet
   !
   !    For conditions on boundaries, bvnat(iboun) is allocated ONLY if 
   !    Neumann or Robin conditions are imposed on it. 
   !    Overmore, its length depends on the condition type. 
   !    The different codes for kfl_fixbo(iboun) are:
   !          Mas  | Mom | Ene
   !    = 1   Dir  | Dir | Dir ... rho  | rho u       | rho E
   !    = 2   Fre  | Fre | Fre ... Free | Free        | Free
   !    = 3   Fre  | Neu | Fre ... Free | sigma.n=val | Free
   !    = 4   Fre  | Neu | Neu ... Free | sigma.n=val | q.n=val
   !    = 5   Fre  | Neu | Rob ... Free | sigma.n=val | a(ti-t).n=val
   !    = 6   Fre  | Fre | Neu ... Free | Free        | q.n=val      
   !    = 7   Fre  | Fre | Rob ... Free | Free        | a(ti-t).n=val
   !    = 8   Fre  | Wlaw| Fre ... Free | sig.n.t=-rho*U_*^2*u/|u| | Free 
   !    Mass can only posses Dirichlet or Free boundary conditions.
   !------------------------------------------------------------------------
   use typre
   use MPI
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip) :: itask
   integer(ip) :: kfl_perio !Periodic boundary conditions
   character(150) :: outstr
   integer :: ierr
   

   !Initializations
   if (itask == 0) then
            
      !Output
      outstr = adjustl(trim(a%exmod))//'_SPECIFIC_REABCS'
      a%delta     = 0.0_rp                               ! Distance to the wall

      a%kfl_confi = 0                                    ! Flow is not confined
      a%kfl_nscnbc = 0
      
      call a%SpecificNSCompReabcs(0)

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
      endif
      
      !Root Reads data
      if (a%MPIrank == a%MPIroot .or. a%kfl_ReadType == 1) then
         if(a%Listener%exists('WALLD')) then
            a%delta=a%Listener%getrea('WALLD',0.0_rp,'#Distance to the wall')
         end if
      endif
      
      !Send information to all
      if (a%kfl_ReadType == 0) then
         CALL MPI_BCAST(a%delta, 1, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%nodpr, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
         CALL MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      endif
      
      call a%SpecificNSCompReabcs(1)

    !Finalization
   elseif(itask == 100) then

      !Send information to all
      if (a%kfl_ReadType == 0) then
         CALL MPI_BCAST(a%kfl_nscnbc, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      endif

      if (a%kfl_nscnbc == 1) then
         call a%CalculateBoundaryVmass
      endif

      call a%SpecificNSCompReabcs(100)

   endif
    
   

end subroutine



subroutine nsc_ReadOnNodes(a)
   use typre
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   logical                    :: isALE
   
   call a%Mesh%GetALE(isALE)   
   if (isALE) a%bvess(:,a%gipoin,2)=a%bvess(:,a%gipoin,1)
   
end subroutine

subroutine nsc_ReadOnBoundaries(a)
   use typre
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   
   integer(ip)                :: idime,ndime

   call a%Mesh%GetNdime(ndime)

   ! Free | Neumann | Free  
   if (a%kfl_fixbo(a%giboun) == 3) then
      allocate(a%bvnat(a%giboun)%a(3))
      !Momentum Neumann value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1) 
      a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2) 
      a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3) 
   
   ! Free | Neumann   | Neumann    
   elseif (a%kfl_fixbo(a%giboun) == 4) then
      allocate(a%bvnat(a%giboun)%a(4))
      !Momentum Neumann value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1) 
      a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2) 
      a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3) 
      !Energy Neumann value
      a%bvnat(a%giboun)%a(4)=a%Listener%param(a%gipsta+4)
               
   ! Free | Neumann   | Robin  
   elseif (a%kfl_fixbo(a%giboun) == 5) then
      allocate(a%bvnat(a%giboun)%a(6))
      !Momentum Neumann value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1) 
      a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2) 
      a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3) 
      !Energy Neumann value
      a%bvnat(a%giboun)%a(4)=a%Listener%param(a%gipsta+4)
      !Temperature value 
      a%bvnat(a%giboun)%a(5)=a%Listener%param(a%gipsta+5)
      !Convective value
      a%bvnat(a%giboun)%a(6)=a%Listener%param(a%gipsta+6)
 
   !Free | Free | Neumann
   elseif (a%kfl_fixbo(a%giboun) == 6) then
      allocate(a%bvnat(a%giboun)%a(1))
      !Energy Neumann value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
   
   !Free | Free | Robin
   elseif (a%kfl_fixbo(a%giboun) == 7) then
      allocate(a%bvnat(a%giboun)%a(3))
      !Energy Neumann value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
      !Temperature value 
      a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2)
      !Convective value
      a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3)

   !Free | Wall_law | Free
   elseif (a%kfl_fixbo(a%giboun) == 8) then 

   !Free | Backflow penalty | Free
   elseif (a%kfl_fixbo(a%giboun) == 9) then
      allocate(a%bvnat(a%giboun)%a(1))
      !Penalty treshold value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)

   !Non Reflecting boundary condition
   elseif (a%kfl_fixbo(a%giboun) == 10) then
      a%kfl_nscnbc = 1
      allocate(a%bvnat(a%giboun)%a(3))
      !NSCBC sigma value
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
      !NSCBC char lenght value
      a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2)
      !NSCBC unperturbed pressure
      a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3)

   !Weak NSCBC imposition
   elseif (a%kfl_fixbo(a%giboun) == 11) then
      allocate(a%bvnat(a%giboun)%a(4))
      !Penalty coefficient
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
      !NSCBC sigma value
      a%bvnat(a%giboun)%a(2)=a%Listener%param(a%gipsta+2)
      !NSCBC char lenght value
      a%bvnat(a%giboun)%a(3)=a%Listener%param(a%gipsta+3)
      !NSCBC unperturbed pressure
      a%bvnat(a%giboun)%a(4)=a%Listener%param(a%gipsta+4)

   !Non Reflecting Pressure boundary condition
   elseif (a%kfl_fixbo(a%giboun) == 12) then
      a%kfl_nscnbc = 1

   !Weak NSCBC Non Reflecting pressure
   elseif (a%kfl_fixbo(a%giboun) == 13) then
      allocate(a%bvnat(a%giboun)%a(1))
      !Penalty coefficient
      a%bvnat(a%giboun)%a(1)=a%Listener%param(a%gipsta+1)
   end if

end subroutine


