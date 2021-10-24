subroutine  ale_Reabcs(a,itask,kflag)
   ! NAME
   !    ale_reabcs
   ! DESCRIPTION
   !    This routine reads the specific boundary conditions for NSTINC.
   !
   !    For conditions on nodes, bvess(idime,ipoin,1) contains the
   !    velocity value.
   !    The different codes for kfl_fixno_nsi(ipoin) are:
   !    = 1 ... Dirichlet
   !    = 0 ... Free Surface
   !    = 2 ... Not Fixed: we transform it to -1
   !
   !   
   !------------------------------------------------------------------------
   use MPI
   use typre
   use Mod_AleMov
   implicit none
  
   class(AlemovProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag
   
   integer(ip) :: iboun,ipoin,idime,ndime
   integer(ip) :: nboun,npoin
   
   character(150) :: outstr
   
   !Periodic boundary conditions
   integer(ip) :: kfl_perio
   
   !MPI
   integer :: ierr

   !Initializations
   if (itask == 0) then
      
   !Header   
   elseif(itask == 1) then
    
    !Finalization
   elseif(itask == 100) then
      
     !We transform fixno 2 to -1
     !For legacy reasons (should be 0 but it was occupied by Free surface (apont))
     call a%Mesh%Getnpoin(npoin) 
     call a%Mesh%GetNdime(ndime)
     do ipoin = 1,npoin
        do idime = 1,ndime    
           if (a%kfl_fixno(idime,ipoin) == 2) a%kfl_fixno(idime,ipoin) = -1
        enddo
     enddo

   endif
    

end subroutine

subroutine ale_ReadOnNodes(a)
   use typre
   use Mod_AleMov
   implicit none
   class(AlemovProblem) :: a
   
   integer(ip) :: ndime
   
end subroutine


subroutine ale_ReadOnFunctions(a,ifunc)
   use typre
   use Mod_AleMov
   implicit none
   class(AlemovProblem) :: a
   integer(ip) :: ifunc
   
   if (a%Listener%words(2) == 'ROTAT') then
      a%kfl_funty(ifunc,1)=-1
      a%kfl_funty(ifunc,2)=6
   endif   
end subroutine

subroutine ale_ReadOnBoundaries(a)
   use typre
   use Mod_Alemov 
   implicit none
   class(AlemovProblem) :: a
   integer(ip)                :: nbody,iboun,nboun,pnodb,ipsta

end subroutine
