subroutine sldale_solite(a)
   !-----------------------------------------------------------------------
   !> This is the main routine of the linear system of mesh displacements.
   !! The ALE flag is set to 0 so the displacement field is computed in the
   !! original coordinates. The currentbvess variable allows the decoupling 
   !! of the system for each dimension of the mesh. Finally, the ALE flag is
   !! reset to 1 so the computed displacements can be added to the coordinates
   !! when loading the mesh, as well as to update the mass matrix and the 
   !! exterior normal.
   !-----------------------------------------------------------------------
   use typre
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem),target :: a
   integer(ip)  :: currentbvess

   !We compute the ALE movement with the deformed mesh
   call a%Mesh%SetALE(0_ip)

   a%itera = a%itera + 1
   if (a%MPIrank == a%MPIroot) then
       if(a%itera==1) write(a%lun_solve,100) a%istep
       write(a%lun_solve,101) a%itera      
   endif

   call a%Timer%BuildLinearSystem%Tic
   !Construct the system matrix and right-hand-side.
   ! Initializations
   call a%LinearSystem%ToZero

   call a%sldale_elmope(currentbvess)
   call a%Timer%BuildLinearSystem%Toc

   !Solve the algebraic system.
   call a%Timer%SolveLinearSystem%Tic

   a%unkno = a%solid%unkno

   !Ghost communicates
   call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%unkno)

   call a%Timer%SolveLinearSystem%Toc

   !!Formats. 
   100 format(/,'SOLVER INFORMATION FOR ISTEP: ',i5)
   101 format('------------------------------------------------------------', &
           /,'   INNER ITERATION NUMBER: ',i5)


end subroutine sldale_solite
