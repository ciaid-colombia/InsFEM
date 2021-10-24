subroutine nsc_begste(a)
   ! DESCRIPTION
   !    This routine prepares for a new time step of the Compressible NS
   !    equations.
   !-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   
!   !!Check if the coupling between modules is done
!   call nsi_coupling_outerr(a)
   
   if(a%kfl_timei==1) then
      a%densf(:,2) = a%densf(:,3)
      a%momen(:,:,2) = a%momen(:,:,3)
      a%energ(:,2) = a%energ(:,3)
      a%densf(:,1) = a%densf(:,2)
      a%momen(:,:,1) = a%momen(:,:,2)
      a%energ(:,1) = a%energ(:,2)
   end if
   

!Dirichlet Boundary Conditions
!call nsc_elmdir(a,e,elmat,elrhs)
!elmope residual must be zero for the dirichlet values
!dirichlet values are fixed at the initial time if steady
!if unsteady, then dirichlet values are changed at
!the begining of each time step
   
   call a%SpecificNSCompBegste

   !Check if flow is confined (only if non-constant boundary conditions)
   if (a%kfl_conbc /= 1) call a%Ifconf

end subroutine nsc_begste
