subroutine opt_getste(a,dtinv)
!-----------------------------------------------------------------------
! NAME 
!    opt_getste
! DESCRIPTION
!    This routine computes the time step size for the incompressible NS
!    equation.
!-----------------------------------------------------------------------
   use typre
   use Mod_Optics
   implicit none 
   class(OpticsProblem) :: a
   real(rp) :: dtinv
   
   
   
   call a%Timer%Total%Tic
   call a%Timer%Getste%Tic

   !Critical time step
   dtinv = 0
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
end subroutine opt_getste
