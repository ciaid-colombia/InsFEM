subroutine php_begite(a)
!-----------------------------------------------------------------------
! NAME 
!    php_begite
! DESCRIPTION
!    This routine starts an internal iteration
!-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   !Initializations
   a%kfl_goite = 1
   a%itera     = 0
   
   if(a%kfl_docoupconv) call a%Updbcs
   call a%SpecificBegite
   
end subroutine php_begite
