subroutine nsc_begite(a)
!-----------------------------------------------------------------------
! DESCRIPTION
!    This routine starts an internal iteration for the compressible NS
!    equations. 
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   
   !Assign var(n,i,0) <-- var(n,i-1,*), initial guess for inner iterations
   a%densf(:,1) = a%densf(:,2)
   a%momen(:,:,1) = a%momen(:,:,2)
   a%energ(:,1) = a%energ(:,2)

   call a%SpecificNSCompBegite

end subroutine nsc_begite
