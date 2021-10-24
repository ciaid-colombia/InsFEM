subroutine GetCheckFoldedALE(a,AreFoldedALE)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   logical :: AreFoldedALE
   
   AreFoldedALE = a%AreFoldedALE
   
end subroutine

subroutine ComputeCheckFoldedALE(a)
   use MPI
   use typre
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: nelem, ielem,auxFoldedALE,auxfoldedALE2,ierr
   
   call a%ElementAlloc(e,a%Memor,'DefaultRule','ComputeCheckFoldedALE')
   call a%GetNelem(nelem)
   
   auxfoldedALE = 0
   do ielem = 1,nelem
      call a%ElementLoad(ielem,e)
      
      call e%elmdcg
      if (e%detjm <= 0.0_rp) then
         auxFoldedALE = 1
      endif
   enddo
   
   !communicate the result of the folded test
   call  MPI_AllREDUCE( auxfoldedALE, auxfoldedALE2, 1, MPI_INTEGER4, MPI_MAX,a%MPIcomm, ierr ) 
   
   if (auxfoldedALE2 == 1) then
      a%AreFoldedALE = .true.
   else
      a%AreFoldedALE = .false.
   endif
   
   call a%ElementDeAlloc(e,a%Memor,'DefaultRule','ComputeCheckFoldedALE')
   
end subroutine
