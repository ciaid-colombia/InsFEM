module Mod_php_AssemblyVectorToSystem 
   use typre
   use Mod_Element
   use Mod_PhysicalProblem
   implicit none

contains

   subroutine php_AssemblyVectorToSystem(a,Vector)
      class(PhysicalProblem) :: a
      real(rp) :: Vector(:,:)
      
      class(FiniteElement), pointer :: e => NULL()
      integer(ip) :: ipoin,npoin
      
      integer(ip), pointer :: aux_lnods(:) => NULL()
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','php')
      aux_Lnods => e%lnods
      
      call a%Memor%palloc(npoin,e%lnods,'lnods','php_assemblyVector')
      do ipoin = 1,npoin
         e%lnods(ipoin) = ipoin
      enddo
      e%pnode = npoin
      e%mnode = npoin
      call a%LinearSystem%AssemblyElrhs(e,Vector)
      
      call a%Memor%pdealloc(npoin,e%lnods,'lnods','php_assemblyVector')
      
      e%lnods => aux_Lnods
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','php')

   end subroutine


end module
