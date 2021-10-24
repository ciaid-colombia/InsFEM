module Mod_lmn_elmdir

contains

   subroutine lmn_elmdir(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      use Mod_php_elmdir
      use Mod_LowMach
      implicit none
      class(LowMachProblem)            :: a
      class(FiniteElement), intent(in) :: e
      real(rp), intent(inout)          :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)
      integer(ip)                      :: aux,aux1
  
 
      aux=a%ndofn
      aux1=a%ndofbcstart
      
      call php_elmdir(a,e,a%ndofn,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
   end subroutine      
        
end module
