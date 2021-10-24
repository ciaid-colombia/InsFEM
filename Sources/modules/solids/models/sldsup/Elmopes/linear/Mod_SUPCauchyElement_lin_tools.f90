module Mod_SUPCauchyElement_lin_tools
   use typre
   use Mod_Element
   implicit none

contains

   subroutine sld_calculate_FdN_tensor_lin(e,nd,tn,inod,jnod,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the elemental material matrix 
       !              (F_jK*d(xi_ij)/dX_K,div(S_iq)_i)
       !
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            :: dxi,dyi,dzi,dxj,dyj,dzj
   
       B=0.0_rp

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           B(1,1) = dxi*dxj
           B(1,2) = 0.0_rp
           B(1,3) = dxi*dyj

           B(2,1) = 0.0_rp
           B(2,2) = dyi*dyj
           B(2,3) = dyi*dxj 

           B(3,1) = dyi*dxj 
           B(3,2) = dxi*dyj 
           B(3,3) = dyi*dyj&
           &     +  dxi*dxj

       else

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1,1) = dxi*dxj
           B(1,2) = 0.0_rp                         
           B(1,3) = 0.0_rp                         
           B(1,4) = 0.0_rp                         
           B(1,5) = dxi*dzj
           B(1,6) = dxi*dyj

           B(2,1) = 0.0_rp                         
           B(2,2) = dyi*dyj
           B(2,3) = 0.0_rp   
           B(2,4) = dyi*dzj
           B(2,5) = 0.0_rp   
           B(2,6) = dyi*dxj 

           B(3,1) = 0.0_rp                         
           B(3,2) = 0.0_rp                         
           B(3,3) = dzi*dzj
           B(3,4) = dzi*dyj
           B(3,5) = dzi*dxj 
           B(3,6) = 0.0_rp 

           B(4,1) = 0.0_rp                         
           B(4,2) = dzi*dyj
           B(4,3) = dyi*dzj                         
           B(4,4) = dzi*dzj &
                & + dyi*dyj
           B(4,5) = dyi*dxj  
           B(4,6) = dzi*dxj 

           B(5,1) = dzi*dxj
           B(5,2) = 0.0_rp 
           B(5,3) = dxi*dzj
           B(5,4) = dxi*dyj
           B(5,5) = dzi*dzj &
                & + dxi*dxj
           B(5,6) = dzi*dyj

           B(6,1) = dyi*dxj
           B(6,2) = dxi*dyj 
           B(6,3) = 0.0_rp
           B(6,4) = dxi*dzj
           B(6,5) = dyi*dzj
           B(6,6) = dyi*dyj &
                & + dxi*dxj

       endif

   end subroutine sld_calculate_FdN_tensor_lin

   subroutine sld_calculate_dqFinv_vector_lin(e,nd,tn,inod,jnod,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes
       !            -d/dX_K((Jp/lam)-1)*F^-1_Ki*q_h)*d(S_ij)/dx_j
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: inod,jnod,tn,nd
       real(rp)            , intent(inout):: B(tn)
       real(rp)            :: dxi,dyi,dzi,dxj,dyj,dzj
   
       B=0.0_rp

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           !-d(q)/dX_K*F^-1_Ki
           B(1) = dxi*dxj
           B(2) = dyi*dyj
           B(3) = dxi*dyj &
              & + dyi*dxj

       elseif (nd.eq.3) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1) = dxi*dxj
           B(2) = dyi*dyj
           B(3) = dzi*dzj
           B(4) = dyi*dzj &
              & + dzi*dyj
           B(5) = dxi*dzj &
              & + dzi*dxj
           B(6) = dxi*dyj &
              & + dyi*dxj

       endif

   end subroutine sld_calculate_dqFinv_vector_lin


end module Mod_SUPCauchyElement_lin_tools
