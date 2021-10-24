module Mod_UPCauchyElement_NH_tools
   use typre
   use Mod_Element
   implicit none

contains

   subroutine sld_calculateVTauTau_JS(e,nd,tn,inod,jnod,tautau,divstr,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !            (v_h,i*(1-tau_k^-1*tau_t), F^-1_Kl(du_l/dX_K)div(S_ij),i)
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: F(nd,nd),divstr(nd),tautau !(1-tau_k^-1*tau_t)
       real(rp)            , intent(inout):: B(nd,nd)
       real(rp)            :: shp,dx,dy,dz
   
       B    = 0.0_rp

       shp  = e%shape(inod,e%igaus)*(1.0_rp-tautau)

       if (nd.eq.2) then

           dx = e%cartd(1,jnod)
           dy = e%cartd(2,jnod)

           B(1,1) = (F(1,1)*dx + F(2,1)*dy)*divstr(1)
           B(1,2) = (F(1,2)*dx + F(2,2)*dy)*divstr(1)

           B(2,1) = (F(1,1)*dx + F(2,1)*dy)*divstr(2)
           B(2,2) = (F(1,2)*dx + F(2,2)*dy)*divstr(2)

       elseif (nd.eq.3) then

           dx = e%cartd(1,jnod)
           dy = e%cartd(2,jnod)
           dz = e%cartd(3,jnod)

           B(1,1) = (F(1,1)*dx + F(2,1)*dy + F(3,1)*dz)*divstr(1)
           B(1,2) = (F(1,2)*dx + F(2,2)*dy + F(3,2)*dz)*divstr(1)
           B(1,3) = (F(1,3)*dx + F(2,3)*dy + F(3,3)*dz)*divstr(1)

           B(2,1) = (F(1,1)*dx + F(2,1)*dy + F(3,1)*dz)*divstr(2)
           B(2,2) = (F(1,2)*dx + F(2,2)*dy + F(3,2)*dz)*divstr(2)
           B(2,3) = (F(1,3)*dx + F(2,3)*dy + F(3,3)*dz)*divstr(2)

           B(3,1) = (F(1,1)*dx + F(2,1)*dy + F(3,1)*dz)*divstr(3)
           B(3,2) = (F(1,2)*dx + F(2,2)*dy + F(3,2)*dz)*divstr(3)
           B(3,3) = (F(1,3)*dx + F(2,3)*dy + F(3,3)*dz)*divstr(3)

       endif

       B = shp*B

   end subroutine 

end module
