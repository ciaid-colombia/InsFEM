module Mod_SUPCauchyElement_NH_tools
   use typre
   use Mod_Element
   implicit none

contains

   subroutine sld_calculate_2Fdu(e,nd,tn,nod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                     F_jK*du_i/dx_K+F_iK*du_j/dx_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,tn,nd
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(tn,nd)
       integer(ip)         :: i,j
       real(rp)            :: dx,dy,dz
   
       B=0.0_rp

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           B(1,1) = dx*F(1,1) + dy*F(1,2)
           B(1,2) = 0.0_rp

           B(2,1) = 0.0_rp
           B(2,2) = dx*F(2,1) + dy*F(2,2)

           B(3,1) = dx*F(2,1) + dy*F(2,2)
           B(3,2) = dx*F(1,1) + dy*F(1,2)

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           B(1,1) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)
           B(1,2) = 0.0_rp
           B(1,3) = 0.0_rp

           B(2,1) = 0.0_rp
           B(2,2) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3)
           B(2,3) = 0.0_rp

           B(3,1) = 0.0_rp
           B(3,2) = 0.0_rp
           B(3,3) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3)

           B(4,1) = 0.0_rp
           B(4,2) = F(3,1)*dx + F(3,2)*dy + F(3,3)*dz
           B(4,3) = F(2,1)*dx + F(2,2)*dy + F(2,3)*dz

           B(5,1) = F(3,1)*dx + F(3,2)*dy + F(3,3)*dz
           B(5,2) = 0.0_rp 
           B(5,3) = F(1,1)*dx + F(1,2)*dy + F(1,3)*dz

           B(6,1) = F(2,1)*dx + F(2,2)*dy + F(2,3)*dz
           B(6,2) = F(1,1)*dx + F(1,2)*dy + F(1,3)*dz
           B(6,3) = 0.0_rp

       endif

       B=2.0_rp*B

   end subroutine sld_calculate_2Fdu


   subroutine sld_calculate_traceFdu(e,nd,tn,nod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                   (1/nd*F_lK*du_l/dX_K*ii_ij)
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,tn,nd
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(tn,nd)
       real(rp)            :: dx,dy,dz

       B=0.0_rp

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           B(1,1) = dx*F(1,1) + dy*F(1,2)
           B(1,2) = dx*F(2,1) + dy*F(2,2)

           B(2,1) = dx*F(1,1) + dy*F(1,2)
           B(2,2) = dx*F(2,1) + dy*F(2,2)

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           B(1,1) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)
           B(1,2) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3)
           B(1,3) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3)
                                                   
           B(2,1) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)
           B(2,2) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3)
           B(2,3) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3)

           B(3,1) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)
           B(3,2) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3)
           B(3,3) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3)

       endif

   end subroutine sld_calculate_traceFdu

   subroutine sld_calculate_traceFduS(e,nd,tn,nod,F,gdev,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                    (J/2mu)F^-1_Kl*du_l/dX_K*S_ij
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,tn,nd
       real(rp)            , intent(in)   :: F(nd,nd),gdev(tn)
       real(rp)            , intent(inout):: B(tn,nd)
       real(rp)            :: dx,dy,dz

       B=0.0_rp

       !!This is the trace

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           B(1,1) = (dx*F(1,1) + dy*F(2,1))*gdev(1)
           B(1,2) = (dx*F(1,2) + dy*F(2,2))*gdev(1)

           B(2,1) = (dx*F(1,1) + dy*F(2,1))*gdev(2)
           B(2,2) = (dx*F(1,2) + dy*F(2,2))*gdev(2)

           B(3,1) = 2.0_rp*(dx*F(1,1) + dy*F(2,1))*gdev(3)
           B(3,2) = 2.0_rp*(dx*F(1,2) + dy*F(2,2))*gdev(3)

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           B(1,1) =        (dx*F(1,1) + dy*F(2,1) + dz*F(3,1))*gdev(1)
           B(1,2) =        (dx*F(1,2) + dy*F(2,2) + dz*F(3,2))*gdev(1)
           B(1,3) =        (dx*F(1,3) + dy*F(2,3) + dz*F(3,3))*gdev(1)

           B(2,1) =        (dx*F(1,1) + dy*F(2,1) + dz*F(3,1))*gdev(2)
           B(2,2) =        (dx*F(1,2) + dy*F(2,2) + dz*F(3,2))*gdev(2)
           B(2,3) =        (dx*F(1,3) + dy*F(2,3) + dz*F(3,3))*gdev(2)

           B(3,1) =        (dx*F(1,1) + dy*F(2,1) + dz*F(3,1))*gdev(3)
           B(3,2) =        (dx*F(1,2) + dy*F(2,2) + dz*F(3,2))*gdev(3)
           B(3,3) =        (dx*F(1,3) + dy*F(2,3) + dz*F(3,3))*gdev(3)

           B(4,1) = 2.0_rp*(dx*F(1,1) + dy*F(2,1) + dz*F(3,1))*gdev(4)
           B(4,2) = 2.0_rp*(dx*F(1,2) + dy*F(2,2) + dz*F(3,2))*gdev(4)
           B(4,3) = 2.0_rp*(dx*F(1,3) + dy*F(2,3) + dz*F(3,3))*gdev(4)
                           
           B(5,1) = 2.0_rp*(dx*F(1,1) + dy*F(2,1) + dz*F(3,1))*gdev(5)
           B(5,2) = 2.0_rp*(dx*F(1,2) + dy*F(2,2) + dz*F(3,2))*gdev(5)
           B(5,3) = 2.0_rp*(dx*F(1,3) + dy*F(2,3) + dz*F(3,3))*gdev(5)
                           
           B(6,1) = 2.0_rp*(dx*F(1,1) + dy*F(2,1) + dz*F(3,1))*gdev(6)
           B(6,2) = 2.0_rp*(dx*F(1,2) + dy*F(2,2) + dz*F(3,2))*gdev(6)
           B(6,3) = 2.0_rp*(dx*F(1,3) + dy*F(2,3) + dz*F(3,3))*gdev(6)

       endif

   end subroutine sld_calculate_traceFduS

   subroutine sld_calculate_traceFduS_tensor(e,nd,tn,inod,jnod,F,gdev,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !       1/2(dv_j/dx_i+dv_i/dx_j,(J/2mu)F^-1_Kl*du_l/dX_K*S_ij)
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,tn,nd
       real(rp)            , intent(in)   :: F(nd,nd),gdev(tn)
       real(rp)            , intent(inout):: B(nd,nd)
       real(rp)            :: dxi,dyi,dzi
       real(rp)            :: dxj,dyj,dzj

       B=0.0_rp

       if (nd.eq.2) then

           dxi=e%cartd(1,inod)
           dyi=e%cartd(2,inod)

           dxj=e%cartd(1,jnod)
           dyj=e%cartd(2,jnod)

           !1/2(dv_1/dx_1+dv_1/dx_1,(J/2mu)F^-1_Kl*du_l/dX_K*S_11)
           B(1,1) = 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(2,1))*gdev(1)
           B(1,2) = 2.0_rp*dxi*(dxj*F(1,2) + dyj*F(2,2))*gdev(1)

           !1/2(dv_2/dx_2+dv_2/dx_2,(J/2mu)F^-1_Kl*du_l/dX_K*S_22)
           B(2,1) = 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(2,1))*gdev(2)
           B(2,2) = 2.0_rp*dyi*(dxj*F(1,2) + dyj*F(2,2))*gdev(2)

           !1/2(dv_1/dx_2+dv_2/dx_1,(J/2mu)F^-1_Kl*du_l/dX_K*S_12)
           !1/2(dv_2/dx_1+dv_1/dx_2,(J/2mu)F^-1_Kl*du_l/dX_K*S_21)
           B(1,1) = B(1,1) + 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(2,1))*gdev(3)
           B(1,2) = B(1,2) + 2.0_rp*dyi*(dxj*F(1,2) + dyj*F(2,2))*gdev(3)

           B(2,1) = B(2,1) + 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(2,1))*gdev(3)
           B(2,2) = B(2,2) + 2.0_rp*dxi*(dxj*F(1,2) + dyj*F(2,2))*gdev(3)

       elseif (nd.eq.3) then

           dxi=e%cartd(1,inod)
           dyi=e%cartd(2,inod)
           dzi=e%cartd(3,inod)

           dxj=e%cartd(1,jnod)
           dyj=e%cartd(2,jnod)
           dzj=e%cartd(3,jnod)

           !1/2(dv_1/dx_1+dv_1/dx_1,(J/2mu)F^-1_Kl*du_l/dX_K*S_11)
           B(1,1) = 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(1)
           B(1,2) = 2.0_rp*dxi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(1)
           B(1,3) = 2.0_rp*dxi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(1)

           !1/2(dv_2/dx_2+dv_2/dx_2,(J/2mu)F^-1_Kl*du_l/dX_K*S_22)
           B(2,1) = 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(2)
           B(2,2) = 2.0_rp*dyi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(2)
           B(2,3) = 2.0_rp*dyi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(2)

           !1/2(dv_3/dx_3+dv_3/dx_3,(J/2mu)F^-1_Kl*du_l/dX_K*S_33)
           B(3,1) = 2.0_rp*dzi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(3)
           B(3,2) = 2.0_rp*dzi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(3)
           B(3,3) = 2.0_rp*dzi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(3)

           !1/2(dv_2/dx_3+dv_3/dx_2,(J/2mu)F^-1_Kl*du_l/dX_K*S_23)
           !1/2(dv_3/dx_2+dv_2/dx_3,(J/2mu)F^-1_Kl*du_l/dX_K*S_32)
           B(2,1) = B(2,1) + 2.0_rp*dzi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(4)
           B(2,2) = B(2,2) + 2.0_rp*dzi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(4)
           B(2,3) = B(2,3) + 2.0_rp*dzi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(4)

           B(3,1) = B(3,1) + 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(4)
           B(3,2) = B(3,2) + 2.0_rp*dyi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(4)
           B(3,3) = B(3,3) + 2.0_rp*dyi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(4)

           !1/2(dv_1/dx_3+dv_3/dx_1,(J/2mu)F^-1_Kl*du_l/dX_K*S_13)
           !1/2(dv_3/dx_1+dv_1/dx_3,(J/2mu)F^-1_Kl*du_l/dX_K*S_31)
           B(1,1) = B(1,1) + 2.0_rp*dzi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(5)
           B(1,2) = B(1,2) + 2.0_rp*dzi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(5)
           B(1,3) = B(1,3) + 2.0_rp*dzi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(5)

           B(3,1) = B(3,1) + 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(5)
           B(3,2) = B(3,2) + 2.0_rp*dxi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(5)
           B(3,3) = B(3,3) + 2.0_rp*dxi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(5)

           !1/2(dv_1/dx_2+dv_2/dx_1,(J/2mu)F^-1_Kl*du_l/dX_K*S_12)
           !1/2(dv_2/dx_1+dv_1/dx_2,(J/2mu)F^-1_Kl*du_l/dX_K*S_21)
           B(1,1) = B(1,1) + 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(6)
           B(1,2) = B(1,2) + 2.0_rp*dyi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(6)
           B(1,3) = B(1,3) + 2.0_rp*dyi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(6)

           B(2,1) = B(2,1) + 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(2,1) + dzj*F(3,1))*gdev(6)
           B(2,2) = B(2,2) + 2.0_rp*dxi*(dxj*F(1,2) + dyj*F(2,2) + dzj*F(3,2))*gdev(6)
           B(2,3) = B(2,3) + 2.0_rp*dxi*(dxj*F(1,3) + dyj*F(2,3) + dzj*F(3,3))*gdev(6)

       endif

   end subroutine sld_calculate_traceFduS_tensor

   subroutine sld_calculate_tracePressFinv(e,nd,nod,C,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, voigt
       !                         F^-1_Kl*du_l/dx_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: nod,nd
       real(rp)            , intent(in)   :: C(nd,nd)
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: dx,dy,dz
   
       B=0.0_rp

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           B(1) = dx*C(1,1) + dy*C(2,1)
           B(2) = dx*C(1,2) + dy*C(2,2)

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           B(1) = dx*C(1,1) + dy*C(2,1) + dz*C(3,1)
           B(2) = dx*C(1,2) + dy*C(2,2) + dz*C(3,2)
           B(3) = dx*C(1,3) + dy*C(2,3) + dz*C(3,3)

       endif

   end subroutine sld_calculate_tracePressFinv

   subroutine sld_calculate_tracePressFmat(e,nd,nod,C,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, voigt
       !                       F_lK*du_l/dx_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: nod,nd
       real(rp)            , intent(in)   :: C(nd,nd)
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: dx,dy,dz
   
       B=0.0_rp

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           B(1) = dx*C(1,1) + dy*C(1,2)
           B(2) = dx*C(2,1) + dy*C(2,2)

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           B(1) = dx*C(1,1) + dy*C(1,2) + dz*C(1,3)
           B(2) = dx*C(2,1) + dy*C(2,2) + dz*C(2,3)
           B(3) = dx*C(3,1) + dy*C(3,2) + dz*C(3,3)

       endif

   end subroutine sld_calculate_tracePressFmat

   !------------------SUP NH SGS FUNCTIONS-------------------------

   subroutine sld_calculate_2Fdu_tensor(e,nd,inod,jnod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !      (dv_j/dx_i+dv_i/dx_j,F_iK*du_j/dX_K + F_jK*du_i/dX_K)
       !-----------------------------------------------------------------------
       implicit none

       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(nd,nd)
       real(rp)            :: dxi,dyi,dzi
       real(rp)            :: dxj,dyj,dzj

       B=0.0_rp

       if (nd.eq.2) then

           dxi=e%cartd(1,inod)
           dyi=e%cartd(2,inod)

           dxj=e%cartd(1,jnod)
           dyj=e%cartd(2,jnod)

           !(dv_i/dx_j,F_iK*du_j/dX_K)
           B(1,1) =          dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(1,2) =          dyi*(dxj*F(1,1) + dyj*F(1,2))

           B(2,1) =          dxi*(dxj*F(2,1) + dyj*F(2,2))
           B(2,2) =          dyi*(dxj*F(2,1) + dyj*F(2,2))

           !(dv_i/dx_j,F_jK*du_i/dX_K)
           B(1,1) = B(1,1) + dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(1,1) = B(1,1) + dyi*(dxj*F(2,1) + dyj*F(2,2))

           B(2,2) = B(2,2) + dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(2,2) = B(2,2) + dyi*(dxj*F(2,1) + dyj*F(2,2))

           !(dv_j/dx_i,F_iK*du_j/dX_K)
           B(1,1) = B(1,1) + dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(1,1) = B(1,1) + dyi*(dxj*F(2,1) + dyj*F(2,2))

           B(2,2) = B(2,2) + dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(2,2) = B(2,2) + dyi*(dxj*F(2,1) + dyj*F(2,2))

           !(dv_j/dx_i,F_jK*du_i/dX_K)
           B(1,1) = B(1,1) + dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(1,2) = B(1,2) + dyi*(dxj*F(1,1) + dyj*F(1,2))

           B(2,1) = B(2,1) + dxi*(dxj*F(2,1) + dyj*F(2,2))
           B(2,2) = B(2,2) + dyi*(dxj*F(2,1) + dyj*F(2,2))

       elseif (nd.eq.3) then

           dxi=e%cartd(1,inod)
           dyi=e%cartd(2,inod)
           dzi=e%cartd(3,inod)

           dxj=e%cartd(1,jnod)
           dyj=e%cartd(2,jnod)
           dzj=e%cartd(3,jnod)

           !(dv_i/dx_j,F_iK*du_j/dX_K)
           B(1,1) =  dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,2) =  dyi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,3) =  dzi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))

           B(2,1) =  dxi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,2) =  dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,3) =  dzi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))

           B(3,1) =  dxi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))
           B(3,2) =  dyi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))
           B(3,3) =  dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           !(dv_i/dx_j,F_jK*du_i/dX_K)
           B(1,1) =B(1,1)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,1) =B(1,1)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(1,1) =B(1,1)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           B(2,2) =B(2,2)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(2,2) =B(2,2)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,2) =B(2,2)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           B(3,3) =B(3,3)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(3,3) =B(3,3)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(3,3) =B(3,3)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           !(dv_j/dx_i,F_iK*du_j/dX_K)
           B(1,1) = B(1,1)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,1) = B(1,1)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(1,1) = B(1,1)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           B(2,2) = B(2,2)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(2,2) = B(2,2)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,2) = B(2,2)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           B(3,3) = B(3,3)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(3,3) = B(3,3)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(3,3) = B(3,3)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           !(dv_j/dx_i,F_jK*du_i/dX_K)
           B(1,1) = B(1,1)+ dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,2) = B(1,2)+ dyi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,3) = B(1,3)+ dzi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))

           B(2,1) = B(2,1)+ dxi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,2) = B(2,2)+ dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,3) = B(2,3)+ dzi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))

           B(3,1) = B(3,1)+ dxi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))
           B(3,2) = B(3,2)+ dyi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))
           B(3,3) = B(3,3)+ dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

       endif

   end subroutine sld_calculate_2Fdu_tensor

   subroutine sld_calculate_traceFdu_tensor(e,nd,inod,jnod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !    1/2(dv_j/dx_i+dv_i/dx_j,F_lK*du_l/dX_K*ii_ij)
       !-----------------------------------------------------------------------
       implicit none

       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(nd,nd)
       real(rp)            :: dxi,dyi,dzi
       real(rp)            :: dxj,dyj,dzj

       B=0.0_rp

       if (nd.eq.2) then

           dxi=e%cartd(1,inod)
           dyi=e%cartd(2,inod)

           dxj=e%cartd(1,jnod)
           dyj=e%cartd(2,jnod)

           B(1,1) = 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(1,2))
           B(1,2) = 2.0_rp*dxi*(dxj*F(2,1) + dyj*F(2,2))

           B(2,1) = 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(1,2))
           B(2,2) = 2.0_rp*dyi*(dxj*F(2,1) + dyj*F(2,2))

       elseif (nd.eq.3) then

           dxi=e%cartd(1,inod)
           dyi=e%cartd(2,inod)
           dzi=e%cartd(3,inod)

           dxj=e%cartd(1,jnod)
           dyj=e%cartd(2,jnod)
           dzj=e%cartd(3,jnod)

           B(1,1) = 2.0_rp*dxi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(1,2) = 2.0_rp*dxi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(1,3) = 2.0_rp*dxi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           B(2,1) = 2.0_rp*dyi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(2,2) = 2.0_rp*dyi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(2,3) = 2.0_rp*dyi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

           B(3,1) = 2.0_rp*dzi*(dxj*F(1,1) + dyj*F(1,2) + dzj*F(1,3))
           B(3,2) = 2.0_rp*dzi*(dxj*F(2,1) + dyj*F(2,2) + dzj*F(2,3))
           B(3,3) = 2.0_rp*dzi*(dxj*F(3,1) + dyj*F(3,2) + dzj*F(3,3))

       endif

   end subroutine sld_calculate_traceFdu_tensor


   subroutine sld_calculate_trace_dNS_F(e,nd,tn,nod,F,grsig,gdev,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, voigt
       !                  F^-1_Ki*d(xi_lm*S_lm)/dX_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement),intent(in)   :: e
       integer(ip)         ,intent(in)   :: nod,tn,nd
       real(rp)            ,intent(in)   :: F(nd,nd),grsig(tn,nd),gdev(tn)
       real(rp)            ,intent(inout):: B(tn,nd)
       real(rp)                          :: shp,dx,dy,dz

       B=0.0_rp
       shp = e%shape(nod,e%igaus)

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           B(1,1)= dx*gdev(1)*F(1,1) + shp*grsig(1,1)*F(1,1) &
               & + dy*gdev(1)*F(2,1) + shp*grsig(1,2)*F(2,1)

           B(1,2)= dx*gdev(1)*F(1,2) + shp*grsig(1,1)*F(1,2) &
               & + dy*gdev(1)*F(2,2) + shp*grsig(1,2)*F(2,2)

           B(2,1)= dx*gdev(2)*F(1,1) + shp*grsig(2,1)*F(1,1) &
               & + dy*gdev(2)*F(2,1) + shp*grsig(2,2)*F(2,1)

           B(2,2)= dx*gdev(2)*F(1,2) + shp*grsig(2,1)*F(1,2) &
               & + dy*gdev(2)*F(2,2) + shp*grsig(2,2)*F(2,2)

           B(3,1)= 2.0_rp*(dx*gdev(3)*F(1,1) + shp*grsig(3,1)*F(1,1)) &
               & + 2.0_rp*(dy*gdev(3)*F(2,1) + shp*grsig(3,2)*F(2,1))
                                                                            
           B(3,2)= 2.0_rp*(dx*gdev(3)*F(1,2) + shp*grsig(3,1)*F(1,2)) &
               & + 2.0_rp*(dy*gdev(3)*F(2,2) + shp*grsig(3,2)*F(2,2))

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           B(1,1)= dx*gdev(1)*F(1,1) + shp*grsig(1,1)*F(1,1) &
               & + dy*gdev(1)*F(2,1) + shp*grsig(1,2)*F(2,1) &
               & + dz*gdev(1)*F(3,1) + shp*grsig(1,3)*F(3,1)

           B(1,2)= dx*gdev(1)*F(1,2) + shp*grsig(1,1)*F(1,2) &
               & + dy*gdev(1)*F(2,2) + shp*grsig(1,2)*F(2,2) &
               & + dz*gdev(1)*F(3,2) + shp*grsig(1,3)*F(3,2)

           B(1,3)= dx*gdev(1)*F(1,3) + shp*grsig(1,1)*F(1,3) &
               & + dy*gdev(1)*F(2,3) + shp*grsig(1,2)*F(2,3) &
               & + dz*gdev(1)*F(3,3) + shp*grsig(1,3)*F(3,3)

           B(2,1)= dx*gdev(2)*F(1,1) + shp*grsig(2,1)*F(1,1) &
               & + dy*gdev(2)*F(2,1) + shp*grsig(2,2)*F(2,1) &
               & + dz*gdev(2)*F(3,1) + shp*grsig(2,3)*F(3,1)

           B(2,2)= dx*gdev(2)*F(1,2) + shp*grsig(2,1)*F(1,2) &
               & + dy*gdev(2)*F(2,2) + shp*grsig(2,2)*F(2,2) &
               & + dz*gdev(2)*F(3,2) + shp*grsig(2,3)*F(3,2)

           B(2,3)= dx*gdev(2)*F(1,3) + shp*grsig(2,1)*F(1,3) &
               & + dy*gdev(2)*F(2,3) + shp*grsig(2,2)*F(2,3) &
               & + dz*gdev(2)*F(3,3) + shp*grsig(2,3)*F(3,3)

           B(3,1)= dx*gdev(3)*F(1,1) + shp*grsig(3,1)*F(1,1) &
               & + dy*gdev(3)*F(2,1) + shp*grsig(3,2)*F(2,1) &
               & + dz*gdev(3)*F(3,1) + shp*grsig(3,3)*F(3,1)

           B(3,2)= dx*gdev(3)*F(1,2) + shp*grsig(3,1)*F(1,2) &
               & + dy*gdev(3)*F(2,2) + shp*grsig(3,2)*F(2,2) &
               & + dz*gdev(3)*F(3,2) + shp*grsig(3,3)*F(3,2)

           B(3,3)= dx*gdev(3)*F(1,3) + shp*grsig(3,1)*F(1,3) &
               & + dy*gdev(3)*F(2,3) + shp*grsig(3,2)*F(2,3) &
               & + dz*gdev(3)*F(3,3) + shp*grsig(3,3)*F(3,3)

           B(4,1)= 2.0_rp*(dx*gdev(4)*F(1,1) + shp*grsig(4,1)*F(1,1)) &
               & + 2.0_rp*(dy*gdev(4)*F(2,1) + shp*grsig(4,2)*F(2,1)) &
               & + 2.0_rp*(dz*gdev(4)*F(3,1) + shp*grsig(4,3)*F(3,1))

           B(4,2)= 2.0_rp*(dx*gdev(4)*F(1,2) + shp*grsig(4,1)*F(1,2)) &
               & + 2.0_rp*(dy*gdev(4)*F(2,2) + shp*grsig(4,2)*F(2,2)) &
               & + 2.0_rp*(dz*gdev(4)*F(3,2) + shp*grsig(4,3)*F(3,2))

           B(4,3)= 2.0_rp*(dx*gdev(4)*F(1,3) + shp*grsig(4,1)*F(1,3)) &
               & + 2.0_rp*(dy*gdev(4)*F(2,3) + shp*grsig(4,2)*F(2,3)) &
               & + 2.0_rp*(dz*gdev(4)*F(3,3) + shp*grsig(4,3)*F(3,3))
                                                                          
           B(5,1)= 2.0_rp*(dx*gdev(5)*F(1,1) + shp*grsig(5,1)*F(1,1)) &
               & + 2.0_rp*(dy*gdev(5)*F(2,1) + shp*grsig(5,2)*F(2,1)) &
               & + 2.0_rp*(dz*gdev(5)*F(3,1) + shp*grsig(5,3)*F(3,1))

           B(5,2)= 2.0_rp*(dx*gdev(5)*F(1,2) + shp*grsig(5,1)*F(1,2)) &
               & + 2.0_rp*(dy*gdev(5)*F(2,2) + shp*grsig(5,2)*F(2,2)) &
               & + 2.0_rp*(dz*gdev(5)*F(3,2) + shp*grsig(5,3)*F(3,2))

           B(5,3)= 2.0_rp*(dx*gdev(5)*F(1,3) + shp*grsig(5,1)*F(1,3)) &
               & + 2.0_rp*(dy*gdev(5)*F(2,3) + shp*grsig(5,2)*F(2,3)) &
               & + 2.0_rp*(dz*gdev(5)*F(3,3) + shp*grsig(5,3)*F(3,3))

           B(6,1)= 2.0_rp*(dx*gdev(6)*F(1,1) + shp*grsig(6,1)*F(1,1)) &
               & + 2.0_rp*(dy*gdev(6)*F(2,1) + shp*grsig(6,2)*F(2,1)) &
               & + 2.0_rp*(dz*gdev(6)*F(3,1) + shp*grsig(6,3)*F(3,1))

           B(6,2)= 2.0_rp*(dx*gdev(6)*F(1,2) + shp*grsig(6,1)*F(1,2)) &
               & + 2.0_rp*(dy*gdev(6)*F(2,2) + shp*grsig(6,2)*F(2,2)) &
               & + 2.0_rp*(dz*gdev(6)*F(3,2) + shp*grsig(6,3)*F(3,2))

           B(6,3)= 2.0_rp*(dx*gdev(6)*F(1,3) + shp*grsig(6,1)*F(1,3)) &
               & + 2.0_rp*(dy*gdev(6)*F(2,3) + shp*grsig(6,2)*F(2,3)) &
               & + 2.0_rp*(dz*gdev(6)*F(3,3) + shp*grsig(6,3)*F(3,3))

endif

   end subroutine sld_calculate_trace_dNS_F

   subroutine sld_calculate_FdN(e,nd,tn,nod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !             (F_jK*d(xi_ij)/dX_K - 1/nd*F_iK*d(xi_ll)/dX_K)
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,nd,tn
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(tn,nd)
       real(rp)            :: aux,dx,dy,dz
   
       B   = 0.0_rp
       aux = 1.0_rp/nd

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           !F_iK*d(xi_ij)/dX_K 
           B(1,1) = dx*F(1,1) + dy*F(1,2)  
           B(1,2) = 0.0_rp

           B(2,1) = 0.0_rp
           B(2,2) = dx*F(2,1) + dy*F(2,2)

           B(3,1) = dx*F(2,1) + dy*F(2,2)  
           B(3,2) = dx*F(1,1) + dy*F(1,2)  

           !- 1/nd*F_lK*d(xi_ll)/dX_K*ii_ij
           B(1,1) = B(1,1)-aux*(dx*F(1,1) + dy*F(1,2))  
           B(1,2) = B(1,2)-aux*(dx*F(2,1) + dy*F(2,2))  
           B(2,1) = B(2,1)-aux*(dx*F(1,1) + dy*F(1,2))  
           B(2,2) = B(2,2)-aux*(dx*F(2,1) + dy*F(2,2))  

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           !F_iK*d(xi_ij)/dX_K 
           B(1,1) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)
           B(1,2) = 0.0_rp
           B(1,3) = 0.0_rp

           B(2,1) = 0.0_rp
           B(2,2) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3)
           B(2,3) = 0.0_rp

           B(3,1) = 0.0_rp
           B(3,2) = 0.0_rp
           B(3,3) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3)

           B(4,1) = 0.0_rp
           B(4,2) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3) 
           B(4,3) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3)

           B(5,1) = dx*F(3,1) + dy*F(3,2) + dz*F(3,3) 
           B(5,2) = 0.0_rp
           B(5,3) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)

           B(6,1) = dx*F(2,1) + dy*F(2,2) + dz*F(2,3) 
           B(6,2) = dx*F(1,1) + dy*F(1,2) + dz*F(1,3)
           B(6,3) = 0.0_rp

           B(1,1) = B(1,1)-aux*(dx*F(1,1) + dy*F(1,2) + dz*F(1,3))  
           B(1,2) = B(1,2)-aux*(dx*F(2,1) + dy*F(2,2) + dz*F(2,3))  
           B(1,3) = B(1,3)-aux*(dx*F(3,1) + dy*F(3,2) + dz*F(3,3))  

           B(2,1) = B(2,1)-aux*(dx*F(1,1) + dy*F(1,2) + dz*F(1,3))  
           B(2,2) = B(2,2)-aux*(dx*F(2,1) + dy*F(2,2) + dz*F(2,3))  
           B(2,3) = B(2,3)-aux*(dx*F(3,1) + dy*F(3,2) + dz*F(3,3))  

           B(3,1) = B(3,1)-aux*(dx*F(1,1) + dy*F(1,2) + dz*F(1,3))  
           B(3,2) = B(3,2)-aux*(dx*F(2,1) + dy*F(2,2) + dz*F(2,3))  
           B(3,3) = B(3,3)-aux*(dx*F(3,1) + dy*F(3,2) + dz*F(3,3))  

       endif
   
   end subroutine sld_calculate_FdN

   subroutine sld_calculate_trace_dNS_F_tensor(e,nd,tn,inod,jnod,F,grsig,gdev,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !            (F^-1_Ki*d(xi_lm)/dX_K*S_ml*d(S_iq)/dx_q
       !            +F^-1_Ki*xi_lm*d(S_ml)/dX_K*S_ml*d(S_iq)/dx_q)
       !
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: F(nd,nd),grsig(tn,nd),gdev(tn)
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            :: shp,dxi,dyi,dzi
       real(rp)            :: dxj,dyj,dzj,aux1,aux2,aux3,aux4,aux5,aux6

       B=0.0_rp
       shp = e%shape(inod,e%igaus)

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           aux1 = (dxi*gdev(1)*F(1,1) + shp*grsig(1,1)*F(1,1))
           aux2 = (dyi*gdev(1)*F(2,1) + shp*grsig(1,2)*F(2,1))
           B(1,1)=   (aux1 + aux2)*dxj 
           B(1,2)=   (aux1 + aux2)*dyj 
           B(1,3)=   (aux1 + aux2)*dyj + (aux1 + aux2)*dxj 
           
           aux3 = (dxi*gdev(2)*F(1,1) + shp*grsig(2,1)*F(1,1))
           aux4 = (dyi*gdev(2)*F(2,1) + shp*grsig(2,2)*F(2,1))
           B(2,1)=   (aux3 + aux4)*dxj 
           B(2,2)=   (aux3 + aux4)*dyj 
           B(2,3)=   (aux3 + aux4)*dyj + (aux3 + aux4)*dxj 

           aux5 = (dxi*gdev(3)*F(1,1) + shp*grsig(3,1)*F(1,1))
           aux6 = (dyi*gdev(3)*F(2,1) + shp*grsig(3,2)*F(2,1))

           B(3,1)=   (aux5 + aux6)*dxj 
           B(3,2)=   (aux5 + aux6)*dyj 
           B(3,3)=   (aux5 + aux6)*dyj + (aux5 + aux6)*dxj 

       elseif (nd.eq.3) then

           !TODO
           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           aux1 = (dxi*gdev(1)*F(1,1) + shp*grsig(1,1)*F(1,1))
           aux2 = (dyi*gdev(1)*F(2,1) + shp*grsig(1,2)*F(2,1))
           aux3 = (dzi*gdev(1)*F(3,1) + shp*grsig(1,3)*F(3,1))
           B(1,1)=   (aux1 + aux2 + aux3)*dxj 
           B(1,2)=   (aux1 + aux2 + aux3)*dyj 
           B(1,3)=   (aux1 + aux2 + aux3)*dzj 
           B(1,3)=   (aux1 + aux2 + aux3)*dyj + (aux1 + aux2)*dxj 
           
           aux3 = (dxi*gdev(2)*F(1,1) + shp*grsig(2,1)*F(1,1))
           aux4 = (dyi*gdev(2)*F(2,1) + shp*grsig(2,2)*F(2,1))
           B(2,1)=   (aux3 + aux4)*dxj 
           B(2,2)=   (aux3 + aux4)*dyj 
           B(2,3)=   (aux3 + aux4)*dyj + (aux3 + aux4)*dxj 

           aux5 = (dxi*gdev(3)*F(1,1) + shp*grsig(3,1)*F(1,1))
           aux6 = (dyi*gdev(3)*F(2,1) + shp*grsig(3,2)*F(2,1))

           B(3,1)=   (aux5 + aux6)*dxj 
           B(3,2)=   (aux5 + aux6)*dyj 
           B(3,3)=   (aux5 + aux6)*dyj + (aux5 + aux6)*dxj 

endif

   end subroutine sld_calculate_trace_dNS_F_tensor

   subroutine sld_calculate_FdN_tensor_trace(e,nd,tn,inod,jnod,F,B)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !       ((F_iK*d(xi_ll)/dX_K),div(S_iq)_i)
      !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            ::dxi,dyi,dzi,dxj,dyj,dzj
   
       B=0.0_rp

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           B(1,1) = (dxi*F(1,1)+dyi*F(1,2))*dxj
           B(1,2) = (dxi*F(2,1)+dyi*F(2,2))*dyj
           B(1,3) = (dxi*F(1,1)+dyi*F(1,2))*dyj &
                & + (dxi*F(2,1)+dyi*F(2,2))*dxj

           B(2,1) = (dxi*F(1,1)+dyi*F(1,2))*dxj 
           B(2,2) = (dxi*F(2,1)+dyi*F(2,2))*dyj 
           B(2,3) = (dxi*F(1,1)+dyi*F(1,2))*dyj & 
                & + (dxi*F(2,1)+dyi*F(2,2))*dxj

       else

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1,1) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dxj
           B(1,2) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dyj
           B(1,3) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dzj
           B(1,4) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dzj &
                & + (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dyj
           B(1,5) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dzj &
                & + (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dxj
           B(1,6) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dyj &
                & + (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dxj

           B(2,1) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dxj
           B(2,2) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dyj
           B(2,3) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dzj
           B(2,4) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dzj &
                & + (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dyj
           B(2,5) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dzj &
                & + (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dxj
           B(2,6) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dyj &
                & + (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dxj

           B(3,1) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dxj
           B(3,2) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dyj
           B(3,3) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dzj
           B(3,4) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dzj &
                & + (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dyj
           B(3,5) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dzj &
                & + (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dxj
           B(3,6) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dyj &
                & + (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dxj

       endif

   end subroutine sld_calculate_FdN_tensor_trace

   subroutine sld_calculate_FdN_tensor(e,nd,tn,inod,jnod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the elemental material matrix 
       !              (F_jK*d(xi_ij)/dX_K,div(S_iq)_i)
       !
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            :: dxi,dyi,dzi,dxj,dyj,dzj
   
       B=0.0_rp

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           B(1,1) = (dxi*F(1,1)+dyi*F(1,2))*dxj
           B(1,2) = 0.0_rp
           B(1,3) = (dxi*F(1,1)+dyi*F(1,2))*dyj

           B(2,1) = 0.0_rp
           B(2,2) = (dxi*F(2,1)+dyi*F(2,2))*dyj
           B(2,3) = (dxi*F(2,1)+dyi*F(2,2))*dxj 

           B(3,1) = (dxi*F(2,1)+dyi*F(2,2))*dxj 
           B(3,2) = (dxi*F(1,1)+dyi*F(1,2))*dyj 
           B(3,3) = (dxi*F(2,1)+dyi*F(2,2))*dyj&
           &     +  (dxi*F(1,1)+dyi*F(1,2))*dxj

       else

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1,1) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dxj
           B(1,2) = 0.0_rp                         
           B(1,3) = 0.0_rp                         
           B(1,4) = 0.0_rp                         
           B(1,5) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dzj
           B(1,6) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dyj

           B(2,1) = 0.0_rp                         
           B(2,2) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dyj
           B(2,3) = 0.0_rp                         
           B(2,4) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dzj
           B(2,5) = 0.0_rp                         
           B(2,6) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dxj 

           B(3,1) = 0.0_rp                         
           B(3,2) = 0.0_rp                         
           B(3,3) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dzj
           B(3,4) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dyj
           B(3,5) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dxj 
           B(3,6) = 0.0_rp 

           B(4,1) = 0.0_rp                         
           B(4,2) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dyj
           B(4,3) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dzj                         
           B(4,4) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dzj &
                & + (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dyj
           B(4,5) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dxj  
           B(4,6) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dxj 

           B(5,1) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dxj
           B(5,2) = 0.0_rp 
           B(5,3) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dzj
           B(5,4) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dyj
           B(5,5) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dzj &
                & + (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dxj
           B(5,6) = (dxi*F(3,1)+dyi*F(3,2)+dzi*F(3,3))*dyj

           B(6,1) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dxj
           B(6,2) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dyj 
           B(6,3) = 0.0_rp
           B(6,4) = (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dzj
           B(6,5) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dzj
           B(6,6) = (dxi*F(2,1)+dyi*F(2,2)+dzi*F(2,3))*dyj &
                & + (dxi*F(1,1)+dyi*F(1,2)+dzi*F(1,3))*dxj

       endif

   end subroutine sld_calculate_FdN_tensor

   subroutine sld_calculate_dqF_vector(e,nd,tn,inod,jnod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                  F_iK*d(q_h)/dX_K*d(S_ij)/dx_j
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: inod,jnod,tn,nd
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(tn)
       integer(ip)         :: i,j
       real(rp)            :: dxi,dyi,dzi,dxj,dyj,dzj
   
       B=0.0_rp

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           B(1) = (F(1,1)*dxi+F(1,2)*dyi)*dxj
           B(2) = (F(2,1)*dxi+F(2,2)*dyi)*dyj
           B(3) = (F(1,1)*dxi+F(1,2)*dyi)*dyj &
              & + (F(2,1)*dxi+F(2,2)*dyi)*dxj

       elseif (nd.eq.3) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1) = (F(1,1)*dxi+F(1,2)*dyi+F(1,3)*dzi)*dxj
           B(2) = (F(2,1)*dxi+F(2,2)*dyi+F(2,3)*dzi)*dyj
           B(3) = (F(3,1)*dxi+F(3,2)*dyi+F(3,3)*dzi)*dzj
           B(4) = (F(2,1)*dxi+F(2,2)*dyi+F(2,3)*dzi)*dzj &
              & + (F(3,1)*dxi+F(3,2)*dyi+F(3,3)*dzi)*dyj
           B(5) = (F(1,1)*dxi+F(1,2)*dyi+F(1,3)*dzi)*dzj &
              & + (F(3,1)*dxi+F(3,2)*dyi+F(3,3)*dzi)*dxj
           B(6) = (F(1,1)*dxi+F(1,2)*dyi+F(1,3)*dzi)*dyj &
              & + (F(2,1)*dxi+F(2,2)*dyi+F(2,3)*dzi)*dxj

       endif

   end subroutine sld_calculate_dqF_vector

   subroutine sld_calculate_dqFinv_vector(e,nd,tn,inod,jnod,J,lam,gp,p,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes
       !            -d/dX_K((Jp/lam)-1)*F^-1_Ki*q_h)*d(S_ij)/dx_j
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: inod,jnod,tn,nd
       real(rp)            , intent(in)   :: F(nd,nd),J,lam,p,gp(nd)
       real(rp)            , intent(inout):: B(tn)
       integer(ip)         :: i,k
       real(rp)            :: dxi,dyi,dzi,dxj,dyj,dzj,shp
       real(rp)            :: aux1,aux2,aux3,Jlam
   
       B=0.0_rp

       Jlam = J/lam
       shp = e%shape(inod,e%igaus)

       if (nd.eq.2) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           !J/lam*d(p)/dX_K*q*F^-1_Ki
           !J/lam*d(q)/dX_K*p*F^-1_Ki ,d(S_ij)/dx_j
           !-d(q)/dX_K*F^-1_Ki
           aux1=(Jlam*gp(1)*shp+Jlam*p*dxi-dxi)
           aux2=(Jlam*gp(2)*shp+Jlam*p*dyi-dyi)

           B(1) = (aux1*F(1,1)+aux2*F(2,1))*dxj
           B(2) = (aux1*F(1,2)+aux2*F(2,2))*dyj
           B(3) = (aux1*F(1,1)+aux2*F(2,1))*dyj &
              & + (aux1*F(1,2)+aux2*F(2,2))*dxj

       elseif (nd.eq.3) then

           dxi = e%cartd(1,inod)
           dyi = e%cartd(2,inod)
           dzi = e%cartd(3,inod)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           aux1=(Jlam*gp(1)*shp+Jlam*p*dxi-dxi)
           aux2=(Jlam*gp(2)*shp+Jlam*p*dyi-dyi)
           aux3=(Jlam*gp(3)*shp+Jlam*p*dzi-dzi)

           B(1) = (aux1*F(1,1)+aux2*F(2,1)+aux3*F(3,1))*dxj
           B(2) = (aux1*F(1,2)+aux2*F(2,2)+aux3*F(3,2))*dyj
           B(3) = (aux1*F(1,3)+aux2*F(2,3)+aux3*F(3,3))*dzj
           B(4) = (aux1*F(1,2)+aux2*F(2,2)+aux3*F(3,2))*dzj &
              & + (aux1*F(1,3)+aux2*F(2,3)+aux3*F(3,3))*dyj
           B(5) = (aux1*F(1,1)+aux2*F(2,1)+aux3*F(3,1))*dzj &
              & + (aux1*F(1,3)+aux2*F(2,3)+aux3*F(3,3))*dxj
           B(6) = (aux1*F(1,1)+aux2*F(2,1)+aux3*F(3,1))*dyj &
              & + (aux1*F(1,2)+aux2*F(2,2)+aux3*F(3,2))*dxj

       endif

   end subroutine sld_calculate_dqFinv_vector

   subroutine sld_calculate_dqF(e,nd,inod,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                         F_iK*dq/dX_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,nd
       real(rp)            , intent(in)   :: F(nd,nd)
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: dx,dy,dz
   
       B=0.0_rp

       if (nd.eq.2) then

           dx = e%cartd(1,inod)
           dy = e%cartd(2,inod)

           B(1) = F(1,1)*dx + F(1,2)*dy
           B(2) = F(2,1)*dx + F(2,2)*dy

       elseif (nd.eq.3) then

           dx = e%cartd(1,inod)
           dy = e%cartd(2,inod)
           dz = e%cartd(3,inod)

           B(1) = F(1,1)*dx + F(1,2)*dy + F(1,3)*dz
           B(2) = F(2,1)*dx + F(2,2)*dy + F(2,3)*dz
           B(3) = F(3,1)*dx + F(3,2)*dy + F(3,3)*dz

       endif

   end subroutine sld_calculate_dqF

   subroutine sld_calculate_dqFinv(e,nd,inod,J,lam,gp,p,F,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                  d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,nd
       real(rp)            , intent(in)   :: F(nd,nd),J,lam,p,gp(nd)
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: shp,dx,dy,dz,aux1,aux2,aux3,Jlam

       B    = 0.0_rp
       shp  = e%shape(inod,e%igaus)
       Jlam = J/lam

       if (nd.eq.2) then

           dx = e%cartd(1,inod)
           dy = e%cartd(2,inod)

           aux1=(Jlam*gp(1)*shp+Jlam*p*dx-dx)
           aux2=(Jlam*gp(2)*shp+Jlam*p*dy-dy)

           B(1) = aux1*F(1,1)+aux2*F(2,1)
           B(2) = aux1*F(1,2)+aux2*F(2,2)

       elseif (nd.eq.3) then

           dx = e%cartd(1,inod)
           dy = e%cartd(2,inod)
           dz = e%cartd(3,inod)

           aux1=(Jlam*gp(1)*shp+Jlam*p*dx-dx)
           aux2=(Jlam*gp(2)*shp+Jlam*p*dy-dy)
           aux3=(Jlam*gp(3)*shp+Jlam*p*dz-dz)

           B(1) = aux1*F(1,1)+aux2*F(2,1)+aux3*F(3,1)
           B(2) = aux1*F(1,2)+aux2*F(2,2)+aux3*F(3,2)
           B(3) = aux1*F(1,3)+aux2*F(2,3)+aux3*F(3,3)

       endif

   end subroutine sld_calculate_dqFinv

   subroutine sld_calculateU_inditial(e,nd,nod,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the vector shape functions U, inditial
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,nd
       real(rp)            , intent(inout):: B(nd)
   
       B=0.0_rp
   
       if (nd .eq. 2) then
   
           B = [ e%shape(nod,e%igaus), e%shape(nod,e%igaus)]
   
       elseif (nd .eq. 3) then
   
           B = [ e%shape(nod,e%igaus), e%shape(nod,e%igaus), e%shape(nod,e%igaus)]
   
       endif
   
   end subroutine sld_calculateU_inditial

   !---------------------------Higher order derivatives----------------------
   subroutine sld_calculate_NdF(e,nd,tn,nod,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !             (xi_ij*d(F_jK)/dX_K - 1/nd*xi_ll*d(F_iK)/dX_K)
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,nd,tn
       real(rp)            , intent(in)   :: gF(nd)
       real(rp)            , intent(inout):: B(tn,nd)
       real(rp)            :: aux,dx,dy,dz
   
       B   = 0.0_rp
       aux = e%shape(nod,e%igaus)/nd

       if (nd.eq.2) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)

           !xi_ij*d(F_iK)/dX_K 
           B(1,1) = gF(1) 
           B(1,2) = 0.0_rp

           B(2,1) = 0.0_rp
           B(2,2) = gF(2)

           B(3,1) = gF(2)
           B(3,2) = gF(1)

           B = e%shape(nod,e%igaus)*B

           B(1,1) = B(1,1)-aux*gF(1)  
           B(1,2) = B(1,2)-aux*gF(2)  

           B(2,1) = B(2,1)-aux*gF(1)  
           B(2,2) = B(2,2)-aux*gF(2)  

       elseif (nd.eq.3) then

           dx = e%cartd(1,nod)
           dy = e%cartd(2,nod)
           dz = e%cartd(3,nod)

           !xi_ij*d(F_iK)/dX_K 
           B(1,1) = gF(1)
           B(1,2) = 0.0_rp
           B(1,3) = 0.0_rp

           B(2,1) = 0.0_rp
           B(2,2) = gF(2)
           B(2,3) = 0.0_rp

           B(3,1) = 0.0_rp
           B(3,2) = 0.0_rp
           B(3,3) = gF(3)

           B(4,1) = 0.0_rp
           B(4,2) = gF(3)
           B(4,3) = gF(2)

           B(5,1) = gF(3)
           B(5,2) = 0.0_rp
           B(5,3) = gF(1)

           B(6,1) = gF(2)
           B(6,2) = gF(1)
           B(6,3) = 0.0_rp

           B = e%shape(nod,e%igaus)*B

           B(1,1) = B(1,1)-aux*gF(1)  
           B(1,2) = B(1,2)-aux*gF(2)  
           B(1,3) = B(1,3)-aux*gF(3)  

           B(2,1) = B(2,1)-aux*gF(1)  
           B(2,2) = B(2,2)-aux*gF(2)  
           B(2,3) = B(2,3)-aux*gF(3)  

           B(3,1) = B(3,1)-aux*gF(1)  
           B(3,2) = B(3,2)-aux*gF(2)  
           B(3,3) = B(3,3)-aux*gF(3)  

       endif
   
   end subroutine sld_calculate_NdF

   subroutine sld_calculate_trace_NS_dF(e,nd,tn,nod,J,gdev,Jder,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, voigt
       !    The inner derivative comes precalculated as Jder
       !                  xi_lm*S_lm*d(J*F^-1_Ki)/dX_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement),intent(in)   :: e
       integer(ip)         ,intent(in)   :: nod,tn,nd
       real(rp)            ,intent(in)   :: Jder(nd),gdev(tn),J,gF(nd)
       real(rp)            ,intent(inout):: B(tn,nd)
       real(rp)                          :: shp,aux1,aux2,aux3

       B=0.0_rp
       shp = J*e%shape(nod,e%igaus)

       if (nd.eq.2) then

           aux1 = gF(1)
           aux2 = gF(2)

                   !S_lm*(d(J)/dX_K  + d(F^-1_Ki)/dX_K)
           B(1,1)= gdev(1)*(Jder(1) + aux1)
           B(1,2)= gdev(1)*(Jder(2) + aux2)

           B(2,1)= gdev(2)*(Jder(1) + aux1)
           B(2,2)= gdev(2)*(Jder(2) + aux2)

           B(3,1)= 2.0_rp*gdev(3)*(Jder(1) + aux1)
           B(3,2)= 2.0_rp*gdev(3)*(Jder(2) + aux2)

       elseif (nd.eq.3) then

           aux1 = gF(1)
           aux2 = gF(2)
           aux3 = gF(3)
                   !S_lm*(d(J)/dX_K  + d(F^-1_Ki)/dX_K)
           B(1,1)= gdev(1)*(Jder(1) + aux1)
           B(1,2)= gdev(1)*(Jder(2) + aux2)
           B(1,3)= gdev(1)*(Jder(3) + aux3)

           B(2,1)= gdev(2)*(Jder(1) + aux1)
           B(2,2)= gdev(2)*(Jder(2) + aux2)
           B(2,2)= gdev(2)*(Jder(3) + aux3)

           B(3,1)= gdev(3)*(Jder(1) + aux1)
           B(3,2)= gdev(3)*(Jder(2) + aux2)
           B(3,3)= gdev(3)*(Jder(3) + aux3)

           B(4,1)= 2.0_rp*gdev(4)*(Jder(1) + aux1)
           B(4,2)= 2.0_rp*gdev(4)*(Jder(2) + aux2)
           B(4,3)= 2.0_rp*gdev(4)*(Jder(3) + aux3)

           B(5,1)= 2.0_rp*gdev(5)*(Jder(1) + aux1)
           B(5,2)= 2.0_rp*gdev(5)*(Jder(2) + aux2)
           B(5,3)= 2.0_rp*gdev(5)*(Jder(3) + aux3)

           B(6,1)= 2.0_rp*gdev(6)*(Jder(1) + aux1)
           B(6,2)= 2.0_rp*gdev(6)*(Jder(2) + aux2)
           B(6,3)= 2.0_rp*gdev(6)*(Jder(3) + aux3)

       endif

       B = shp*B

   end subroutine sld_calculate_trace_NS_dF

   subroutine sld_calculate_qdF(e,nd,tn,nod,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                         q*dF_iK/dX_K
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,nd,tn
       real(rp)            , intent(in)   :: gF(nd)
       real(rp)            , intent(inout):: B(nd)
   
       B=0.0_rp

       if (nd.eq.2) then

           B(1) = gF(1)
           B(2) = gF(2)

       elseif (nd.eq.3) then

           B(1) = gF(1)
           B(2) = gF(2)
           B(3) = gF(3)

       endif

       B = B*e%shape(nod,e%igaus)

   end subroutine sld_calculate_qdF

   subroutine sld_calculate_qdFinv(e,nd,tn,inod,J,lam,p,Jder,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                  d/dX_K((Jp/lam)-1)q*F^(-1)_Ki
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,nd,tn
       real(rp)            , intent(in)   :: gF(nd),J,lam
       real(rp)            , intent(in)   :: Jder(nd),p
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: shp,aux1,aux2,aux3,pqlam

       B    = 0.0_rp
       shp  = e%shape(inod,e%igaus)
       pqlam = p*shp/lam

       if (nd.eq.2) then

           aux1 = gF(1)
           aux2 = gF(2)

           B(1) = pqlam*Jder(1) + J*pqlam*aux1 - shp*aux1
           B(2) = pqlam*Jder(2) + J*pqlam*aux2 - shp*aux2

       elseif (nd.eq.3) then

           aux1 = gF(1)
           aux2 = gF(2)
           aux3 = gF(3)

           B(1) = pqlam*Jder(1) + J*pqlam*aux1 - shp*aux1
           B(2) = pqlam*Jder(2) + J*pqlam*aux2 - shp*aux2
           B(3) = pqlam*Jder(3) + J*pqlam*aux3 - shp*aux3

       endif

   end subroutine sld_calculate_qdFinv

   subroutine sld_calculate_qdF_vector(e,nd,tn,inod,jnod,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !                  q_h*d(F_iK)/dX_K*d(S_ij)/dx_j
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,tn,nd
       real(rp)            , intent(in)   :: gF(nd)
       real(rp)            , intent(inout):: B(tn)
       real(rp)            :: dxj,dyj,dzj,aux1,aux2,aux3,shp
   
       B=0.0_rp
       shp  = e%shape(inod,e%igaus)

       if (nd.eq.2) then

           aux1 = gF(1)
           aux2 = gF(2)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           B(1) = aux1*dxj
           B(2) = aux2*dyj
           B(3) = aux1*dyj + aux2*dxj


       elseif (nd.eq.3) then

           aux1 = gF(1)
           aux2 = gF(2)
           aux3 = gF(3)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1) = aux1*dxj
           B(2) = aux2*dyj
           B(3) = aux3*dzj
           B(4) = aux2*dzj + aux3*dyj
           B(5) = aux1*dzj + aux3*dxj
           B(6) = aux1*dyj + aux2*dxj

       endif

       B=B*shp

   end subroutine sld_calculate_qdF_vector

   subroutine sld_calculate_qdFinv_vector(e,nd,tn,inod,jnod,J,lam,p,Jder,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes (nonlinear part) dJ/X_K and dF/dX_K
       !            -d/dX_K((Jp/lam)-1)*F^-1_Ki*q_h)*d(S_ij)/dx_j
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,tn,nd
       real(rp)            , intent(in)   :: gF(nd),J,lam,p,Jder(nd)
       real(rp)            , intent(inout):: B(tn)
       real(rp)            :: dxj,dyj,dzj,shp
       real(rp)            :: aux1,aux2,aux3,pqlam
   
       B     = 0.0_rp
       shp   = e%shape(inod,e%igaus)
       pqlam = p*shp/lam

       if (nd.eq.2) then

           !J/lam*d(F^-1_Ki)/dX_K*q*p
           !J/lam*d(F^-1_Ki)/dX_K*p*q ,d(S_ij)/dx_j
           !-d(F^-1_Ki)/dX_K*q
           aux1 = pqlam*Jder(1) + J*pqlam*gF(1) - shp*gF(1)
           aux2 = pqlam*Jder(2) + J*pqlam*gF(2) - shp*gF(2)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           B(1) = aux1*dxj
           B(2) = aux2*dyj
           B(3) = aux1*dyj + aux2*dxj

       elseif (nd.eq.3) then

           !J/lam*d(F^-1_Ki)/dX_K*q*p
           !J/lam*d(F^-1_Ki)/dX_K*p*q ,d(S_ij)/dx_j
           !-d(F^-1_Ki)/dX_K*q
           aux1 = pqlam*Jder(1) + J*pqlam*gF(1) - shp*gF(1)
           aux2 = pqlam*Jder(2) + J*pqlam*gF(2) - shp*gF(2)
           aux3 = pqlam*Jder(3) + J*pqlam*gF(3) - shp*gF(3)

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           B(1) = aux1*dxj
           B(2) = aux2*dyj
           B(3) = aux3*dzj
           B(4) = aux2*dzj + aux3*dyj
           B(5) = aux1*dzj + aux3*dxj
           B(6) = aux1*dyj + aux2*dxj

       endif

   end subroutine sld_calculate_qdFinv_vector

   subroutine sld_calculate_NdF_tensor(e,nd,tn,inod,jnod,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the elemental material matrix 
       !              (xi_ij*d(F_jK)/dX_K,div(S_iq)_i)
       !
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: gF(nd)
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            :: dxj,dyj,dzj,aux1,aux2,aux3
   
       B=0.0_rp

       if (nd.eq.2) then

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           !xi_ij*d(F_iK)/dX_K 
           aux1 = gF(1)
           aux2 = gF(2)

           B(1,1) = aux1*dxj
           B(1,2) = 0.0_rp                         
           B(1,3) = aux1*dyj

           B(2,1) = 0.0_rp                         
           B(2,2) = aux2*dyj
           B(2,3) = aux2*dxj 

           B(3,1) = aux2*dxj
           B(3,2) = aux1*dyj 
           B(3,3) = aux2*dyj + aux1*dxj

       else
           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           !xi_ij*d(F_iK)/dX_K 
           aux1 = gF(1)
           aux2 = gF(2)
           aux3 = gF(3)

           B(1,1) = aux1*dxj
           B(1,2) = 0.0_rp                         
           B(1,3) = 0.0_rp                         
           B(1,4) = 0.0_rp                         
           B(1,5) = aux1*dzj
           B(1,6) = aux1*dyj

           B(2,1) = 0.0_rp                         
           B(2,2) = aux2*dyj
           B(2,3) = 0.0_rp                         
           B(2,4) = aux2*dzj
           B(2,5) = 0.0_rp                         
           B(2,6) = aux2*dxj 

           B(3,1) = 0.0_rp                         
           B(3,2) = 0.0_rp                         
           B(3,3) = aux3*dzj
           B(3,4) = aux3*dyj
           B(3,5) = aux3*dxj 
           B(3,6) = 0.0_rp 

           B(4,1) = 0.0_rp                         
           B(4,2) = aux3*dyj
           B(4,3) = aux2*dzj                         
           B(4,4) = aux3*dzj + aux2*dyj
           B(4,5) = aux2*dxj  
           B(4,6) = aux3*dxj 

           B(5,1) = aux3*dxj
           B(5,2) = 0.0_rp 
           B(5,3) = aux1*dzj
           B(5,4) = aux1*dyj
           B(5,5) = aux3*dzj + aux1*dxj
           B(5,6) = aux3*dyj

           B(6,1) = aux2*dxj
           B(6,2) = aux1*dyj 
           B(6,3) = 0.0_rp
           B(6,4) = aux1*dzj
           B(6,5) = aux2*dzj
           B(6,6) = aux2*dyj + aux1*dxj

       endif

       B = B*e%shape(inod,e%igaus)

   end subroutine sld_calculate_NdF_tensor

   subroutine sld_calculate_NdF_tensor_trace(e,nd,tn,inod,jnod,gF,B)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !       ((xi_ll*d(F_iK)/dX_K),div(S_iq)_i)
      !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: gF(nd)
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            :: dxj,dyj,dzj,aux1,aux2,aux3
   
       B=0.0_rp

       if (nd.eq.2) then

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           aux1 = gF(1)
           aux2 = gF(2)

           B(1,1) = aux1*dxj
           B(1,2) = aux2*dyj
           B(1,3) = aux1*dyj + aux2*dxj

           B(2,1) = aux1*dxj
           B(2,2) = aux2*dyj
           B(2,3) = aux1*dyj + aux2*dxj

       else

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           aux1 = gF(1)
           aux2 = gF(2)
           aux3 = gF(3)

           B(1,1) = aux1*dxj
           B(1,2) = aux2*dyj
           B(1,3) = aux3*dzj
           B(1,4) = aux2*dzj + aux3*dyj
           B(1,5) = aux1*dzj + aux3*dxj
           B(1,6) = aux1*dyj + aux2*dxj

           B(2,1) = aux1*dxj
           B(2,2) = aux2*dyj
           B(2,3) = aux3*dzj
           B(2,4) = aux2*dzj + aux3*dyj
           B(2,5) = aux1*dzj + aux3*dxj
           B(2,6) = aux1*dyj + aux2*dxj

           B(3,1) = aux1*dxj
           B(3,2) = aux2*dyj
           B(3,3) = aux3*dzj
           B(3,4) = aux2*dzj + aux3*dyj
           B(3,5) = aux1*dzj + aux3*dxj
           B(3,6) = aux1*dyj + aux2*dxj

       endif

       B = B*e%shape(inod,e%igaus)

   end subroutine sld_calculate_NdF_tensor_trace

   subroutine sld_calculate_trace_NS_dF_tensor(e,nd,tn,inod,jnod,gdev,J,Jder,gF,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes 
       !            (xi_lm*S_ml*d(J)/dX_K*F^-1_Ki
       !            +J*xi_lm*S_ml*d(F^-1_Ki)/dX_K    ,d(S_iq)/dx_q)
       !
       !-----------------------------------------------------------------------
       implicit none
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: gF(nd),Jder(nd),gdev(tn),J
       real(rp)            , intent(inout):: B(tn,tn)
       real(rp)            :: dxj,dyj,dzj,aux1,aux2,aux3

       B=0.0_rp

       if (nd.eq.2) then

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)

           !( Jd(F^-1_Ki)/dX_K + d(J)/dX_K )
           aux1 = J*gF(1) + Jder(1)
           aux2 = J*gF(2) + Jder(2)

           B(1,1) = gdev(1)*aux1*dxj
           B(1,2) = gdev(1)*aux2*dyj
           B(1,3) = gdev(1)*aux1*dyj + aux2*dxj

           B(2,1) = gdev(2)*aux1*dxj
           B(2,2) = gdev(2)*aux2*dyj
           B(2,3) = gdev(2)*aux1*dyj + aux2*dxj

           B(3,1) = gdev(3)*aux1*dxj
           B(3,2) = gdev(3)*aux2*dyj
           B(3,3) = gdev(3)*aux1*dyj + aux2*dxj


       elseif (nd.eq.3) then

           dxj = e%cartd(1,jnod)
           dyj = e%cartd(2,jnod)
           dzj = e%cartd(3,jnod)

           !( Jd(F^-1_Ki)/dX_K + d(J)/dX_K )
           aux1 = J*gF(1) + Jder(1)
           aux2 = J*gF(2) + Jder(2)
           aux3 = J*gF(3) + Jder(3)

           B(1,1) = gdev(1)*aux1*dxj
           B(1,2) = gdev(1)*aux2*dyj
           B(1,3) = gdev(1)*aux3*dzj
           B(1,4) = gdev(1)*aux2*dzj + aux3*dyj
           B(1,5) = gdev(1)*aux1*dzj + aux3*dxj
           B(1,6) = gdev(1)*aux1*dyj + aux2*dxj

           B(2,1) = gdev(2)*aux1*dxj
           B(2,2) = gdev(2)*aux2*dyj
           B(2,3) = gdev(2)*aux3*dzj
           B(2,4) = gdev(2)*aux2*dzj + aux3*dyj
           B(2,5) = gdev(2)*aux1*dzj + aux3*dxj
           B(2,6) = gdev(2)*aux1*dyj + aux2*dxj

           B(3,1) = gdev(3)*aux1*dxj
           B(3,2) = gdev(3)*aux2*dyj
           B(3,3) = gdev(3)*aux3*dzj
           B(3,4) = gdev(3)*aux2*dzj + aux3*dyj
           B(3,5) = gdev(3)*aux1*dzj + aux3*dxj
           B(3,6) = gdev(3)*aux1*dyj + aux2*dxj

           B(4,1) = gdev(4)*aux1*dxj
           B(4,2) = gdev(4)*aux2*dyj
           B(4,3) = gdev(4)*aux3*dzj
           B(4,4) = gdev(4)*aux2*dzj + aux3*dyj
           B(4,5) = gdev(4)*aux1*dzj + aux3*dxj
           B(4,6) = gdev(4)*aux1*dyj + aux2*dxj

           B(5,1) = gdev(5)*aux1*dxj
           B(5,2) = gdev(5)*aux2*dyj
           B(5,3) = gdev(5)*aux3*dzj
           B(5,4) = gdev(5)*aux2*dzj + aux3*dyj
           B(5,5) = gdev(5)*aux1*dzj + aux3*dxj
           B(5,6) = gdev(5)*aux1*dyj + aux2*dxj

           B(6,1) = gdev(6)*aux1*dxj
           B(6,2) = gdev(6)*aux2*dyj
           B(6,3) = gdev(6)*aux3*dzj
           B(6,4) = gdev(6)*aux2*dzj + aux3*dyj
           B(6,5) = gdev(6)*aux1*dzj + aux3*dxj
           B(6,6) = gdev(6)*aux1*dyj + aux2*dxj


       endif

       B = B*e%shape(inod,e%igaus)

   end subroutine sld_calculate_trace_NS_dF_tensor

   subroutine sld_calculateJacobianDer_Fspat_ki(e,nd,tn,J,Finv,gF,derJFki)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the term (from Jacobi's Formula)
       !          d/dX_K(J)F^-1_Ki = J*F^-1_Pq*d/dX_K(F_qP)*F^-1_Ki
       ! Used in the non-linear terms of the adjoint
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nd,tn
       real(rp)            , intent(in)   :: gF(nd,tn),J
       real(rp)            , intent(in)   :: Finv(nd,nd)
       real(rp)            , intent(inout):: derJFki(nd)
       real(rp)            :: auxx,auxy,auxz
   
       derJFki=0.0_rp
   
       if (nd .eq. 2) then
   
        auxx       = gF(1,1)*Finv(1,1) + gF(1,3)*Finv(2,1) &
               &   + gF(2,1)*Finv(1,2) + gF(2,3)*Finv(2,2)

        auxy       = gF(1,3)*Finv(1,1) + gF(1,2)*Finv(2,1) &
               &   + gF(2,3)*Finv(1,2) + gF(2,2)*Finv(2,2)

        derJFki(1) = Finv(1,1)*auxx &
               &   + Finv(2,1)*auxy

        derJFki(2) = Finv(1,2)*auxx &
               &   + Finv(2,2)*auxy 

       elseif (nd .eq. 3) then

        auxx       = gF(1,1)*Finv(1,1) + gF(1,4)*Finv(2,1) + gF(1,5)*Finv(3,1)&
               &   + gF(2,1)*Finv(1,2) + gF(2,4)*Finv(2,2) + gF(2,5)*Finv(3,2)&
               &   + gF(3,1)*Finv(1,3) + gF(3,4)*Finv(2,3) + gF(3,5)*Finv(3,3)

        auxy       = gF(1,4)*Finv(1,1) + gF(1,2)*Finv(2,1) + gF(1,6)*Finv(3,1)&
               &   + gF(2,4)*Finv(1,2) + gF(2,2)*Finv(2,2) + gF(2,6)*Finv(3,2)&
               &   + gF(3,4)*Finv(1,3) + gF(3,2)*Finv(2,3) + gF(3,6)*Finv(3,3)

        auxz       = gF(1,5)*Finv(1,1) + gF(1,6)*Finv(2,1) + gF(1,3)*Finv(3,1)&
               &   + gF(2,5)*Finv(1,2) + gF(2,6)*Finv(2,2) + gF(2,3)*Finv(3,2)&
               &   + gF(3,5)*Finv(1,3) + gF(3,6)*Finv(2,3) + gF(3,3)*Finv(3,3)
   

        derJFki(1) = Finv(1,1)*auxx &
               &   + Finv(2,1)*auxy &
               &   + Finv(3,1)*auxz

        derJFki(2) = Finv(1,2)*auxx &
               &   + Finv(2,2)*auxy &
               &   + Finv(3,2)*auxz

        derJFki(3) = Finv(1,3)*auxx &
               &   + Finv(2,3)*auxy &
               &   + Finv(3,3)*auxz
       endif

       derJFki = J*derJFki
   
   end subroutine sld_calculateJacobianDer_Fspat_ki

   !---------------------------Dynamic subscales----------------------

   subroutine sld_calculateVTauTau_U(e,nd,tn,inod,jnod,tautau,densi,dtinv,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the elemental matrix 
       !            (v_h,i*(1-tau_k^-1*tau_t),rho*U^n+1,i/dt2)
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: densi,dtinv,tautau !(1-tau_k^-1*tau_t)
       real(rp)            , intent(inout):: B(nd,nd)
       real(rp)            :: shp_i,shp_j
   
       B    = 0.0_rp

       shp_i  = e%shape(inod,e%igaus)*(1.0_rp-tautau)
       shp_j  = densi*e%shape(jnod,e%igaus)*dtinv

       if (nd.eq.2) then

           B(1,1) = shp_i*shp_j
           B(2,2) = shp_i*shp_j

       elseif (nd.eq.3) then

           B(1,1) = shp_i*shp_j
           B(2,2) = shp_i*shp_j
           B(3,3) = shp_i*shp_j

       endif

   end subroutine sld_calculateVTauTau_U

   subroutine sld_calculateVTauTau_S(e,nd,tn,inod,jnod,tautau,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the elemental matrix 
       !            (v_h,i*(1-tau_k^-1*tau_t), -div(S_ij),i)
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: tautau !(1-tau_k^-1*tau_t)
       real(rp)            , intent(inout):: B(nd,tn)
       real(rp)            :: shp,dx,dy,dz
   
       B    = 0.0_rp

       !Negative sign to account for -div(S_ij)
       shp  = -e%shape(inod,e%igaus)*(1.0_rp-tautau)

       if (nd.eq.2) then

           dx = e%cartd(1,jnod)
           dy = e%cartd(2,jnod)

           B(1,1) = dx*shp
           B(1,2) = 0.0_rp
           B(1,3) = dy*shp

           B(2,1) = 0.0_rp
           B(2,2) = dy*shp
           B(2,3) = dx*shp


       elseif (nd.eq.3) then

           dx = e%cartd(1,jnod)
           dy = e%cartd(2,jnod)
           dz = e%cartd(3,jnod)

           B(1,1) = dx*shp
           B(1,2) = 0.0_rp
           B(1,3) = 0.0_rp
           B(1,4) = 0.0_rp
           B(1,5) = dz*shp
           B(1,6) = dy*shp

           B(2,1) = 0.0_rp
           B(2,2) = dy*shp
           B(2,3) = 0.0_rp
           B(2,4) = dz*shp
           B(2,5) = 0.0_rp
           B(2,6) = dx*shp

           B(3,1) = 0.0_rp
           B(3,2) = 0.0_rp
           B(3,3) = dz*shp
           B(3,4) = dy*shp
           B(3,5) = dx*shp
           B(3,6) = 0.0_rp

       endif

   end subroutine sld_calculateVTauTau_S

   subroutine sld_calculateVTauTau_P(e,nd,tn,inod,jnod,tautau,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the elemental matrix 
       !            (v_h,i*(1-tau_k^-1*tau_t), -grad(P),i)
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: inod,jnod,nd,tn
       real(rp)            , intent(in)   :: tautau !(1-tau_k^-1*tau_t)
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: shp,dx,dy,dz
   
       B    = 0.0_rp

       !Negative sign to account for -grad(P),i
       shp  = -e%shape(inod,e%igaus)*(1.0_rp-tautau)

       if (nd.eq.2) then
           dx = e%cartd(1,jnod)
           dy = e%cartd(2,jnod)

           B(1) = dx*shp
           B(2) = dy*shp

       elseif (nd.eq.3) then

           dx = e%cartd(1,jnod)
           dy = e%cartd(2,jnod)
           dz = e%cartd(3,jnod)

           B(1) = dx*shp
           B(2) = dy*shp
           B(3) = dz*shp

       endif

   end subroutine sld_calculateVTauTau_P

   subroutine sld_calculateVTauTauRHS_U(e,nd,nod,asgs,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the vector shape functions U, inditial
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,nd
       real(rp)            , intent(in)   :: asgs(nd)
       real(rp)            , intent(inout):: B(nd)
       real(rp)            :: shp
   
       B=0.0_rp
       shp  = e%shape(nod,e%igaus)
   
       if (nd .eq. 2) then
   
           B(1) = asgs(1)*shp
           B(2) = asgs(2)*shp
   
       elseif (nd .eq. 3) then
   
           B(1) = asgs(1)*shp
           B(2) = asgs(2)*shp
           B(3) = asgs(3)*shp
   
       endif
   
   end subroutine sld_calculateVTauTauRHS_U

end module Mod_SUPCauchyElement_NH_tools
