module Mod_CauchyElement_tools
   use typre
   use Mod_Element
   implicit none

contains

   !---------------------------- Subroutines Form function gradients ------------------------------
   
   subroutine sld_calculateB(e,ndime,tn,nod,B) 
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, voigt
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: nod,ndime,tn
       real(rp)            , intent(inout):: B(tn,ndime)
   
       B=0.0_rp
   
       if (ndime.eq.2) then
   
           B = reshape([ e%cartd(1,nod),      0.0_rp    ,e%cartd(2,nod) ,&
                            0.0_rp     , e%cartd(2,nod) ,e%cartd(1,nod)],[3,2])
   
       elseif (ndime.eq.3) then
   
B = reshape([e%cartd(1,nod),0.0_rp,0.0_rp,0.0_rp,e%cartd(3,nod),e%cartd(2,nod), &
             0.0_rp,e%cartd(2,nod),0.0_rp,e%cartd(3,nod),0.0_rp,e%cartd(1,nod),&
             0.0_rp,0.0_rp,e%cartd(3,nod),e%cartd(2,nod),e%cartd(1,nod),0.0_rp],[6,3])
   
       endif
   
   end subroutine sld_calculateB

   subroutine sld_calculateBt(e,ndime,tn,nod,B) 
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, voigt
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in) :: e
       integer(ip)         , intent(in)   :: nod,ndime,tn
       real(rp)            , intent(inout):: B(ndime,tn)
   
       B=0.0_rp
   
       if (ndime.eq.2) then
   
           B = reshape([ e%cartd(1,nod),     0.0_rp    ,&
                         0.0_rp        , e%cartd(2,nod),&
                         e%cartd(2,nod), e%cartd(1,nod)],[2,3])
   
       elseif (ndime.eq.3) then
   
   B = reshape([ e%cartd(1,nod),0.0_rp        ,0.0_rp,& 
                 0.0_rp        ,e%cartd(2,nod),0.0_rp,& 
                 0.0_rp        ,0.0_rp        ,e%cartd(3,nod),& 
                 0.0_rp        ,e%cartd(3,nod),e%cartd(2,nod),&
                 e%cartd(3,nod),0.0_rp        ,e%cartd(1,nod),&
                 e%cartd(2,nod),e%cartd(1,nod),0.0_rp],[3,6])
   
       endif

   end subroutine sld_calculateBt
   
   subroutine sld_calculateB_inditial(e,ndime,nod,B)
       !-----------------------------------------------------------------------
       !
       ! This routine computes the gradient of the shape functions B, inditial
       !
       !-----------------------------------------------------------------------
       implicit none
   
       class(FiniteElement), intent(in)   :: e
       integer(ip)         , intent(in)   :: nod,ndime
       real(rp)            , intent(inout):: B(ndime)
   
       B=0.0_rp
   
       if (ndime.eq.2) then
   
           B = [ e%cartd(1,nod), e%cartd(2,nod)]
   
       elseif (ndime.eq.3) then
   
           B = [ e%cartd(1,nod), e%cartd(2,nod), e%cartd(3,nod) ]
   
       endif
   
   end subroutine sld_calculateB_inditial 

end module Mod_CauchyElement_tools

