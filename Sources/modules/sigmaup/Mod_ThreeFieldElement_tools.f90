module Mod_ThreeFieldElement_tools
   use typre
   use Mod_Element
   implicit none

contains

   subroutine elmbst(e,nd,tn,aux,aux2,elmat)

      class(FiniteElement)       :: e
      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: aux,aux2
      real(rp),    intent(inout) :: elmat(tn,e%mnode,tn,e%mnode)
      integer(ip)                :: inode,jnode

      elmat=0.0_rp

      if (nd.eq.2) then

      do jnode=1,e%pnode    
         do inode=1,e%pnode
      
            elmat(1,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)+aux2*(e%cartd(1,inode)*e%cartd(1,jnode)) + elmat(1,inode,1,jnode)
            elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
            elmat(1,inode,3,jnode) = aux2*(e%cartd(1,inode)*e%cartd(2,jnode)) + elmat(1,inode,3,jnode)
            elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
            elmat(2,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)+aux2*e%cartd(2,inode)*e%cartd(2,jnode) + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = aux2*e%cartd(2,inode)*e%cartd(1,jnode) + elmat(2,inode,3,jnode)
            elmat(3,inode,1,jnode) = aux2*e%cartd(2,inode)*e%cartd(1,jnode) + elmat(3,inode,1,jnode)
            elmat(3,inode,2,jnode) = aux2*e%cartd(1,inode)*e%cartd(2,jnode) + elmat(3,inode,2,jnode)                                         
            elmat(3,inode,3,jnode) = 2.0_rp*aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) + &
                  aux2*(e%cartd(1,inode)*e%cartd(1,jnode) + e%cartd(2,inode)*e%cartd(2,jnode)) + elmat(3,inode,3,jnode)
            

         end do
      end do   

       elseif (nd.eq.3) then

      do jnode=1,e%pnode
         do inode=1,e%pnode
       
            elmat(1,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)+aux2*e%cartd(1,inode)*e%cartd(1,jnode) + elmat(1,inode,1,jnode)
            elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
            elmat(1,inode,3,jnode) = 0.0_rp + elmat(1,inode,3,jnode)
            elmat(1,inode,4,jnode) = 0.0_rp + elmat(1,inode,4,jnode)
            elmat(1,inode,5,jnode) = aux2*e%cartd(1,inode)*e%cartd(3,jnode) + elmat(1,inode,5,jnode)
            elmat(1,inode,6,jnode) = aux2*e%cartd(1,inode)*e%cartd(2,jnode) + elmat(1,inode,6,jnode)
            
            elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)                                  
            elmat(2,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)+aux2*e%cartd(2,inode)*e%cartd(2,jnode) + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = 0.0_rp + elmat(2,inode,3,jnode)
            elmat(2,inode,4,jnode) = aux2*e%cartd(2,inode)*e%cartd(3,jnode) + elmat(2,inode,4,jnode)
            elmat(2,inode,5,jnode) = 0.0_rp + elmat(2,inode,5,jnode)
            elmat(2,inode,6,jnode) = aux2*e%cartd(2,inode)*e%cartd(1,jnode) + elmat(2,inode,6,jnode)

            elmat(3,inode,1,jnode) = 0.0_rp + elmat(3,inode,1,jnode)
            elmat(3,inode,2,jnode) = 0.0_rp + elmat(3,inode,2,jnode)
            elmat(3,inode,3,jnode) = aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)+aux2*e%cartd(3,inode)*e%cartd(3,jnode) + elmat(3,inode,3,jnode)
            elmat(3,inode,4,jnode) = aux2*e%cartd(3,inode)*e%cartd(2,jnode) + elmat(3,inode,4,jnode)
            elmat(3,inode,5,jnode) = aux2*e%cartd(3,inode)*e%cartd(1,jnode) + elmat(3,inode,5,jnode)
            elmat(3,inode,6,jnode) = 0.0_rp + elmat(3,inode,6,jnode)
            
            elmat(4,inode,1,jnode) = 0.0_rp + elmat(4,inode,1,jnode) 
            elmat(4,inode,2,jnode) = aux2*e%cartd(3,inode)*e%cartd(2,jnode) + elmat(4,inode,2,jnode)
            elmat(4,inode,3,jnode) = aux2*e%cartd(2,inode)*e%cartd(3,jnode) + elmat(4,inode,3,jnode)
            elmat(4,inode,4,jnode) = 2.0_rp*aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) + &
                                       aux2*(e%cartd(3,inode)*e%cartd(3,jnode) + e%cartd(2,inode)*e%cartd(2,jnode)) &
                                       + elmat(4,inode,4,jnode)
            elmat(4,inode,5,jnode) = aux2*e%cartd(2,inode)*e%cartd(1,jnode) + elmat(4,inode,5,jnode)
            elmat(4,inode,6,jnode) = aux2*e%cartd(3,inode)*e%cartd(1,jnode) + elmat(4,inode,6,jnode)
                              
            elmat(5,inode,1,jnode) = aux2*e%cartd(3,inode)*e%cartd(1,jnode) + elmat(5,inode,1,jnode)
            elmat(5,inode,2,jnode) = 0.0_rp + elmat(5,inode,2,jnode)
            elmat(5,inode,3,jnode) = aux2*e%cartd(1,inode)*e%cartd(3,jnode) + elmat(5,inode,3,jnode)
            elmat(5,inode,4,jnode) = aux2*e%cartd(1,inode)*e%cartd(2,jnode) + elmat(5,inode,4,jnode)
            elmat(5,inode,5,jnode) = 2.0_rp*aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) + &
                                       aux2*(e%cartd(1,inode)*e%cartd(1,jnode) + e%cartd(3,inode)*e%cartd(3,jnode)) &
                                       + elmat(5,inode,5,jnode)    
            elmat(5,inode,6,jnode) = aux2*e%cartd(3,inode)*e%cartd(2,jnode) + elmat(5,inode,6,jnode)  
            
            elmat(6,inode,1,jnode) = aux2*e%cartd(2,inode)*e%cartd(1,jnode) + elmat(6,inode,1,jnode)
            elmat(6,inode,2,jnode) = aux2*e%cartd(1,inode)*e%cartd(2,jnode) + elmat(6,inode,2,jnode)
            elmat(6,inode,3,jnode) = 0.0_rp + elmat(6,inode,3,jnode)
            elmat(6,inode,4,jnode) = aux2*e%cartd(1,inode)*e%cartd(3,jnode) + elmat(6,inode,4,jnode)
            elmat(6,inode,5,jnode) = aux2*e%cartd(2,inode)*e%cartd(3,jnode) + elmat(6,inode,5,jnode)
            elmat(6,inode,6,jnode) = 2.0_rp*aux*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) + &
                                    aux2*(e%cartd(1,inode)*e%cartd(1,jnode) + e%cartd(2,inode)*e%cartd(2,jnode)) &
                                    + elmat(6,inode,6,jnode)                                  

         end do
      end do   

       endif

   end subroutine elmbst

   subroutine elmbut(e,nd,tn,tmp,aux,aux2,acden,gpadv,elmat)  
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in constitutive equation
    ! -(gra_sym(u),T) +tau3(gra_sym(u),1/2mu*T) -tau1(Div(T),rho*a*grad(u)+rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: gpadv(e%pnode),acden
      real(rp),    intent(in)    :: aux,tmp,aux2
      real(rp),    intent(inout) :: elmat(tn,e%mnode,nd,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1

      elmat = 0.0_rp

      if (nd.eq.2) then

          do jnode=1,e%pnode  
              aux1=acden*gpadv(jnode) 
              do inode=1,e%pnode     

                  elmat(1,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(1,jnode) &              
                      + aux2*e%cartd(1,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &  
                      + elmat(1,inode,1,jnode)                                      

                  elmat(1,inode,2,jnode) = 0.0_rp &
                      + elmat(1,inode,2,jnode)

                  elmat(2,inode,1,jnode) = 0.0_rp &
                      + elmat(2,inode,1,jnode) 

                  elmat(2,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(2,jnode) &
                      + aux2*e%cartd(2,inode)*(aux1+e%shape(jnode,e%igaus)*tmp) &
                      + elmat(2,inode,2,jnode)

                  elmat(3,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(2,jnode) &
                      + aux2*e%cartd(2,inode)*(aux1+e%shape(jnode,e%igaus)*tmp) &
                      + elmat(3,inode,1,jnode)

                  elmat(3,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(1,jnode) & 
                      + aux2*e%cartd(1,inode)*(aux1+e%shape(jnode,e%igaus)*tmp) &              
                      + elmat(3,inode,2,jnode)                                         

              end do  
          end do


      elseif (nd.eq.3) then

          do jnode=1,e%pnode
              aux1=acden*gpadv(jnode)       
              do inode=1,e%pnode

                  elmat(1,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(1,jnode) &
                      + aux2*e%cartd(1,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &            
                      + elmat(1,inode,1,jnode)
                  elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
                  elmat(1,inode,3,jnode) = 0.0_rp + elmat(1,inode,3,jnode)

                  elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
                  elmat(2,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(2,jnode) & 
                      + aux2*e%cartd(2,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(2,inode,2,jnode)  
                  elmat(2,inode,3,jnode) = 0.0_rp + elmat(2,inode,3,jnode)

                  elmat(3,inode,1,jnode) = 0.0_rp + elmat(3,inode,1,jnode)
                  elmat(3,inode,2,jnode) = 0.0_rp + elmat(3,inode,2,jnode)
                  elmat(3,inode,3,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(3,jnode) &
                      + aux2*e%cartd(3,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(3,inode,3,jnode) 

                  elmat(4,inode,1,jnode) = 0.0_rp + elmat(4,inode,1,jnode)
                  elmat(4,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(3,jnode) &
                      + aux2*e%cartd(3,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(4,inode,2,jnode)
                  elmat(4,inode,3,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(2,jnode) & 
                      + aux2*e%cartd(2,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(4,inode,3,jnode)

                  elmat(5,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(3,jnode) & 
                      + aux2*e%cartd(3,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(5,inode,1,jnode)
                  elmat(5,inode,2,jnode) = 0.0_rp + elmat(5,inode,2,jnode)
                  elmat(5,inode,3,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(1,jnode) &
                      + aux2*e%cartd(1,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(5,inode,3,jnode)

                  elmat(6,inode,1,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(2,jnode) &
                      + aux2*e%cartd(2,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(6,inode,1,jnode)
                  elmat(6,inode,2,jnode) = aux*e%shape(inode,e%igaus)*e%cartd(1,jnode) &
                      + aux2*e%cartd(1,inode)*(aux1 + e%shape(jnode,e%igaus)*tmp) &
                      + elmat(6,inode,2,jnode)
                  elmat(6,inode,3,jnode) = 0.0_rp + elmat(6,inode,3,jnode)                                  

              end do  
          end do

  endif

   end subroutine elmbut

   subroutine elmbpt(e,nd,tn,aux,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the third lhs term for ASGS in constitutive equation
    !    -tau1*(gra(p),div(t))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: aux
      real(rp),    intent(inout) :: elmat(tn,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode

      elmat = 0.0_rp

      if (nd.eq.2) then

          do jnode=1,e%pnode  
              do inode=1,e%pnode    

                  elmat(1,inode,1,jnode) = e%cartd(1,inode)*aux*e%cartd(1,jnode) + elmat(1,inode,1,jnode)
                  elmat(2,inode,1,jnode) = e%cartd(2,inode)*aux*e%cartd(2,jnode) + elmat(2,inode,1,jnode)
                  elmat(3,inode,1,jnode) = (e%cartd(2,inode)*e%cartd(1,jnode) + e%cartd(1,inode)*e%cartd(2,jnode))*aux &
                      + elmat(3,inode,1,jnode)

              end do  
   end do

          elseif (nd.eq.3) then

              do jnode=1,e%pnode
                  do inode=1,e%pnode   

                      elmat(1,inode,1,jnode) = e%cartd(1,inode)*aux*e%cartd(1,jnode) + elmat(1,inode,1,jnode)  
                      elmat(2,inode,1,jnode) = e%cartd(2,inode)*aux*e%cartd(2,jnode) + elmat(2,inode,1,jnode)
                      elmat(3,inode,1,jnode) = e%cartd(3,inode)*aux*e%cartd(3,jnode) + elmat(3,inode,1,jnode)   
                      elmat(4,inode,1,jnode) = (e%cartd(3,inode)*e%cartd(2,jnode) + e%cartd(2,inode)*e%cartd(3,jnode))*aux &
                          + elmat(4,inode,1,jnode) 
                      elmat(5,inode,1,jnode) = (e%cartd(3,inode)*e%cartd(1,jnode) + e%cartd(1,inode)*e%cartd(3,jnode))*aux &
                          + elmat(5,inode,1,jnode)  
                      elmat(6,inode,1,jnode) = (e%cartd(2,inode)*e%cartd(1,jnode) + e%cartd(1,inode)*e%cartd(2,jnode))*aux &
                          + elmat(6,inode,1,jnode)   

                  end do  
              end do
  endif

   end subroutine elmbpt   

   subroutine elmbsv(e,nd,tn,aux,timom,dvolu,acden,gpadv,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in momentum equation
    !    (gra_sym(v),S) -tau3(gra_sym(v),1/2mu*S) -tau1(rho*a*grad(v),div(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: acden,gpadv(e%pnode),aux,timom,dvolu
      real(rp),    intent(inout) :: elmat(nd,e%mnode,tn,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1

      elmat = 0.0_rp

      if(nd .eq. 2) then

          do inode=1,e%pnode
              aux1 = -acden*gpadv(inode)*timom*dvolu        
              do jnode=1,e%pnode

                  elmat(1,inode,1,jnode) = e%cartd(1,inode)*aux*e%shape(jnode,e%igaus) + e%cartd(1,jnode)*aux1 & 
                      + elmat(1,inode,1,jnode)                                      
                  elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
                  elmat(1,inode,3,jnode) = e%cartd(2,inode)*aux*e%shape(jnode,e%igaus) + e%cartd(2,jnode)*aux1 & 
                      + elmat(1,inode,3,jnode)

                  elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
                  elmat(2,inode,2,jnode) = e%cartd(2,inode)*aux*e%shape(jnode,e%igaus) + e%cartd(2,jnode)*aux1 & 
                      + elmat(2,inode,2,jnode)
                  elmat(2,inode,3,jnode) = e%cartd(1,inode)*aux*e%shape(jnode,e%igaus) + e%cartd(1,jnode)*aux1 & 
                      + elmat(2,inode,3,jnode)                                    

              end do  
          end do

  elseif(nd .eq. 3) then

          do inode=1,e%pnode
              aux1 = -acden*gpadv(inode)*timom*dvolu         
              do jnode=1,e%pnode

                  elmat(1,inode,1,jnode) = e%cartd(1,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(1,jnode)*aux1 & 
                      + elmat(1,inode,1,jnode)
                  elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
                  elmat(1,inode,3,jnode) = 0.0_rp + elmat(1,inode,3,jnode)
                  elmat(1,inode,4,jnode) = 0.0_rp + elmat(1,inode,4,jnode)   
                  elmat(1,inode,5,jnode) = e%cartd(3,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(3,jnode)*aux1 &
                      + elmat(1,inode,5,jnode)
                  elmat(1,inode,6,jnode) = e%cartd(2,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(2,jnode)*aux1 &
                      + elmat(1,inode,6,jnode)

                  elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
                  elmat(2,inode,2,jnode) = e%cartd(2,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(2,jnode)*aux1 & 
                      + elmat(2,inode,2,jnode)
                  elmat(2,inode,3,jnode) = 0.0_rp + elmat(2,inode,3,jnode)
                  elmat(2,inode,4,jnode) = e%cartd(3,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(3,jnode)*aux1 & 
                      + elmat(2,inode,4,jnode)  
                  elmat(2,inode,5,jnode) = 0.0_rp + elmat(2,inode,5,jnode)
                  elmat(2,inode,6,jnode) = e%cartd(1,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(1,jnode)*aux1 &
                      + elmat(2,inode,6,jnode)  

                  elmat(3,inode,1,jnode) = 0.0_rp + elmat(3,inode,1,jnode)
                  elmat(3,inode,2,jnode) = 0.0_rp + elmat(3,inode,2,jnode)
                  elmat(3,inode,3,jnode) = e%cartd(3,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(3,jnode)*aux1 &
                      + elmat(3,inode,3,jnode) 
                  elmat(3,inode,4,jnode) = e%cartd(2,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(2,jnode)*aux1 &
                      + elmat(3,inode,4,jnode)
                  elmat(3,inode,5,jnode) = e%cartd(1,inode)*(aux)*e%shape(jnode,e%igaus) + e%cartd(1,jnode)*aux1 &
                      + elmat(3,inode,5,jnode)
                  elmat(3,inode,6,jnode) = 0.0_rp + elmat(3,inode,6,jnode)                                 

              end do  
          end do

  endif
 
   end subroutine elmbsv  

   subroutine elmbuv2(e,nd,aux,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the second lhs term for ASGS in momentum equation
    !    tau3*(gra_sym(v),gra_sym(u))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip)                :: nd
      real(rp),    intent(in)    :: aux
      real(rp),    intent(inout) :: elmat(nd,e%mnode,nd,e%mnode)
      integer(ip)                :: inode,jnode

      elmat = 0.0_rp
      
  if(nd.eq.2) then
      do inode=1,e%pnode         
          do jnode=1,e%pnode

              elmat(1,inode,1,jnode)= aux*(e%cartd(1,inode)*e%cartd(1,jnode)+0.5_rp*e%cartd(2,inode)*e%cartd(2,jnode)) &
                  + elmat(1,inode,1,jnode)            
              elmat(1,inode,2,jnode)= aux*(0.5_rp*e%cartd(2,inode)*e%cartd(1,jnode)) + elmat(1,inode,2,jnode)         
              elmat(2,inode,1,jnode)= aux*(0.5_rp*e%cartd(1,inode)*e%cartd(2,jnode)) + elmat(2,inode,1,jnode)
              elmat(2,inode,2,jnode)= aux*(0.5_rp*e%cartd(1,inode)*e%cartd(1,jnode)+e%cartd(2,inode)*e%cartd(2,jnode)) &
                  + elmat(2,inode,2,jnode)
          end do  
      end do

  elseif(nd.eq.3) then

      do jnode=1,e%pnode
          do inode=1,e%pnode

              elmat(1,inode,1,jnode)= aux*(e%cartd(1,inode)*e%cartd(1,jnode)+0.5_rp*e%cartd(2,inode)*e%cartd(2,jnode) &
                  + 0.5_rp*e%cartd(3,inode)*e%cartd(3,jnode)) + elmat(1,inode,1,jnode)                                  
              elmat(1,inode,2,jnode)= aux*(0.5_rp*e%cartd(2,inode)*e%cartd(1,jnode)) + elmat(1,inode,2,jnode)
              elmat(1,inode,3,jnode)= aux*(0.5_rp*e%cartd(3,inode)*e%cartd(1,jnode)) + elmat(1,inode,3,jnode)

              elmat(2,inode,1,jnode)= aux*(0.5_rp*e%cartd(1,inode)*e%cartd(2,jnode)) + elmat(2,inode,1,jnode)
              elmat(2,inode,2,jnode)= aux*(0.5_rp*e%cartd(1,inode)*e%cartd(1,jnode) + e%cartd(2,inode)*e%cartd(2,jnode) &
                  + 0.5_rp*e%cartd(3,inode)*e%cartd(3,jnode)) + elmat(2,inode,2,jnode)                                
              elmat(2,inode,3,jnode)= aux*(0.5_rp*e%cartd(3,inode)*e%cartd(2,jnode)) + elmat(2,inode,3,jnode)

              elmat(3,inode,1,jnode)= aux*(0.5_rp*e%cartd(1,inode)*e%cartd(3,jnode)) + elmat(3,inode,1,jnode)
              elmat(3,inode,2,jnode)= aux*(0.5_rp*e%cartd(2,inode)*e%cartd(3,jnode)) + elmat(3,inode,2,jnode)                                
              elmat(3,inode,3,jnode)= aux*(0.5_rp*e%cartd(1,inode)*e%cartd(1,jnode) + 0.5_rp*e%cartd(2,inode)*e%cartd(2,jnode) &
                  + e%cartd(3,inode)*e%cartd(3,jnode)) + elmat(3,inode,3,jnode)            


          end do  
      end do

  endif
 
   end subroutine elmbuv2   

   subroutine elmbsq(e,nd,tn,aux,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs term for ASGS in momentum equation
    !    tau1*(gra(q),div(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: aux
      real(rp),    intent(inout) :: elmat(1,e%mnode,tn,e%mnode)
      integer(ip)                :: inode,jnode

      elmat = 0.0_rp

      if(nd.eq.2) then
          do jnode=1,e%pnode      
              do inode=1,e%pnode   

                  elmat(1,inode,1,jnode) = e%cartd(1,inode)*(aux)*e%cartd(1,jnode) + elmat(1,inode,1,jnode)
                  elmat(1,inode,2,jnode) = e%cartd(2,inode)*(aux)*e%cartd(2,jnode) + elmat(1,inode,2,jnode)
                  elmat(1,inode,3,jnode) = (e%cartd(1,inode)*e%cartd(2,jnode) + e%cartd(2,inode)*e%cartd(1,jnode))*aux &
                      + elmat(1,inode,3,jnode)                      

              end do  
          end do

  elseif(nd.eq.3) then
          do jnode=1,e%pnode
              do inode=1,e%pnode   

                  elmat(1,inode,1,jnode) = e%cartd(1,inode)*(aux)*e%cartd(1,jnode) + elmat(1,inode,1,jnode)
                  elmat(1,inode,2,jnode) = e%cartd(2,inode)*(aux)*e%cartd(2,jnode) + elmat(1,inode,2,jnode)
                  elmat(1,inode,3,jnode) = e%cartd(3,inode)*(aux)*e%cartd(3,jnode) + elmat(1,inode,3,jnode) 
                  elmat(1,inode,4,jnode) = (e%cartd(2,inode)*e%cartd(3,jnode) + e%cartd(3,inode)*e%cartd(2,jnode))*aux & 
                      + elmat(1,inode,4,jnode)  
                  elmat(1,inode,5,jnode) = (e%cartd(1,inode)*e%cartd(3,jnode) + e%cartd(3,inode)*e%cartd(1,jnode))*aux &
                      + elmat(1,inode,5,jnode)                     
                  elmat(1,inode,6,jnode) = (e%cartd(1,inode)*e%cartd(2,jnode) + e%cartd(2,inode)*e%cartd(1,jnode))*aux &
                      + elmat(1,inode,6,jnode)                    

              end do  
          end do

  endif

   end subroutine elmbsq     
  
   subroutine elmrhc(e,nd,tn,aux,tisig,acvis,dvolu,elext,elextS,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: aux,dvolu,tisig,acvis
      real(rp),    intent(in)    :: elext(nd),elextS(tn)
      real(rp),    intent(inout) :: elrhs(tn,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux1

      elrhs = 0.0_rp

      if(nd.eq.2) then

          do inode=1,e%pnode
              aux1= (e%shape(inode,e%igaus)*(1.0_rp - tisig/(2.0_rp*acvis))*dvolu)

              elrhs(1,inode) = e%cartd(1,inode)*elext(1)*aux + aux1*elextS(1) &
                  + elrhs(1,inode) 
              elrhs(2,inode) = e%cartd(2,inode)*elext(2)*aux + aux1*elextS(2) &
                  + elrhs(2,inode)
              elrhs(3,inode) = (e%cartd(1,inode)*elext(2) + e%cartd(2,inode)*elext(1))*aux & 
                  + 2.0_rp*aux1*elextS(3) &
                  + elrhs(3,inode)

          end do

  elseif(nd.eq.3) then

          do inode=1,e%pnode

              elrhs(1,inode) = e%cartd(1,inode)*(elext(1))*aux + elrhs(1,inode) 
              elrhs(2,inode) = e%cartd(2,inode)*(elext(2))*aux + elrhs(2,inode)
              elrhs(3,inode) = e%cartd(3,inode)*(elext(3))*aux + elrhs(3,inode)
              elrhs(4,inode) = (e%cartd(2,inode)*elext(3) + e%cartd(3,inode)*elext(2))*aux & 
                  + elrhs(4,inode)
              elrhs(5,inode) = (e%cartd(1,inode)*elext(3) + e%cartd(3,inode)*elext(1))*aux & 
                  + elrhs(5,inode)
              elrhs(6,inode) = (e%cartd(1,inode)*elext(2) + e%cartd(2,inode)*elext(1))*aux & 
                  + elrhs(6,inode)

          end do
  endif

   end subroutine elmrhc    

   subroutine elmrhu(e,nd,tn,aux,aux1,elextC,elextS,elrhs)  
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: nd,tn
      real(rp),    intent(in)    :: aux,aux1,elextC(1),elextS(tn)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode

      elrhs = 0.0_rp

      if(nd.eq.2) then

          do inode=1,e%pnode

              elrhs(1,inode) = e%cartd(1,inode)*elextC(1)*aux &
                  - (e%cartd(1,inode)*elextS(1) + e%cartd(2,inode)*elextS(3))*aux1 &
                  + elrhs(1,inode) 
              elrhs(2,inode) = e%cartd(2,inode)*elextC(1)*aux &
                  - (e%cartd(2,inode)*elextS(2) + e%cartd(1,inode)*elextS(3))*aux1 &           
                  + elrhs(2,inode)

          end do

      elseif(nd.eq.3) then

          do inode=1,e%pnode

              elrhs(1,inode) = e%cartd(1,inode)*elextC(1)*aux &
                  - (e%cartd(1,inode)*elextS(1) + e%cartd(2,inode)*elextS(6) + e%cartd(3,inode)*elextS(5))*aux1 &
                  + elrhs(1,inode) 
              elrhs(2,inode) = e%cartd(2,inode)*elextC(1)*aux &
                  - (e%cartd(1,inode)*elextS(6) + e%cartd(2,inode)*elextS(2) + e%cartd(3,inode)*elextS(4))*aux1 &         
                  + elrhs(2,inode)               
              elrhs(3,inode) = e%cartd(3,inode)*elextC(1)*aux &
                  - (e%cartd(1,inode)*elextS(5) + e%cartd(2,inode)*elextS(4) + e%cartd(3,inode)*elextS(3))*aux1 &           
                  + elrhs(3,inode)               

          end do

  endif

   end subroutine elmrhu    

end module Mod_ThreeFieldElement_tools
