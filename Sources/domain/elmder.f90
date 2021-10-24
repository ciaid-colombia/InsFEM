subroutine elmder(nnode,ndime,deriv,elcod,cartd,detjm,xjacm,xjaci)
!-----------------------------------------------------------------------
! NAME
!    elmder
! DESCRIPTION
!    This routine calculates the Cartesian derivatives cartd and detjm.
! USES
!    invmtx
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: nnode,ndime
  real(rp),    intent(in)  :: deriv(ndime,nnode),elcod(ndime,nnode)
  real(rp),    intent(out) :: detjm,cartd(ndime,nnode)
  real(rp),    intent(out) :: xjacm(ndime,ndime)
  real(rp),    intent(out) :: xjaci(ndime,ndime)
  
  xjacm=matmul(elcod,transpose(deriv))
  call invmtx(xjacm,xjaci,detjm,ndime)
  cartd=matmul(transpose(xjaci),deriv)
  
end subroutine elmder
!-----------------------------------------------------------------------
! NOTES
! Try to remember... when life was so tender...
!        _        _           _                    _
!       | x1 x2 x3 |         | dN1/ds dN2/ds dN3/ds |
! elcod=|          |,  deriv=|                      |
!       |_y1 y2 y3_|         |_dN1/dt dN2/dt dN3/dt_|
!        _           _                        
!       | dx/ds dx/dt |            t
! xjacm=|             |=elcod*deriv
!       |_dy/ds dy/dt_|
!        _           _
!       | ds/dx ds/dy |      -1
! xjaci=|             |=xjacm  .This is because:
!       |_dt/dx dt/dy_|
!        _     _   _           _   _     _ 
!       | du/dx | | ds/dx dt/dx | | du/ds |    
!       |       |=|             | |       | and
!       |_du/dy_| |_ds/dy dt/dy_| |_du/dt_|     
!        _     _   _           _   _     _      _     _   _           _ -1  _     _ 
!       | du/ds | | dx/ds dy/ds | | du/dx |    | du/dx | | dx/ds dy/ds |   | du/ds |
!       |       |=|             | |       | => |       |=|             |   |       |
!       |_du/dt_| |_dx/dt dy/dt_| |_du/dy_|    |_du/dy_| |_dx/dt dy/dt_|   |_du/dt_| 
!                  _           _   _           _ -1
!                 | ds/dx dt/dx | | dx/ds dy/ds |      -1
!       Therefore |             |=|             |=xjacm  
!                 |_ds/dy dt/dy_| |_dx/dt dy/dt_|  
!        _                                                                               _    
!       | dN1/ds*ds/dx+dN1/dt*dt/dx  dN2/ds*ds/dx+dN2/dt*dt/dx dN3/ds*ds/dx+dN3/dt*dt/dx  |   
! cartd=|                                                                                 |   
!       |_dN1/ds*ds/dy+dN1/dt*dt/dy  dN2/ds*ds/dy+dN2/dt*dt/dy dN3/ds*ds/dy+dN3/dt*dt/dy _|
!            t
!      =xjaci *deriv
!           _     _
!          | du/dx |                 t
! Example: |       |=cartd*[U1 U2 U3]
!          |_du/dy_|        
!     
!***
!-----------------------------------------------------------------------
