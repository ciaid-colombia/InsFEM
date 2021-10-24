subroutine shape1(s,nnode,shape,deriv)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions and their first derivates 
! for 1-d continuos with 2 3 4 & 5 nodes
!
!-----------------------------------------------------------------------

  use      typre
  implicit none
  integer(ip), intent(in)  :: nnode
  real(rp),    intent(in)  :: s
  real(rp),    intent(out) :: deriv(nnode), shape(nnode)
  real                     :: s1, s2, s3, s4, a, c

  if(nnode==2) then     
     shape(1)= 0.5_rp*(1.0_rp-s)              !  1    2
     shape(2)= 0.5_rp*(1.0_rp+s)
     deriv(1)=-0.5_rp
     deriv(2)= 0.5_rp
  else if(nnode==3) then
     shape(1)= 0.5_rp*s*(s-1.0_rp)            !  1    3    2
     shape(2)= 0.5_rp*s*(s+1.0_rp)
     shape(3)=-(s+1.0_rp)*(s-1.0_rp)
     deriv(1)= s-0.5_rp
     deriv(2)= s+0.5_rp
     deriv(3)=-2.0_rp*s
  else if(nnode==4) then
     a = 9.0_rp/16.0_rp
     s1= s + 1.0_rp
     s2= s - 1.0_rp     
     s3= s + 1.0_rp/3.0_rp
     s4= s - 1.0_rp/3.0_rp     
     c = 3.0_rp*a
     shape(1)= -a*s2*s3*s4                    !  1    3    4    2
     shape(2)=  a*s1*s3*s4
     shape(3)=  c*s1*s2*s4
     shape(4)= -c*s1*s2*s3
     deriv(1)= - a*(s2*s3+s2*s4+s3*s4)
     deriv(2)=   a*(s1*s3+s1*s4+s3*s4)
     deriv(3)=   c*(s1*s2+s1*s4+s2*s4)
     deriv(4)= - c*(s1*s2+s1*s3+s2*s3)
! * * * * * * * * * * * * * * * * * * 
! *Valores calculados para nnode=5  *
! * * * * * * * * * * * * * * * * * * 
  else if(nnode==5) then  
     s1 = s + 1.0_rp     
     s2 = s - 1.0_rp
     s3 = s + 1.0_rp/2.0_rp
     s4 = s - 1.0_rp/2.0_rp
     shape(1) = 2.0_rp/3.0_rp*s*s3*s4*s2     !  1    3    4    5    2
     shape(2) =-8.0_rp/3.0_rp*s*s1*s4*s2 
     shape(3) = 4.0_rp*s1*s3*s4*s2 
     shape(4) =-8.0_rp/3.0_rp*s*s1*s3*s2
     shape(5) = 2.0_rp/3.0_rp*s*s1*s3*s4
     deriv(1) = 2.0_rp/3.0_rp*(s*s3*s4+s*s3*s2+s*s4*s2+s3*s4*s2)
     deriv(2) =-8.0_rp/3.0_rp*(s*s1*s4+s*s1*s2+s*s4*s2+s1*s4*s2)
     deriv(3) =        4.0_rp*(s1*s3*s4+s1*s3*s2+s1*s4*s2+s3*s4*s2)
     deriv(4) =-8.0_rp/3.0_rp*(s*s1*s3+s*s1*s2+s*s3*s2+s1*s3*s2)
     deriv(5) = 2.0_rp/3.0_rp*(s*s1*s3+s*s1*s4+s*s3*s4+s1*s3*s4)

  end if
!  a= 0.0_rp 
!  c= 0.0_rp
!  a = a + shape(1)+ shape(2)+ shape(3)+ shape(4)
!  c = c + deriv(1)+ deriv(2)+ deriv(3)+ deriv(4)
 
end subroutine shape1
      
