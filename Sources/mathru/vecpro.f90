subroutine vecpro(v1,v2,v3,n)

!-----------------------------------------------------------------------
!
! Two and three-dimensional vectorial product of two vectors  v3 = v1 x v2.
! The same pointer as for v1 or v2 may be used for v3. If N = 2, it is
!  assumed that v1 = (0,0,v1_3) and v2 = (v2_1,v2_2,0).      
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  :: n
  real(rp),    intent(in)  :: v2(n),v1(3)
  real(rp),    intent(out) :: v3(n)
  real(rp)                 :: c1,c2,c3

  if(n==2) then
     c1=-v1(3)*v2(2)
     c2= v1(3)*v2(1)
     v3(1)=c1
     v3(2)=c2
  else if(n==3) then
     c1=v1(2)*v2(3)-v1(3)*v2(2)
     c2=v1(3)*v2(1)-v1(1)*v2(3)
     c3=v1(1)*v2(2)-v1(2)*v2(1)
     v3(1)=c1
     v3(2)=c2
     v3(3)=c3
  end if
  
end subroutine vecpro
