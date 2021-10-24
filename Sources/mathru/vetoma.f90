subroutine vetoma(vecto,xmatr,ndime,ntens)
!**********************************************************************
!                                      
!**** This routine stores a vector VECTO as a symmetric matrix XMATR
!
!**********************************************************************
  use      typre
  implicit none
  integer(ip),intent(in)   :: ndime,ntens
  real(rp),intent(out)     ::xmatr(ndime,ndime)
  real(rp),intent(in)      ::vecto(ntens)

  if(ndime.eq.2) then
        xmatr(1,1)=vecto(1)
        xmatr(1,2)=vecto(3)
        xmatr(2,1)=vecto(3)
        xmatr(2,2)=vecto(2)
   else
        xmatr(1,1)=vecto(1)
        xmatr(1,2)=vecto(4)
        xmatr(1,3)=vecto(5)
        xmatr(2,1)=vecto(4)
        xmatr(2,2)=vecto(2)
        xmatr(2,3)=vecto(6)
        xmatr(3,1)=vecto(5)
        xmatr(3,2)=vecto(6)
        xmatr(3,3)=vecto(3)
   end if

end subroutine vetoma