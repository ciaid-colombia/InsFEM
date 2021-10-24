subroutine matove(xmatr,vecto,ndime,ntens)
!**********************************************************************
!                                      
!**** This routine stores a symmetric matrix XMATR into a vector VECTO
!
!**********************************************************************
  use      typre
  implicit none
  integer(ip),intent(in)   :: ndime,ntens
  real(rp),intent(in)      ::xmatr(ndime,ndime)
  real(rp),intent(out)     ::vecto(ntens)

      if(ndime.eq.2) then
        vecto(1)=xmatr(1,1)
        vecto(3)=xmatr(1,2)
        vecto(2)=xmatr(2,2)
      else
        vecto(1)=xmatr(1,1)
        vecto(4)=xmatr(1,2)
        vecto(2)=xmatr(2,2)
        vecto(5)=xmatr(1,3)
        vecto(6)=xmatr(2,3)
        vecto(3)=xmatr(3,3)
      end if

end subroutine matove
