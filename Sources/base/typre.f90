module typre
!-----------------------------------------------------------------------
!****f* typre
! NAME
!    typre
! DESCRIPTION
!    This module contains kind and type definitions and some commonly
!    used constants and parameters.
!***
!-----------------------------------------------------------------------
  integer, parameter  :: ip = 4    ! Integer precision
  integer, parameter  :: rp = 8    ! Real precision
  integer, parameter  :: lg = 1    ! logical precision

  ! Types
  type i1p
     integer(ip), pointer :: l(:) => NULL()
  end type i1p
  type i2p
     integer(ip), pointer :: l(:,:) => NULL()
  end type i2p
  type r1p
     real(rp),    pointer, contiguous :: a(:) => NULL()
  end type r1p
  type r2p
     real(rp),    pointer, contiguous :: a(:,:) => NULL()
  end type r2p
  type r3p
     real(rp),    pointer, contiguous :: a(:,:,:) => NULL()
  end type r3p
  type r4p
     real(rp),    pointer, contiguous :: a(:,:,:,:) => NULL()
  end type r4p
  type l1p
     logical, pointer :: l(:) => NULL()
  end type l1p 
  type p1p
      type(r2p),    pointer :: a(:) => NULL()
  end type p1p
  type p2p
      type(r3p),    pointer :: a(:) => NULL()
  end type p2p

end module typre
