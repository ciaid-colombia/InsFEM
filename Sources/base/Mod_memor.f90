module Mod_memor
!-----------------------------------------------------------------------
!
!    This module contains memory allocation class
!
!-----------------------------------------------------------------------
  use typre
  implicit none
  private
  public :: MemoryMan, TotalMemo, TotalMax!, alloc
  
  integer(8) :: TotalMemo = 0, TotalMax
  
  type MemoryMan 
      private
      !public :: alloc
      integer(8) :: MaxMemo = 0
      integer(8) :: CurrentMemo = 0
      integer(ip) :: lumem = 0
      integer(ip) :: luerr = 0  
      
contains
     !private
     !public :: alloc
     procedure :: allocrp1
     procedure :: allocrp2
     procedure :: allocrp3
     procedure :: allocrp4
     procedure :: allocip1
     procedure :: allocip2
     procedure :: allocip3          ! integers
     procedure :: pallocip1
     procedure :: pallocip2
     procedure :: pallocip3          ! integers pointers
     procedure :: pallocip4
     procedure :: alloci1p1                        ! i1p types
     procedure :: alloci2p1                        ! i2p types
     procedure :: allocr1p1
     procedure :: allocr1p2                 ! r1p types
     procedure :: allocr2p1
     procedure :: allocr2p2                 ! r2p types
     procedure :: allocr3p1
     procedure :: allocr4p1
     procedure :: allocr3p2                   ! r3p types 
     procedure :: pallocrp1
     procedure :: pallocrp2
     procedure :: pallocrp3
     procedure :: pallocrp4
     procedure :: alloclg1
     procedure :: alloclg2
     procedure :: allocch1
!      procedure :: allocPostprFile
     
     generic :: alloc => allocrp1, allocrp2,allocrp3,allocrp4, &
            allocip1,allocip2,allocip3,&          ! integers			  
            alloci1p1,&                        ! i1p types
            alloci2p1,&                        ! i2p types
            allocr1p1,allocr1p2,&                 ! r1p types
            allocr2p1,allocr2p2,&                 ! r2p types
            allocr3p1,allocr3p2,&                   ! r3p types 
            alloclg1,alloclg2,& ! logicals
            allocch1,allocr4p1 ! characters

     procedure :: allocObj 

     generic :: palloc => pallocip1,pallocip2,pallocip3,pallocip4,&        ! integers
                          pallocrp1,pallocrp2,pallocrp3,pallocrp4,&
                          palloclp1
     
     procedure :: deallocrp1
     procedure :: deallocrp2
     procedure :: deallocrp3
     procedure :: deallocrp4
     procedure :: deallocip1
     procedure :: deallocip2
     procedure :: deallocip3          ! integers
     procedure :: pdeallocip1
     procedure :: pdeallocip2
     procedure :: pdeallocip3          ! integers
     procedure :: pdeallocip4
     procedure :: dealloci1p1                        ! i1p types
     procedure :: dealloci2p1                        ! i2p types
     procedure :: deallocr1p1
     procedure :: deallocr1p2                 ! r1p types
     procedure :: deallocr2p1
     procedure :: deallocr2p2                 ! r2p types
     procedure :: deallocr3p1
     procedure :: deallocr4p1
     
     procedure :: deallocr3p2                   ! r3p types 
     procedure :: pdeallocrp1
     procedure :: pdeallocrp2
     procedure :: pdeallocrp3
     procedure :: pdeallocrp4
     
     procedure :: palloclp1
     procedure :: pdealloclp1
     
     procedure :: dealloclg1
     procedure :: dealloclg2
     
     procedure :: deallocch1
     
     generic :: dealloc => deallocrp1, deallocrp2,deallocrp3,deallocrp4, &
         deallocip1,deallocip2,deallocip3,&          ! integers
         dealloci1p1,&                        ! i1p types
         dealloci2p1,&                        ! i2p types
         deallocr1p1,deallocr1p2,&            ! r1p types
         deallocr2p1,deallocr2p2,&            ! r2p types
         deallocr3p1,deallocr3p2,&            ! r3p types 
         dealloclg1,dealloclg2,&              ! logicals
         deallocch1,deallocr4p1                           ! characters
     
     generic :: pdealloc => pdeallocip1,pdeallocip2,pdeallocip3,pdeallocip4, &
                            pdeallocrp1,pdeallocrp2,pdeallocrp3,pdeallocrp4,&
                            pdealloclp1
     
     procedure :: deallocObj
     
     !Realloc
     procedure :: reallocip1
     procedure :: reallocrp1
     procedure :: reallocrp2
     procedure :: reallocip2
     procedure :: reallocip3
     procedure :: reallocrp3
     procedure :: reallocrp4
     
     generic :: realloc => reallocip1,reallocrp1,reallocrp2,reallocip2,reallocip3,reallocrp3,reallocrp4
     
     procedure :: preallocip1
     procedure :: preallocrp1
     procedure :: preallocrp2
     generic :: prealloc => preallocip1, preallocrp1,preallocrp2
     
     
     procedure :: init
     procedure :: GetValue
     procedure :: GetLumem
  
  end type


contains


  subroutine init(a,lumem,luerr)
    implicit none
    class(MemoryMan) :: a
    integer(ip), intent(in) :: lumem, luerr
      
    a%MaxMemo = 0
    a%CurrentMemo = 0
    a%lumem = lumem
    a%luerr = luerr
  end subroutine
  
  subroutine GetValue(a,Current,MaxMemo,TotalMem,TotalMa)
    implicit none
    class(MemoryMan) :: a
    integer(8) :: Current, Maxmemo, TotalMem,TotalMa
      
    MaxMemo = a%MaxMemo
    Current = a%CurrentMemo
    TotalMem = TotalMemo
    TotalMa  = TotalMax
  end subroutine
  
  subroutine memctr(CurrentMemo,MaxMemo,lbyts,vanam,lumem,vacal)
      use typre
      implicit none
      integer(ip), intent(in)    :: lbyts
      character(*), intent(in)   :: vanam
      integer(ip), intent(in)    :: lumem
      character(*), intent(in)   :: vacal       ! Calling routine
      integer(8), intent(inout) :: MaxMemo,CurrentMemo
      character(6)               :: lbyte, lbyte2, lbyte3
      real(rp) :: rMaxMemo, rCurrentMemo,rTotalMemo,rMemory,rdummy

      CurrentMemo = CurrentMemo + lbyts
      TotalMemo   = TotalMemo   + lbyts
      MaxMemo = max(MaxMemo,CurrentMemo)
      TotalMax = max(TotalMemo,TotalMax)

      if(lumem>0) then
         call Mem2String(CurrentMemo,rCurrentMemo,rdummy,lbyte)
         call Mem2String(TotalMemo,rTotalMemo,rdummy,lbyte2)
         call Mem2String(int(lbyts,8),rMemory,rdummy,lbyte3)
         
         write(lumem,'(a12,2x,a8,f8.1,1x,a6,2x,a13,f8.1,1x,a6,2x,a12,2x,f8.1,1x,a6,1x,a10)') &
            trim(vanam),'Memory: ', rMemory,lbyte3,'Module memo: ', rCurrentMemo,lbyte, 'Total memo: ', rTotalMemo,lbyte2,vacal
      else
         !write(*,*) 'hey check this: betterh than runend, some buffers do not have output channel'
      
      endif

end subroutine memctr
  
  
!-------------------------------------------------------------------
! General Objects (unknown but fixed size)
!-------------------------------------------------------------------

  subroutine allocObj(a,istat,vanam,vacal,lbyts)
    implicit none
    class(MemoryMan)            :: a
    integer(ip)                 :: istat
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    integer(ip), intent(in)     :: lbyts
    
    if(istat==0) then
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  
  
  end subroutine
  
  subroutine deallocObj(a,istat,vanam,vacal,elbyts)
    implicit none
    class(MemoryMan)            :: a
    integer(ip)                 :: istat
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    integer(ip), intent(in)     :: elbyts
    
    integer(ip) :: lbyts
    
    if(istat==0) then
      lbyts=-elbyts
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
!-----------------------------------------------------------------------
! Real matrices
!---------------------------------------------------------------------
 
  subroutine allocrp1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    real(rp), allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      varia=0.0_rp
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocrp1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    real(rp), allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine reallocrp1(a,newndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1
    
    real(rp), allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    real(rp), allocatable :: aux_varia(:)
    
    integer(ip) :: oldndim1
    
    allocate(aux_varia(newndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*newndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia)
    
    aux_varia(1:min(newndim1,oldndim1)) = varia(1:min(newndim1,oldndim1))
    call move_alloc(aux_varia,varia)
    
    lbyts = -rp*oldndim1
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine  
    
  subroutine allocrp2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    real(rp), allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      varia=0.0_rp
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocrp2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    real(rp), allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine reallocrp2(a,newndim1,newndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1,newndim2
    
    real(rp), allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    real(rp), allocatable :: aux_varia(:,:)
    
    integer(ip) :: oldndim1,oldndim2
    
    allocate(aux_varia(newndim1,newndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*newndim1*newndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia,1)
    oldndim2 = size(varia,2)
    
    lbyts = -rp*size(varia)
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    aux_varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2)) = varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2))
    call move_alloc(aux_varia,varia)
    
    
  end subroutine    
    
  subroutine allocrp3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    real(rp), allocatable           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2*ndim3
      varia=0.0_rp
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocrp3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    real(rp), allocatable           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2*ndim3
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine reallocrp3(a,newndim1,newndim2,newndim3,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1,newndim2,newndim3
    
    real(rp), allocatable           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    real(rp), allocatable :: aux_varia(:,:,:)
    
    integer(ip) :: oldndim1,oldndim2,oldndim3
    
    allocate(aux_varia(newndim1,newndim2,newndim3),stat=istat)
    if(istat==0) then
      lbyts=rp*newndim1*newndim2*newndim3
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia,1)
    oldndim2 = size(varia,2)
    oldndim3 = size(varia,3)
    
    lbyts = -rp*oldndim1*oldndim2*oldndim3
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    aux_varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2),1:min(newndim3,oldndim3)) = varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2),1:min(newndim3,oldndim3))
    call move_alloc(aux_varia,varia)
    
    
  end subroutine    
  
  subroutine allocrp4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    real(rp), allocatable           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3,ndim4),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2*ndim3*ndim4
      varia=0.0_rp
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocrp4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    real(rp), allocatable           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2*ndim3*ndim4
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine reallocrp4(a,newndim1,newndim2,newndim3,newndim4,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1,newndim2,newndim3,newndim4
    
    real(rp), allocatable           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    real(rp), allocatable :: aux_varia(:,:,:,:)
    
    integer(ip) :: oldndim1,oldndim2,oldndim3,oldndim4
    
    allocate(aux_varia(newndim1,newndim2,newndim3,newndim4),stat=istat)
    if(istat==0) then
      lbyts=rp*newndim1*newndim2*newndim3*newndim4
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia,1)
    oldndim2 = size(varia,2)
    oldndim3 = size(varia,3)
    oldndim4 = size(varia,4)
    
    lbyts = -rp*oldndim1*oldndim2*oldndim3*oldndim4
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    aux_varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2),1:min(newndim3,oldndim3),1:min(newndim4,oldndim4)) = varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2),1:min(newndim3,oldndim3),1:min(newndim4,oldndim4))
    call move_alloc(aux_varia,varia)
    
    
  end subroutine      
    
  
!-----------------------------------------------------------------------
! Integer Matrices
!-----------------------------------------------------------------------
  
  subroutine allocip1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    integer(ip), allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocip1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    integer(ip), allocatable    :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine reallocip1(a,newndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1
    
    integer(ip), allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    integer(ip), allocatable :: aux_varia(:)
    
    integer(ip) :: oldndim1
    
    allocate(aux_varia(newndim1),stat=istat)
    if(istat==0) then
      lbyts=ip*newndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia)
    
    aux_varia(1:min(newndim1,oldndim1)) = varia(1:min(newndim1,oldndim1))
    call move_alloc(aux_varia,varia)
    
    lbyts = -ip*oldndim1
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine  
  
  subroutine allocip2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    integer(ip), allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1*ndim2
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocip2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    integer(ip), allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
    
  subroutine reallocip2(a,newndim1,newndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1,newndim2
    
    integer(ip), allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    integer(ip), allocatable :: aux_varia(:,:)
    
    integer(ip) :: oldndim1,oldndim2
    
    allocate(aux_varia(newndim1,newndim2),stat=istat)
    if(istat==0) then
      lbyts=ip*newndim1*newndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia,1)
    oldndim2 = size(varia,2)
    
    aux_varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2)) = varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2))
    call move_alloc(aux_varia,varia)
    
    lbyts = -ip*oldndim1*oldndim2
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine     
    
  subroutine allocip3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    integer(ip), allocatable           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1*ndim2*ndim3
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocip3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    integer(ip), allocatable           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1*ndim2*ndim3
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine reallocip3(a,newndim1,newndim2,newndim3,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1,newndim2,newndim3
    
    integer(ip), allocatable           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    integer(ip), allocatable :: aux_varia(:,:,:)
    
    integer(ip) :: oldndim1,oldndim2,oldndim3
    
    allocate(aux_varia(newndim1,newndim2,newndim3),stat=istat)
    if(istat==0) then
      lbyts=ip*newndim1*newndim2*newndim3
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia,1)
    oldndim2 = size(varia,2)
    oldndim3 = size(varia,3)
    
    aux_varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2),1:min(newndim3,oldndim3)) = varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2),1:min(newndim3,oldndim3))
    call move_alloc(aux_varia,varia)
    
    lbyts = -ip*oldndim1*oldndim2*oldndim3
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine   
  
  
  
  
  subroutine allocip4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    integer(ip), allocatable           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3,ndim4),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1*ndim2*ndim3*ndim4
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocip4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    integer(ip), allocatable           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1*ndim2*ndim3*ndim4
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
!----------------------------------------------------------------------
! Integer Pointers
!----------------------------------------------------------------------
  subroutine pallocip1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    integer(ip), pointer        :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocip1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    integer(ip), pointer        :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine preallocip1(a,newndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1
    
    integer(ip), pointer        :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    integer(ip), allocatable, target :: aux_varia(:)
    
    integer(ip) :: oldndim1
    
    allocate(aux_varia(newndim1),stat=istat)
    if(istat==0) then
      lbyts=ip*newndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia)
    
    aux_varia(1:min(newndim1,oldndim1)) = varia(1:min(newndim1,oldndim1))
    deallocate(varia)
    allocate(varia(newndim1))
    varia = aux_varia
    deallocate(aux_varia)
    
    lbyts = -ip*oldndim1
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine  
  
  subroutine pallocip2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    integer(ip), pointer        :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1*ndim2
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocip2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    integer(ip), pointer        :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
    
  subroutine pallocip3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    integer(ip), pointer        :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1*ndim2*ndim3
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocip3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    integer(ip), pointer        :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1*ndim2*ndim3
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine pallocip4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    integer(ip), pointer        :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3,ndim4),stat=istat)
    if(istat==0) then
      lbyts=ip*ndim1*ndim2*ndim3*ndim4
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocip4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    integer(ip), pointer        :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ip*ndim1*ndim2*ndim3*ndim4
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
!----------------------------------------------------------------------
! Real Pointers
!----------------------------------------------------------------------
  subroutine pallocrp1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    real(rp), pointer           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocrp1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    real(rp), pointer           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine preallocrp1(a,newndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1
    
    real(rp), pointer           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    real(rp), allocatable, target :: aux_varia(:)
    
    integer(ip) :: oldndim1
    
    allocate(aux_varia(newndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*newndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia)
    
    aux_varia(1:min(newndim1,oldndim1)) = varia(1:min(newndim1,oldndim1))
    deallocate(varia)
    allocate(varia(newndim1))
    varia = aux_varia
    deallocate(aux_varia)
    
    lbyts = -rp*oldndim1
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine  
  
  subroutine pallocrp2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    real(rp), pointer           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocrp2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    real(rp), pointer           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
    
  subroutine preallocrp2(a,newndim1,newndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip), intent(in)     :: newndim1,newndim2
    
    real(rp), pointer           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    real(rp), pointer :: aux_varia(:,:) => NULL()
    
    integer(ip) :: oldndim1,oldndim2
    
    allocate(aux_varia(newndim1,newndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*newndim1*newndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    oldndim1 = size(varia,1)
    oldndim2 = size(varia,2)
    
    lbyts = -rp*size(varia)
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
    
    aux_varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2)) = varia(1:min(newndim1,oldndim1),1:min(newndim2,oldndim2))
    deallocate(varia)
    varia => aux_varia
    
    
  end subroutine    
    
  subroutine pallocrp3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    real(rp), pointer        :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2*ndim3
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocrp3(a,ndim1,ndim2,ndim3,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3       ! Size of the variable
    real(rp), pointer           :: varia(:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2*ndim3
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine pallocrp4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    real(rp), pointer           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2,ndim3,ndim4),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2*ndim3*ndim4
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdeallocrp4(a,ndim1,ndim2,ndim3,ndim4,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2,ndim3,ndim4       ! Size of the variable
    real(rp), pointer           :: varia(:,:,:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2*ndim3*ndim4
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  
  
!-----------------------------------------------------------------------
! ip matrices
!---------------------------------------------------------------------
 
  subroutine alloci1p1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(i1p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      do idime = 1,ndim1
         varia(idime)%l => NULL()
      enddo   
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine dealloci1p1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(i1p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine alloci1p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(i1p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime,jdime

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      do idime = 1,ndim2
         do jdime = 1,ndim1
            varia(jdime,idime)%l => NULL()
         enddo
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine dealloci1p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(i1p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
    
 
!-----------------------------------------------------------------------
! i2p matrices
!---------------------------------------------------------------------
 
  subroutine alloci2p1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(i2p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      do idime = 1,ndim1
         varia(idime)%l => NULL()
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine dealloci2p1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(i2p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine alloci2p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(i2p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime,jdime

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      do idime = 1,ndim2
         do jdime = 1,ndim1
            varia(jdime,idime)%l => NULL()
         enddo
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine dealloci2p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(i2p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  
!-----------------------------------------------------------------------
! r1p matrices
!---------------------------------------------------------------------
 
  subroutine allocr1p1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r1p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      do idime = 1,ndim1
         varia(idime)%a => NULL()
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr1p1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r1p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine allocr1p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(r1p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime,jdime

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      do idime = 1,ndim2
         do jdime = 1,ndim1
            varia(jdime,idime)%a => NULL()
         enddo
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr1p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(r1p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
    

!-----------------------------------------------------------------------
! r2p matrices
!---------------------------------------------------------------------
 
  subroutine allocr2p1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r2p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,jdime

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      do jdime = 1,ndim1
         varia(jdime)%a => NULL()
      enddo
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr2p1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r2p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine allocr2p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(r2p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,idime,jdime

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      do idime = 1,ndim2
         do jdime = 1,ndim1
            varia(jdime,idime)%a => NULL()
         enddo
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr2p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(r2p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
    

!-----------------------------------------------------------------------
! r3p matrices
!---------------------------------------------------------------------
 
  subroutine allocr3p1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r3p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,jdime

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      do jdime = 1,ndim1
         varia(jdime)%a => NULL()
      enddo
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr3p1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r3p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine allocr4p1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r4p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,jdime

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1
      do jdime = 1,ndim1
         varia(jdime)%a => NULL()
      enddo
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr4p1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    type(r4p), allocatable      :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine allocr3p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(r3p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts,jdime,idime

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=rp*ndim1*ndim2
      do idime = 1,ndim2
         do jdime = 1,ndim1
            varia(jdime,idime)%a => NULL()
         enddo
      enddo 
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocr3p2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    type(r3p), allocatable      :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-rp*ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
!----------------------------------------------------------------------  
! logicals
!----------------------------------------------------------------------
  
  subroutine alloclg1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    logical, allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=ndim1
      varia=.false.
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine dealloclg1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    logical, allocatable           :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  subroutine alloclg2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    logical, allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1,ndim2),stat=istat)
    if(istat==0) then
      lbyts=ndim1*ndim2
      varia=.false.
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine dealloclg2(a,ndim1,ndim2,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1,ndim2       ! Size of the variable
    logical, allocatable           :: varia(:,:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ndim1*ndim2
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  !Logical pointers
    subroutine palloclp1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    logical, pointer            :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=1*ndim1
      varia=0_ip
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine pdealloclp1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    logical, pointer            :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-1*ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
  
!----------------------------------------------------------------------  
! characters
!----------------------------------------------------------------------
  
  subroutine allocch1(a,ndim1,varia,vanam,vacal)
    implicit none
    
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    character(*), allocatable :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts

    allocate(varia(ndim1),stat=istat)
    if(istat==0) then
      lbyts=ndim1
      varia=''
    else
      call memerr(0,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine

  subroutine deallocch1(a,ndim1,varia,vanam,vacal)
    implicit none
    class(MemoryMan)            :: a
    integer(ip),  intent(in)    :: ndim1       ! Size of the variable
    character(*), allocatable :: varia(:)
    character(*), intent(in)    :: vanam       ! Variable name
    character(*), intent(in)    :: vacal       ! Calling routine
    
    integer(ip)                 :: istat,lbyts
    
    deallocate(varia,stat=istat)
    if(istat==0) then
      lbyts=-ndim1
    else
      call memerr(1,vanam,vacal,a%luerr)
    end if
    call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
  end subroutine
  
! !---------------------------------------------------------------------------
! ! ListenFile
! !---------------------------------------------------------------------------
! 
!   subroutine allocListenFile(a,varia,vanam,vacal)
!     use Mod_Listen
!     implicit none
!     
!     class(MemoryMan)            :: a
!     type(ListenFile), allocatable           :: varia
!     character(*), intent(in)    :: vanam       ! Variable name
!     character(*), intent(in)    :: vacal       ! Calling routine
!     
!     integer(ip)                 :: istat,lbyts
! 
!     allocate(varia,stat=istat)
!     if(istat==0) then
!       lbyts=maxwp*(5+5*rp)+8*ip+151
!     else
!       call memerr(0,vanam,vacal,a%luerr)
!     end if
!     call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
!   end subroutine
! 
!   subroutine deallocListenFile(a,varia,vanam,vacal)
!     use Mod_Listen
!     implicit none
!     class(MemoryMan)            :: a
!     type(ListenFile), allocatable           :: varia(:)
!     character(*), intent(in)    :: vanam       ! Variable name
!     character(*), intent(in)    :: vacal       ! Calling routine
!     
!     integer(ip)                 :: istat,lbyts
!     
!     deallocate(varia,stat=istat)
!     if(istat==0) then
!       lbyts=-maxwp*(5+5*rp)+8*ip+151
!     else
!       call memerr(1,vanam,vacal,a%luerr)
!     end if
!     call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
!   end subroutine

!---------------------------------------------------------------------------
! PostprFile
!---------------------------------------------------------------------------

!   subroutine allocPostprFile(a,varia,vanam,vacal)
!     use Mod_postpr
!     implicit none
!     
!     class(MemoryMan)            :: a
!     type(PostprFile), pointer           :: varia
!     character(*), intent(in)    :: vanam       ! Variable name
!     character(*), intent(in)    :: vacal       ! Calling routine
!     
!     integer(ip)                 :: istat,lbyts
! 
!     allocate(varia,stat=istat)
!     if(istat==0) then
!       lbyts=6*ip
!     else
!       call memerr(0,vanam,vacal,a%luerr)
!     end if
!     call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
!   end subroutine
! 
!   subroutine deallocPostprFile(a,varia,vanam,vacal)
!     use Mod_postpr
!     implicit none
!     class(MemoryMan)            :: a
!     type(PostprFile), allocatable           :: varia(:)
!     character(*), intent(in)    :: vanam       ! Variable name
!     character(*), intent(in)    :: vacal       ! Calling routine
!     
!     integer(ip)                 :: istat,lbyts
!     
!     deallocate(varia,stat=istat)
!     if(istat==0) then
!       lbyts=-6*ip
!     else
!       call memerr(1,vanam,vacal,a%luerr)
!     end if
!     call memctr(a%CurrentMemo,a%MaxMemo,lbyts,vanam,a%lumem,vacal)
!   end subroutine

  
!-----------------------------------------------------------------------
! Count and error functions
!-----------------------------------------------------------------------
   

   subroutine memerr(itask,vanam,vacal,luerr)
      implicit none
      integer(ip), intent(in) :: itask,luerr
      character*(*)           :: vanam,vacal
      if(itask==0) then
         if(luerr>0) then
            write(luerr,'(a)') trim(vacal)//': MEMORY FOR '//trim(vanam)//' CANNOT BE ALLOCATED'
         else
            write(6,'(a)') trim(vacal)//': MEMORY FOR '//trim(vanam)//' CANNOT BE ALLOCATED'
         end if
         call runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' CANNOT BE ALLOCATED')
         stop
      else
         if(luerr>0) then
            write(luerr,'(a)') trim(vacal)//': MEMORY FOR '//trim(vanam)//' CANNOT BE DEALLOCATED'
         else
            write(6,'(a)') trim(vacal)//': MEMORY FOR '//trim(vanam)//' CANNOT BE DEALLOCATED'
         end if
         call runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' CANNOT BE DEALLOCATED')
         stop
      end if
   end subroutine memerr
   
   subroutine GetLumem(a,lumem)
      implicit none
      class(MemoryMan) :: a
      integer(ip) :: lumem
      
      lumem = a%lumem
   end subroutine
      

end module Mod_memor

