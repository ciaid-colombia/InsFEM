module Mod_Listen
 use typre
 use Mod_iofile
 private
 public ListenFile, maxwp, upcase
!-----------------------------------------------------------------------
!    
!    This common block contains the a%parameters needed for some input-
!    output operations
! 
!-----------------------------------------------------------------------

!
! Listen files
!

 integer(ip), parameter :: maxwp=500

 type ListenFile
   
   integer(ip)            :: nwopa,nnwor,nnpar
   real(rp)               :: param(maxwp)
   character(5)           :: words(maxwp)

   integer(ip)            :: lisda,lisre,lispa
   integer(ip)            :: lisin, lisi1
   character(151)         :: ccard = ' '             ! it must be 80+1
   character(150)         :: DataFolder=' ', fil_include=' '       
   integer(ip)            :: nunit
   logical(lg)            :: echo=.false. ! default echo off. to change it use: echo on
   
contains 
   procedure :: getint
   procedure :: getrea
   procedure :: exists
   procedure :: decod1
   procedure :: SetLunits
   procedure :: getarrs
   procedure :: listen
   procedure :: rewind
   procedure :: listennumbers
   procedure :: SetReadingFolder
   procedure :: CloseIncludeFiles
   
 end type  
 
 
 
 
contains

  subroutine SetLunits(a,lunin,lunout)
   implicit none
   
   class(ListenFile) :: a
   integer(ip), intent(in) :: lunin, lunout
   
   a%lispa = 0
   a%lisda = lunin
   a%lisre = lunout
  end subroutine
  
  subroutine rewind(a)
   implicit none
   class(ListenFile) :: a
   integer(ip) :: stat
   
   rewind(a%lisda, iostat=stat)
  end subroutine
  
  subroutine getarrs(a,pwords,pparam,pnnpar,pnnwor)
   implicit none
   class(ListenFile), target     :: a
   real(rp), pointer     :: pparam(:)
   character(5), pointer :: pwords(:)
   integer(ip), pointer  :: pnnpar,pnnwor
   
   pparam=>a%param
   pwords=>a%words
   pnnpar=>a%nnpar
   pnnwor=>a%nnwor
  end subroutine
  
  subroutine CloseIncludeFiles(a)
   implicit none
   class(ListenFile), target :: a
   
   if (a%nunit == a%lisi1) then
      call CloseFileDom(a%lisi1)
      call CloseFileDom(a%lisin)
   elseif (a%nunit == a%lisin) then
      call CloseFileDom(a%lisin)
   endif   
  end subroutine

  function getint(a,fword,defau,textvari)    
    !-----------------------------------------------------------------------
    !     
    !     This routine gets the integer value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory a%parameter
    !      - texts(1:1).eq.'*' => compulsory a%parameter
    !
    !-----------------------------------------------------------------------
    implicit         none    
    class(ListenFile) :: a
    character(5)  :: fword
    integer       :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    integer       :: value
    integer       :: iword,lengh,first,ndots
    logical       :: found
    character(1)  :: markc
    integer(ip)   :: getint
    !
    ! Formats
    !
1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')
    !
    ! Begin
    !
    value=defau 
    found=.false.
    iword=0
    markc=textvari(1:1)
    
    do while((iword<a%nwopa).and.(.not.found))
       iword=iword+1
       if(a%words(iword)==fword) then
          found=.true.
          value=int(a%param(iword))
       end if
    end do
    
    if(((markc=='!').or.(markc=='*'))&
         .and.(.not.found)) then
       write(a%lisre,1) fword,texts(2:35)
       call runend('GETINT: COMPULSORY PARAM. NOT FOUND')
    end if
    
    getint=value
    
  end function getint
  
  function getrea(a,fword,defau,textvari)
    !-----------------------------------------------------------------------
    !     
    !     This routine gets the real value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory a%parameter
    !      - texts(1:1).eq.'*' => compulsory a%parameter
    !
    !-----------------------------------------------------------------------
    implicit         none
    class(ListenFile) :: a
    character(5)  :: fword
    real(kind=8)  :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    integer       :: iword
    logical       :: found
    real(kind=8)  :: value
    character(1)  :: markc
    real(rp)      :: getrea
    !
    ! Formats
    !
1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')
    !
    ! Begin
    !
    value=defau 
    found=.false.
    iword=0
    markc=textvari(1:1)
    
    do while((iword<a%nwopa).and.(.not.found))
       iword=iword+1
       if(a%words(iword).eq.fword) then
          found=.true.
          value=a%param(iword)
       end if
    end do
    
    if(((markc=='!').or.(markc=='*'))&
         .and.(.not.found)) then
       write(a%lisre,1) fword,texts(2:35)
       call runend('GETREA: COMPULSORY PARAM. NOT FOUND')
    end if
    
    getrea=value
    
  end function getrea
  
  function exists(a,fword)
    !-----------------------------------------------------------------------
    !     
    !     This routine verifys if fword exists in a%words
    !
    !-----------------------------------------------------------------------
    implicit none
    class(ListenFile) :: a
    character(5) :: fword
    integer      :: iword
    logical(lg)  :: exists

    exists=.false.  
    iword=0
    do while((iword<a%nwopa).and.(.not.exists))
       iword=iword+1
       if(a%words(iword)==fword) then
          exists=.true.
       end if
    end do
    
  end function exists
  
  subroutine upcase(word)
    !-----------------------------------------------------------------------
    !
    !     This routine converts word to upper case 
    !
    !-----------------------------------------------------------------------
    implicit none
    character(5), intent(inout) :: word
    integer(ip)                 :: iposi,ioctv
    
    do iposi=1,5                                   ! process all positions
       ioctv=ichar(word(iposi:iposi))               ! octal value
       if(o'141'<=ioctv.and.o'172'>=ioctv) then ! it is a lower case
          ioctv=ioctv-o'40'                          ! equivalent upper case
          word(iposi:iposi)=char(ioctv)              ! convert it to upcase
       end if
    end do ! iposi=1,5
    
  end subroutine upcase

  subroutine decod1(a,nstri,forma,strin,lflag,digit)
    !-----------------------------------------------------------------------
    !
    !     This routine decodified a string
    !
    !-----------------------------------------------------------------------
    implicit none
    class(ListenFile) :: a
    integer(ip)    :: nstri,lflag,istri,decim
    real(kind=8)   :: digit
    character(7)   :: forma
    character(*)   :: strin                           ! (nstri+1)
    character(151) :: stri1=' '
    integer(ip)    :: stat
    
    lflag=0                                           ! It is a a%parameter
    istri=1                                           
    do while(istri<=nstri)                          
       decim=ichar(strin(istri:istri))                ! decimal value
       if (decim<48.or.decim>57) then                 ! It is not a num.
          if (&
               decim/= 43.and.decim/= 45.and. &       ! discard + -
               decim/= 68.and.decim/= 69.and. &       ! discard D E
               decim/=100.and.decim/=101.and. &       ! discard d e
               decim/= 46) then                       ! discard .
             lflag=1
             istri=nstri
          end if
       end if
       istri=istri+1
    end do
    
    if (lflag==0) then                                ! It's a number
       do istri=1,nstri
          stri1(istri:istri)=strin(istri:istri)
       end do
       stri1(nstri+1:nstri+1)=' '
       stri1=adjustl(trim(stri1))
       read(stri1,*,IOSTAT=stat) digit                             ! DIGIT<-'STRIN'
       
    end if
    
  end subroutine decod1
 
  subroutine listen(a,subna)

!-----------------------------------------------------------------------
!
! Reads a string and interprets it as a%words and a%parameters.
!
!   - Maximum number of a%words and a%parameters = maxwp.
!   - Only the first five characters of each word are decoded.
!   - The underline characters '_' are discarted.
!   - Lower case letters are converted to upper case.
!   - Each word or a%parameter must be separated by ' ', '=', ':' or ','
!   - A "comment" begins with '$', '!', '/'  or '/' in any place.
!   - A line that has a comment beginning with '/' or '/'
!     continuates in the next line.
!   - A line beginning with title is not decoded. title directive.
!   - A line beginning with include is an include directive.
!   - A line beginning with echo turns 'on' or 'off' the echo.
!
!-----------------------------------------------------------------------

  implicit none
!
! Var
!
  class(ListenFile) :: a
  character(6)      :: subna
  
  real(kind=8)      :: digit
  integer(ip)       :: first,firsp,i,last,lastp,ptrwo,npptr,nwptr
  integer(ip)       :: leng,flag,resum,istat
  logical(lg)       :: newline
  
!
! Data.
!      
  if(a%lispa==0) then
     a%nunit = a%lisda                                  !  initial data file
     a%lispa = 1
  end if
!
! Format
!
10 format(a150)
20 format(1x,a6,' <-- ',a)
!
! Begin
!
  a%ccard=' '
  a%nnwor=0                                           ! initialize.
  a%nnpar=0
  nwptr=0
  npptr=0
  resum=0
  do i=1,maxwp
     a%words(i)=' '
     a%param(i)=0.0_rp
  end do
100 continue
  do while(((a%nnwor==0).and.(a%nnpar==0)&              ! don't return without answer
       .or.newline).or.resum==1)                    ! continue reading if / or \
     
     if (resum==0) then
        newline=.false.                             ! initialize.
        last=0
        lastp=0
     end if
     resum=0
     
     read(a%nunit,10,end=101,err=1,iostat=istat) a%ccard             ! read a card
!     leng=lnbln1(a%ccard)                             ! calculate the length.
     leng=len_trim(a%ccard)                             ! calculate the length.
     
     decode_card: do while(last<leng)               ! decode all the card.
        first=last+1
        do while(first<=150        .and.(&
             first<leng            .and. &              
             a%ccard(first:first)=='_'.or. &          ! jump null character (_)
             a%ccard(first:first)==' '.or. &          ! jump separators ( =:,)
             a%ccard(first:first)=='='.or. &
             a%ccard(first:first)==':'.or. &
             a%ccard(first:first)==','    ))
           first=first+1
        end do
        if(last==0) firsp=first                     ! save first to print card
        last=first
        do while(last<=151        .and.(&
             (last<=leng)         .and. &           ! look for end (last=leng).
             a%ccard(last:last)/=' '.and. &           ! look for separator ( =:,).
             a%ccard(last:last)/='='.and. &
             a%ccard(last:last)/=':'.and. &
             a%ccard(last:last)/=','.and. &
             a%ccard(last:last)/='$'.and. &           ! look for coment ($!/\).
             a%ccard(last:last)/='!'.and. &
             a%ccard(last:last)/='\'.and. &
             a%ccard(last:last)/='\'     ))
           last=last+1
        end do
        if(last<=151            .and.(&
             a%ccard(last:last)=='$'.or.&
             a%ccard(last:last)=='!')) leng=last-1
        if(last<=151            .and.(&
             a%ccard(last:last)=='\'.or.&
             a%ccard(last:last)=='\')) then
           leng=last-1                              ! deal with continuation
           newline=.true.                           ! set new line flag.
        end if
        last=last-1
        if(last>=first) then                        ! is it a word or a a%parameter
           lastp=last                               ! save last to print card
           call a%decod1(last-first+1_ip,'(f20.0)',&
                a%ccard(first:last),flag,digit)
           if(flag==0) then                         ! it is a a%parameter.
              a%nnpar=a%nnpar+1                         ! # of a%parameters
              npptr=npptr+1                         ! integer :: to next a%parameter
              if(npptr>maxwp) go to 4               ! error.
              if(nwptr>npptr) npptr=nwptr
              a%param(npptr)=digit
           else                                     ! it is a word.
              a%nnwor=a%nnwor+1                         ! # of a%words
              nwptr=nwptr+1                         ! integer :: to next word
              if(nwptr>maxwp) go to 5               ! error.
              if(npptr>=nwptr) nwptr=npptr+1
              ptrwo=1
              do while ((first<=last).and.(ptrwo<=5))
                 a%words(nwptr)(ptrwo:ptrwo)=a%ccard(first:first)
                 ptrwo=ptrwo+1
                 first=first+1
                 do while (a%ccard(first:first)=='_') ! jump null character
                    first=first+1
                 end do
              end do
              call upcase(a%words(nwptr)) ! convert to upper case.
           end if
        end if ! (last>=first)
        if((a%nnwor==1).and.(a%nnpar==0).and.&          ! deal with title or include
             ((a%words(1)=='TITLE').or.((a%words(1)=='INCLU').and.&
             (subna/='NOREAD')))) then
           if(a%echo.and.(subna/='noecho'))&
                write(a%lisre,20) 'listen',a%ccard(firsp:leng)
           last=last+2
           do while(a%ccard(last:last)==' ')          ! remove blank spaces
              last=last+1
           end do
           if(leng<last) go to 6                    ! error
           a%ccard=a%ccard(last:leng)                   ! remove a%words(1) from a%ccard
           leng=leng-last+1
           if(a%words(1)=='TITLE') then               ! deal with titles
              firsp=1
              lastp=leng
           else                                     ! deal with include directive
              !if(a%nunit==a%lisin) go to 3              ! error
              last=1                                ! remove tail comments
              do while((last<=leng).and.&           ! look for end (last=leng)
                 a%ccard(last:last)/=' ')           ! look for separator ( )
                 last=last+1
              end do
              a%ccard=a%ccard(1:last-1)                 ! remove tail comments
              a%fil_include = adjustl(trim(a%DataFolder))//'/'//a%ccard
              if(a%nunit==a%lisin) then
                 call iofile(0,a%nunit,a%fil_include,'LISTEN_INCLUDE','old','formatted')
                 a%lisi1 = a%nunit
              else
                 call iofile(0,a%nunit,a%fil_include,'LISTEN_INCLUDE','old','formatted')
                 a%lisin = a%nunit
              end if
              !open(unit=a%nunit,file=a%ccard,&
              !     form='formatted',status='old')
              lastp=0                               ! to ignore the echo
              a%nnwor=0                               ! forget all
              nwptr=0
              a%words(1)=' '
           end if
           last=leng                                ! to break the do while
        else if(subna=='NOREAD') then
           newline=.false.
           a%nnwor=1
        end if
     end do decode_card 
     
     if((a%words(1)=='ECHO') .and.&
          (a%nnpar==0).and.(a%nnwor==2)) then           ! deal with echo
        if(a%words(2)=='OFF') then
           a%echo=.false.
           if(subna/='noecho') write(a%lisre,20) 'LISTEN','ECHO OFF'
        else
           a%echo=.true.
           if(subna/='noecho') write(a%lisre,20) 'LISTEN','ECHO ON'
        endif
        a%nnwor=0                                     ! forget all
        nwptr=0
        do i=1,maxwp
           a%words(i)=' '
        end do
     else                                           ! print card
        if((a%echo).and.(firsp<=lastp).and.(subna/='noecho')) then
           if(newline) then
              lastp=lastp+2
              a%ccard(lastp-1:lastp)=' /'
           end if
           write(a%lisre,20) subna,a%ccard(firsp:lastp)
        end if
     end if
     
  end do ! while ((a%nnwor==0).and.(a%nnpar==0).or.newline)
  
  a%nwopa=max(npptr,nwptr) 
  return
!
! End of include
!
101 if(a%nunit/=a%lisin.and.a%nunit/=a%lisi1) go to 2 ! error
  close(unit=a%nunit,status='keep')
  if(a%nunit==a%lisi1) then
     a%nunit=a%lisin
  else
     a%nunit=a%lisda
  end if
  if(a%echo.and.(subna/='noecho'))&
       write(a%lisre,20) 'LISTEN','END OF INCLUDE FILE'
  resum=1
  go to 100 ! for resume the error return to the same place.
!
! Errors:
!
1 call runend('LISTEN: ERROR DETECTED WHEN READING')
2 call runend('LISTEN: END OF FILE DETECTED       ')
3 call runend('LISTEN: ERROR: INCLUDE FROM INCLUDE')
4 call runend('LISTEN: TOO MANY PARAM IN COMMAND  ')
5 call runend('LISTEN: TOO MANY WORDS IN COMMAND  ')
6 call runend('LISTEN: BLANK IS ILEGAL HEAR       ')
  
end subroutine listen

subroutine listennumbers(a,nnumbers,subna)
   implicit none
   class(ListenFile) :: a
   character(6)      :: subna
   integer(ip)       :: nnumbers
   
   integer(ip)       :: myIOStat1,auxcount
   real(rp), parameter :: nan = 1.458968237e12
   
   a%nnpar = nnumbers
   read(a%nunit,*) a%param(1:nnumbers)
end subroutine

subroutine SetReadingFolder(a,DataFolder)
   implicit none
   class(ListenFile) :: a
   character(150)      :: DataFolder

   a%DataFolder = Datafolder
end subroutine
  
!  function lnbln1(a%ccard)
!    !-----------------------------------------------------------------------
!    !
!    !     Compute the length of a%ccard
!    !
!    !-----------------------------------------------------------------------
!    implicit none
!    character(1) :: a%ccard(151)
!    integer(ip)  :: lnbln1
!    
!    lnbln1=151
!    do while (a%ccard(lnbln1)==' '.and.lnbln1/=1)
!       lnbln1=lnbln1-1
!    end do
!    
!    if (lnbln1==1.and.a%ccard(1)==' ') lnbln1=0
!    
!  end function lnbln1
  
  
end module Mod_Listen

