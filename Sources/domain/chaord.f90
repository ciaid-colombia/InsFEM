subroutine chaord(inoga,ninte)
!-----------------------------------------------------------------------
!
!     Change the ordering of integration points for closed integration
!     rules, according to : INOGA(IGAUS)=INODE
!
!-----------------------------------------------------------------------
  use      typre
  implicit none
  integer(ip), intent(in)  ::   ninte
  integer(ip), intent(out) :: inoga(ninte)

  if(ninte==4) then                                   ! 2 x 2 
     inoga(1)= 1
     inoga(2)= 4
     inoga(3)= 2
     inoga(4)= 3
  else if(ninte==9) then                              ! 3 x 3
     inoga(1)= 1
     inoga(2)= 8
     inoga(3)= 4
     inoga(4)= 5
     inoga(5)= 9
     inoga(6)= 7
     inoga(7)= 2
     inoga(8)= 6
     inoga(9)= 3
  else if(ninte==16) then                             ! 4 x 4
     inoga( 1)= 1
     inoga( 2)=12
     inoga( 3)=11
     inoga( 4)= 4
     inoga( 5)= 5                                        ! 4  10   9   3
     inoga( 6)=13                                        !
     inoga( 7)=16                                        ! 11 16   15  8
     inoga( 8)=10                                        ! 
     inoga( 9)= 6                                        ! 12 13   14  7
     inoga(10)=14                                        ! 
     inoga(11)=15                                        ! 1   5   6   2
     inoga(12)= 9                                        
     inoga(13)= 2
     inoga(14)= 7
     inoga(15)= 8
     inoga(16)= 3
! Añadido
  else if(ninte==25) then                             ! 5 x 5
     inoga( 1)= 1
     inoga( 2)=16
     inoga( 3)=15
     inoga( 4)=14
     inoga( 5)= 4                                        !  4  13  12  11   3
     inoga( 6)= 5                                        !
     inoga( 7)=17                                        ! 14  23  22  21  10
     inoga( 8)=24                                        ! 
     inoga( 9)=23                                        ! 15  24  25  20   9
     inoga(10)=13                                        ! 
     inoga(11)= 6                                        ! 16  17  18  19   8
     inoga(12)=18                                        !                                         
     inoga(13)=25                                        !  1   5   6   7   2
     inoga(14)=22
     inoga(15)=12
     inoga(16)= 7
     inoga(17)=19
     inoga(18)=20
     inoga(19)=21
     inoga(20)=11
     inoga(21)= 2                                        
     inoga(22)= 8                                        
     inoga(23)= 9                                        
     inoga(24)=10                                        
     inoga(25)= 3  
! Hasta aquí 
  else if(ninte==8) then                              ! 2 x 2 x 2 
     inoga(1)= 1
     inoga(2)= 5
     inoga(3)= 4
     inoga(4)= 8
     inoga(5)= 2
     inoga(6)= 6
     inoga(7)= 3
     inoga(8)= 7
  else if(ninte==27) then                             ! 3 x 3 x 3 
     inoga( 1)= 1
     inoga( 2)=13
     inoga( 3)= 5
     inoga( 4)=12
     inoga( 5)=25
     inoga( 6)=20
     inoga( 7)= 4
     inoga( 8)=16
     inoga( 9)= 8
     inoga(10)= 9
     inoga(11)=22
     inoga(12)=17
     inoga(13)=21
     inoga(14)=27
     inoga(15)=26
     inoga(16)=11
     inoga(17)=24
     inoga(18)=19
     inoga(19)= 2
     inoga(20)=14
     inoga(21)= 6
     inoga(22)=10
     inoga(23)=23
     inoga(24)=18
     inoga(25)= 3
     inoga(26)=15
     inoga(27)= 7
  else if(ninte==64) then
     inoga( 1)= 1
     inoga( 2)=12
     inoga( 3)=21
     inoga( 4)= 5
     inoga( 5)=16
     inoga( 6)=44    
     inoga( 7)=52
     inoga( 8)=32
     inoga( 9)=15
     inoga(10)=43
     inoga(11)=51
     inoga(12)=31     
     inoga(13)= 4
     inoga(14)=20
     inoga(15)=24
     inoga(16)= 8
     inoga(17)= 9
     inoga(18)=37      
     inoga(19)=45
     inoga(20)=25
     inoga(21)=33
     inoga(22)=57
     inoga(23)=61
     inoga(24)=53     
     inoga(25)=36
     inoga(26)=60
     inoga(27)=64
     inoga(28)=56
     inoga(29)=14
     inoga(30)=42     
     inoga(31)=50
     inoga(32)=30
     inoga(33)=10
     inoga(34)=38
     inoga(35)=46
     inoga(36)=26     
     inoga(37)=34
     inoga(38)=58
     inoga(39)=62
     inoga(40)=54
     inoga(41)=35
     inoga(42)=59    
     inoga(43)=63
     inoga(44)=55
     inoga(45)=13
     inoga(46)=41
     inoga(47)=49
     inoga(48)=29     
     inoga(49)= 2
     inoga(50)=18
     inoga(51)=22
     inoga(52)= 6
     inoga(53)=11
     inoga(54)=39     
     inoga(55)=47
     inoga(56)=27
     inoga(57)=12
     inoga(58)=40
     inoga(59)=48
     inoga(60)=28     
     inoga(61)= 3
     inoga(62)=19
     inoga(63)=23
     inoga(64)= 7
  end if

end subroutine chaord
