subroutine sortb(n,a)
use      typre
implicit none
integer(ip)   n,a(n),i,j,t
do i=1,n-1
   do j=i+1,n
      if (a(i)>a(j)) then
         t   =a(i)
         a(i)=a(j) 
         a(j)=t
      endif
   enddo
enddo

end subroutine
   
