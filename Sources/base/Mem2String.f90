 subroutine Mem2String(memory,outmemo,rbyte,outstr)
      use typre
      implicit none
      integer(8) :: memory
      real(rp) :: outmemo
      character(6) :: outstr
      real(rp) :: rbyte
      
       if(abs(memory)>=1024*1024*1024) then
         rbyte=1024.0_rp*1024.0_rp*1024.0_rp
         outstr='Gbytes'
      else if(abs(memory)>=1024*1024) then 
         rbyte=1024.0_rp*1024.0_rp
         outstr='Mbytes'     
      else if(abs(memory)>=1024) then 
         rbyte=1024.0_rp
         outstr='kbytes'          
      else  
         rbyte=1.0_rp
         outstr=' bytes'     
      end if
      outmemo = real(memory)/rbyte
         
   
end subroutine
  