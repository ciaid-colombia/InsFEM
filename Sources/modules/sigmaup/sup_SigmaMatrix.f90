subroutine sup_SigmaMatrix(ndime,tens,sig,sigmatrix)
   !This subroutine convert the vector sigma=(sig_x, sig_y, sig_xy)
   !in a symmetric matrix.
   use Mod_Element
   use typre
   implicit none
   integer(ip) :: idime, jdime, cont,ndime,tens
   real(rp) :: sigmatrix(ndime,ndime), sig(tens)
   
   
   if (ndime==2) then
      sigmatrix(1,1)=sig(1)
      sigmatrix(2,2)=sig(2)
      sigmatrix(1,2)=sig(3)
      sigmatrix(2,1)=sig(3)
   else if (ndime==3) then
      sigmatrix(1,1)=sig(1)
      sigmatrix(2,2)=sig(2)
      sigmatrix(3,3)=sig(3)
      sigmatrix(2,3)=sig(4)
      sigmatrix(3,2)=sig(4)
      sigmatrix(1,3)=sig(5)
      sigmatrix(3,1)=sig(5)
      sigmatrix(1,2)=sig(6)
      sigmatrix(2,1)=sig(6)
   end if   
      
    
end subroutine
