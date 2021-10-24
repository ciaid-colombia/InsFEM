module Mod_tem_SigmaTerm
contains
subroutine tem_SigmaTerm(e,sig,grvel,dens,sph,elext,term)
   !Computes the source term from the three field problem sigma:grad_sym(u)
   use Mod_Element
   use typre
   implicit none
   class(FiniteElement) :: e
   real(rp) :: sig((e%ndime-1)*(e%ndime-1)+2), grvel(e%ndime,e%ndime), visco         
   real(rp) :: sigcomp(e%ndime,e%ndime)
   real(rp), intent(out) :: elext,term, dens,sph
   integer(ip) :: idime, jdime
  
   term=0.0_rp
   
   call tem_SigmaMatrix(e,sig,sigcomp)
   
   do idime=1,e%ndime
     do jdime=1,e%ndime
         term= term+sigcomp(idime,jdime)*0.5*(grvel(idime,jdime) &
               +grvel(jdime,idime))
     end do
   end do
   
   term=term/(dens*sph)
   
   elext=term+elext
   
end subroutine


subroutine tem_SigmaMatrix(e,sig,sigmatrix)
   !This subroutine convert the vector sigma=(sig_x, sig_y, sig_xy)
   !in a symmetric matrix.
   use Mod_Element
   use typre
   implicit none
   class(FiniteElement) :: e
   real(rp) :: sigmatrix(e%ndime,e%ndime), sig((e%ndime-1)*(e%ndime-1)+2)
   integer(ip) :: idime, jdime, cont
   
   if(e%ndime==2) then 
      cont=0
   else
       cont=e%ndime
   endif
   
   do idime=1,e%ndime
     do jdime=1,e%ndime
        if (idime==jdime) then
            sigmatrix(idime,jdime)=sig(idime)
         elseif(idime<jdime) then
                sigmatrix(idime,jdime)=sig(3+cont)
                sigmatrix(jdime,idime)=sig(3+cont)
             if(e%ndime==3) cont=cont-1
         endif   
     end do
   end do
   
    
end subroutine

end module
