subroutine sldsup_GetPhysicalParameters_lin(a,ndime,sz,densi,K,G,C_dev,D,D_dev,Dv_scalar,id,ielem)
   use typre
   use Mod_SUPSolids_lin
   implicit none
   class(SUPSolidsProblem_lin), target :: a
   integer(ip),intent(in) :: ndime,sz
   real(rp),   intent(out):: densi,K,G,Dv_scalar
   real(rp),   intent(out):: C_dev(sz,sz),D(sz,sz),D_dev(sz,sz),id(sz,sz)
   integer(ip),optional   :: ielem
   real(rp)               :: D_vol(sz,sz)
   real(rp)               :: D_aux(sz,1)
   real(rp)               :: ii_tt(sz,sz),iv(sz,sz)
   
   densi = a%densi
   K     = a%bulkMod
   G     = a%mu

   call get4thIITensor(ndime,sz,ii_tt)
   call get4thIIVolumetricTensor(ndime,sz,iv)

   !Also known as P
   id = ii_tt - iv

   C_dev = (2.0_rp*G)*id

   D_vol = (ndime/K)*iv
   D_dev = (1.0_rp/(2.0_rp*G))*id

   D     = D_vol + D_dev 

   if(ndime==2) then
       Dv_scalar= 16.0_rp/(27.0_rp*K)
   else 
       Dv_scalar= 1.0_rp/K
   endif

end subroutine
