module Mod_NavierStokesElement_Residuals
   use typre
   use Mod_Element
   
contains

   !This subroutine assemblies the contribution of the residual to the RHS for projecting it
   subroutine nsm_elmrep(e,ndime,dvol,gprep,elrep)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dvol,gprep(ndime)
      real(rp),    intent(inout) :: elrep(ndime,e%mnode)
      integer(ip)                :: ndime,inode
      
      do inode = 1,e%pnode
         elrep(:,inode) = elrep(:,inode) + e%shape(inode,e%igaus)*gprep(:)*dvol
      enddo
   end subroutine  

   !This suborutine computes the residual at a Gauss point
   subroutine nsm_elmrfe_oto(e,dtnsi,acden,veadv,vegau,grapr,grave,elext,eltemp,resim)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dtnsi,acden,veadv(e%ndime),vegau(e%ndime)
      real(rp),    intent(in)    :: grapr(e%ndime),grave(e%ndime,e%ndime)
      real(rp),    intent(in)    :: elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(inout) :: resim(e%ndime+1)
      integer(ip)                :: idime,jdime
      
      resim = 0.0_rp
      resim(1:e%ndime) = resim(1:e%ndime) + acden*dtnsi*vegau !Temporal derivative LHS
      resim(1:e%ndime) = resim(1:e%ndime) - elext - eltemp !Temporal derivative + External force RHS
      resim(1:e%ndime) = resim(1:e%ndime) + grapr !Pressure gradient
      do jdime = 1,e%ndime
         resim(1:e%ndime) = resim(1:e%ndime) + acden*veadv(jdime)*grave(:,jdime) !Convective term
      end do
      do idime = 1,e%ndime
         resim(e%ndime+1) = resim(e%ndime+1) + grave(idime,idime) !Divergence term
      end do
   end subroutine nsm_elmrfe_oto
   
   subroutine nsm_elmrfe_oto_nonlinear(e,acvis,elvel,gprep)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: acvis,elvel(e%ndime,e%pnode)
      real(rp),    intent(inout) :: gprep(e%ndime)
      integer(ip)                :: inode,idime
      real(rp)                   :: aux1(e%pnode)
      
      do inode = 1,e%pnode
         aux1(inode) = sum(e%hessi(1:e%ndime,inode))
      enddo   
      do idime = 1,e%ndime
         gprep(idime) = gprep(idime) - acvis*dot_product(aux1,elvel(idime,1:e%pnode))
      end do
   end subroutine
   
   !This subroutine computes the residual at a Gauss Point, but it does not assembly the
   !terms which are on the finite element space (elext, temporal derivatives)
   subroutine nsm_elmrfe_trm(e,acden,veadv,grapr,grave,resim)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: acden,veadv(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime)
      real(rp),    intent(inout) :: resim(e%ndime+1)
      integer(ip)                :: idime,jdime
      
      resim = 0.0_rp
      resim(1:e%ndime) = resim(1:e%ndime) + grapr  !Pressure gradient
      do jdime = 1,e%ndime
         resim(1:e%ndime) = resim(1:e%ndime) + acden*veadv(jdime)*grave(:,jdime) !Convective term 
      end do
      do idime = 1,e%ndime
         resim(e%ndime+1) = resim(e%ndime+1) + grave(idime,idime) !Divergence term
      end do
   end subroutine nsm_elmrfe_trm
   
   subroutine nsm_elmrfe_split(e,acden,veadv,grapr,grave,resim)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: acden,veadv(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime)
      real(rp),    intent(inout) :: resim(2*e%ndime+1)
      integer(ip)                :: idime,jdime
      
      resim = 0.0_rp
      do jdime = 1,e%ndime
         resim(1:e%ndime) = resim(1:e%ndime) + acden*veadv(jdime)*grave(:,jdime) 
      end do
      do idime = 1,e%ndime
         resim(e%ndime+1) = resim(e%ndime+1) + grave(idime,idime)
      end do
      resim(e%ndime+2:2*e%ndime+1) = resim(e%ndime+2:2*e%ndime+1) + grapr
   end subroutine nsm_elmrfe_split
      
   subroutine nsm_elmrfe_oto_coriolis(e,acden,AxesAngularVeloc,vegau,gprep)
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: acden,AxesAngularVeloc(e%ndime),vegau(e%ndime)
      real(rp),    intent(inout) :: gprep(e%ndime)
      
      gprep(1) = gprep(1) +  acden*(AxesAngularVeloc(2)*vegau(3) - AxesAngularVeloc(3)*vegau(2))
      gprep(2) = gprep(2) +  acden*(AxesAngularVeloc(3)*vegau(1) - AxesAngularVeloc(1)*vegau(3))
      gprep(3) = gprep(3) +  acden*(AxesAngularVeloc(1)*vegau(2) - AxesAngularVeloc(2)*vegau(1))
   end subroutine

end module
