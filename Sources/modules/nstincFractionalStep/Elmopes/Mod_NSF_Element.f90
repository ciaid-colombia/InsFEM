module Mod_NSF_Element
   use typre
   use Mod_Element
   use Mod_NSFractionalStep
   implicit none

contains

   subroutine nsm_elmrhuf(e,dvolu,testf,grpre,gppre,elext,eltemp,elrhs)
		!-----------------------------------------------------------------------------------
      !rhs in nsf1_elmope, additional terms involving the pressure for OSS
      !                  (div v_h, p^n) - tau_m (a * grad v_h, grad p^n)
      !----------------------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: testf(e%pnode),elext(e%ndime),eltemp(e%ndime)
      real(rp),    intent(in)    :: dvolu,gppre,grpre(e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      
      real(rp) :: aux2
      integer(ip) :: inode
 
      do inode=1,e%pnode
         aux2 = testf(inode)*dvolu
         elrhs(:,inode) = -grpre(:)*aux2 + gppre*e%cartd(:,inode)*dvolu & 
               + elrhs(:,inode)
      end do
   end subroutine
   
   subroutine nsm_elmrhuf_split(e,dvolu,gppre,elrhs)
		!-------------------------------------------------------------------------------------
		! rhs in nsf1_elmope, additional terms involving the pressure for split OSS
		!                             (div v_h, p^n)
		!------------------------------------------------------------------------------------
		
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu,gppre
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      
      integer(ip) :: inode
 
      do inode=1,e%pnode
         elrhs(:,inode) = gppre*e%cartd(:,inode)*dvolu + elrhs(:,inode)
      end do

   end subroutine
   
   subroutine nsf_elmdir_press(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      class(NSFractionalStepProblem) :: a
      real(rp) :: elmat(1,e%mnode,1,e%mnode), elrhs(1,e%mnode)
      
      real(rp) :: adiag
      integer(ip) :: inode,ipoin,jnode
      real(rp)    :: prepr
      
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         
         if (ipoin == a%nodpr .or. a%kfl_fixpr(ipoin) /= 0) then
            prepr = a%prepr
            adiag=elmat(1,inode,1,inode)
            elmat(1,inode,1,1:e%pnode)=0.0_rp
            do jnode = 1,e%pnode
               elrhs(1,jnode) = elrhs(1,jnode)- elmat(1,jnode,1,inode)*prepr
            enddo
            elmat(1,1:e%pnode,1,inode) = 0.0_rp
            elmat(1,inode,1,inode) = adiag
            elrhs(1,inode) = adiag*prepr
         endif
      enddo
   end subroutine 
    
   subroutine nsm_elmrhf_pre(e,dvolu,acden,timom,dtinv,divvel,gprep,grpre,&
       grvel,gpadv,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for p_n+1 for Fractional step
      ! -(div u, q)-tau*(u · grad u, grad q)+tau*(proj, grad q)+(1/rho)*dt*(grad p_n, grad q)
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: acden,gprep(e%ndime+1),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: gpadv(e%ndime),grpre(e%ndime)
      real(rp),    intent(in)    :: dvolu,timom,dtinv,divvel
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: tmp,tmp2
      real(rp)                   :: aux1

      tmp = timom*dvolu
      
      aux1 = dtinv*acden
      tmp2= dvolu*(1.0_rp/aux1)

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = (tmp*gprep(idime)+tmp2*grpre(idime))*e%cartd(idime,inode) &
                  + elrhs(1,inode)
            do jdime=1,e%ndime
               elrhs(1,inode) = -tmp*acden*gpadv(jdime)*grvel(idime,jdime)*e%cartd(idime,inode) &
                     + elrhs(1,inode)
            end do
         end do
         elrhs(1,inode) = - dvolu*divvel*e%shape(inode,e%igaus) + elrhs(1,inode)
      end do

   end subroutine nsm_elmrhf_pre
   
	subroutine nsm_elmrhf_pre_split(e,dvolu,acden,timom,dtinv,divvel,gprep,grpre,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for p_n+1 for Fractional step in split OSS
      ! -(div u, q)+tau*(proj, grad q)+(1/rho)*dt*(grad p_n, grad q)
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: acden,gprep(e%ndime)
      real(rp),    intent(in)    :: grpre(e%ndime)
      real(rp),    intent(in)    :: dvolu,timom,dtinv,divvel
      real(rp),    intent(inout) :: elrhs(1,e%pnode)

      integer(ip)                :: inode,idime
      real(rp)                   :: tmp,tmp2
      real(rp)                   :: aux1

      tmp = timom*dvolu
      
      aux1 = dtinv*acden
      tmp2= dvolu*(1.0_rp/aux1)

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = (tmp*gprep(idime)+tmp2*grpre(idime))*e%cartd(idime,inode) &
                  + elrhs(1,inode)
         end do
         elrhs(1,inode) = - dvolu*divvel*e%shape(inode,e%igaus) + elrhs(1,inode)
      end do

   end subroutine nsm_elmrhf_pre_split
   

   subroutine nsf_hydro_rhs(e,dvolu,elext,elrhs)
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: elext(e%ndime)
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      
      integer(ip) :: inode
      
      do inode = 1,e%pnode
         elrhs(1,inode) = elrhs(1,inode) + dot_product(elext(1:e%ndime),e%cartd(1:e%ndime,inode))*dvolu
      enddo
   end subroutine
   
   subroutine nsf_elm2nd_reslap(e,dvolu,timom,acvis,gplapu,elrhs)
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: dvolu
      real(rp)                   :: timom
      real(rp)                   :: acvis
      real(rp)                   :: gplapu(e%ndime)
      real(rp), intent(inout)    :: elrhs(1,e%mnode)
      
      real(rp) :: aux1
      integer(ip) :: inode
      
      aux1 = timom*acvis*dvolu
      do inode = 1,e%pnode
         elrhs(1,inode) = elrhs(1,inode) + aux1*dot_product(e%cartd(1:e%ndime,inode),gplapu(1:e%ndime))
      enddo   
   end subroutine  

   subroutine nsf_elmcomp(e,dvolt0,visac,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for the volumetric viscosity 
    !     (1/3)*mu*grad(div(u))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement)   :: e
    real(rp)                   :: dvolt0
    real(rp),    intent(in)    :: visac
    real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)

    integer(ip)                :: inode,jnode,idime,jdime
    real(rp)                   :: tmp

    tmp = dvolt0*visac
    do jnode=1,e%pnode
      do jdime=1,e%ndime
         do inode=1,e%pnode
            do idime=1,e%ndime
               elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) + (1.0/3.0)*tmp*e%cartd(jdime,inode)*e%cartd(idime,jnode)
            enddo
         enddo
      enddo
    enddo
   end subroutine 

end module
