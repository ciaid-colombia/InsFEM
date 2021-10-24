module Mod_LowMachElement_Newton

   use typre
   use Mod_Element
   
contains
   
   subroutine lmn_elmbuv_NR(e,dvolu,acden,gradu,testf_mom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for block U,V 
      !   ((grad U, N)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gradu(e%ndime,e%ndime)
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
      real(rp),    intent(in)    :: dvolu,acden
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)

      integer(ip)                :: jnode,idime,inode,jdime

      do jnode=1,e%pnode
         do idime=1,e%ndime
            do jdime=1,e%ndime
               do inode=1,e%pnode
                  elmat(idime,inode,jdime,jnode) = elmat(idime,inode,idime,jnode) + &
                  acden*(e%shape(inode,e%igaus) + testf_mom(idime,inode)) * &
                  gradu(jdime,idime)*e%shape(jnode,e%igaus)*dvolu
               end do
            end do
         end do
      end do

   end subroutine lmn_elmbuv_NR
   
   subroutine lmn_elmbuq_NR(e,timom,dvolu,acden,gradu,elmat)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gradu(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,timom,acden
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jnode,jdime

      do jnode=1,e%pnode
         do inode=1,e%pnode
            do jdime=1,e%ndime
               do idime=1,e%ndime
                  elmat(1,inode,jdime,jnode) = elmat(1,inode,jdime,jnode) + &
                     e%cartd(idime,inode)*e%shape(jnode,e%igaus)*gradu(jdime,idime)*timom*dvolu*acden*acden
               enddo
            enddo
         end do
      end do

   end subroutine lmn_elmbuq_NR

   subroutine lmn_elmbtv_NR(e,acden,actex,dvol,gradt,ticon,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block T,V 
      !    rho*ticon*alpha*(a.gradT,divv)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)        :: e  
      real(rp),    intent(in)    :: gradt(e%ndime)
      real(rp)                   :: dvol,acden,actex,ticon
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,idime

      do jnode=1,e%pnode
         do inode=1,e%pnode
            do idime=1,e%ndime
               elmat(idime,inode,1,jnode) =  elmat(idime,inode,1,jnode) - &
                  acden*actex*ticon*sum(gradt)*e%shape(jnode,e%igaus)*e%cartd(idime,inode)*dvol
            end do
         end do
      end do
   end subroutine lmn_elmbtv_NR

   subroutine lmn_elmbuw_NR(e,acden,dvol,grat,testf_ene,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for block U,W 
      !    grad(T) *N *N
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)        :: e  
      real(rp),    intent(in)    :: grat(e%ndime)
      real(rp)                   :: dvol, acden, testf_ene(e%pnode)
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode

      do jnode=1,e%pnode
         do inode=1,e%pnode
            elmat(1,inode,1:e%ndime,jnode) = elmat(1,inode,1:e%ndime,jnode) + &
               acden*grat(1:e%ndime)*e%shape(jnode,e%igaus)*(e%shape(inode,e%igaus)+testf_ene(inode))*dvol
         end do
      end do

   end subroutine lmn_elmbuw_NR

   subroutine lmn_elmbuv_NR_rhs(e,dvolu,acden,gradu,gpvel,testf_mom,elrhs)
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvel(e%ndime,*),gradu(e%ndime,e%ndime)
      real(rp),    intent(in)    :: testf_mom(e%ndime,e%pnode)
      real(rp),    intent(in)    :: dvolu,acden
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: idime,inode,jdime

      do idime=1,e%ndime
         do inode=1,e%pnode
            do jdime=1,e%ndime
               elrhs(idime,inode) = elrhs(idime,inode) + &
                  acden*(e%shape(inode,e%igaus) + testf_mom(idime,inode)) * &
                  gradu(jdime,idime)*gpvel(jdime,1)*dvolu
            end do
         end do
      end do

   end subroutine lmn_elmbuv_NR_rhs
   
   subroutine lmn_elmbuq_NR_rhs(e,timom,dvolu,acden,gradu,gpvel,elrhs)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpvel(e%ndime,*),gradu(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,timom,acden
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,idime,jdime

      do inode=1,e%pnode
         do idime=1,e%ndime
            do jdime=1,e%ndime
               elrhs(1,inode) = elrhs(1,inode) + &
                  e%cartd(idime,inode)*gpvel(jdime,1)*gradu(jdime,idime)*timom*dvolu*acden*acden
            enddo
         enddo
      end do

   end subroutine lmn_elmbuq_NR_rhs

   subroutine lmn_elmbtv_NR_rhs(e,acden,actex,dvol,gradt,gpvel,ticon,elrhs)
      implicit none
      class(FiniteElement)        :: e  
      real(rp),    intent(in)    :: gpvel(e%ndime,*),gradt(e%ndime)
      real(rp)                   :: dvol,acden,actex,ticon
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime

      do inode=1,e%pnode
         do idime=1,e%ndime
            do jdime=1,e%ndime
               elrhs(idime,inode) =  elrhs(idime,inode) - &
                  acden*actex*ticon*gradt(jdime)*gpvel(jdime,1)*e%cartd(idime,inode)*dvol
            end do
         end do
      end do
   end subroutine lmn_elmbtv_NR_rhs

   subroutine lmn_elmbuw_NR_rhs(e,acden,dvol,grat,gpvel,testf_ene,elrhs)
      implicit none
      class(FiniteElement)        :: e  
      real(rp),    intent(in)    :: gpvel(e%ndime,*),grat(e%ndime)
      real(rp)                   :: dvol, acden, testf_ene(e%pnode)
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,jdime

      do inode=1,e%pnode
         do jdime=1,e%ndime
            elrhs(1,inode) = elrhs(1,inode) + &
               acden*grat(jdime)*gpvel(jdime,1)*(e%shape(inode,e%igaus)+testf_ene(inode))*dvol
         end do
      end do

   end subroutine lmn_elmbuw_NR_rhs

end module
