module Mod_php_elmdir
   use typre
   use Mod_Element
   use Mod_PhysicalProblem
   implicit none
contains

   subroutine php_elmdir(a,e,ndofn,ndofbc,ndofbcstart,currentbvess,elmat,elrhs)
      !------------------------------------------------------------------------
      ! This routine modifies the element stiffness to impose the correct 
      ! boundary conditions
      !------------------------------------------------------------------------
      implicit none
      class(PhysicalProblem)     :: a
      class(FiniteElement)       :: e
      integer(ip)                :: ndofn,ndofbcstart,ndofbc
      real(rp)                   :: elmat(:,:,:,:),elrhs(:,:)
      
      integer(ip), intent(in)    :: currentbvess
      real(rp)                   :: adiag
      integer(ip)                :: inode,ipoin,idofn
      integer(ip)                :: iffix_aux

      !Dirichlet boundary conditions
      do inode=1,e%pnode
         ipoin=e%lnods(inode)
         do idofn=1,ndofbc
            iffix_aux = a%kfl_fixno(idofn+currentbvess-1,ipoin)
            if(iffix_aux==1) then
               adiag=elmat(idofn+ndofbcstart,inode,idofn+ndofbcstart,inode)
               
               !Only if the Dirichlet columns are to be deleted
               if (a%kfl_DeleteDirichletColumns) then
                  elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                        - elmat(:,1:e%pnode,idofn+ndofbcstart,inode)*a%bvess(currentbvess+idofn-1,ipoin,1)
                  elmat(:,1:e%pnode,ndofbcstart+idofn,inode)=0.0_rp
               endif
               
               elrhs(idofn+ndofbcstart,inode)=adiag*a%bvess(currentbvess+idofn-1,ipoin,1)
               elmat(ndofbcstart+idofn,inode,:,1:e%pnode)=0.0_rp
               elmat(ndofbcstart+idofn,inode,ndofbcstart+idofn,inode)=adiag
            end if
         end do
      end do
   end subroutine php_elmdir
end module
