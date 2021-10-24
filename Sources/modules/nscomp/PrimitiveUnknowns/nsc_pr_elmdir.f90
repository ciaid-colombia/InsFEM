module Mod_nsc_pr_elmdir

contains

   subroutine nsc_pr_elmdir(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      use Mod_NSCompressiblePrimitive
      implicit none
      class(NSCompressiblePrimitiveProblem)  :: a
      class(FiniteElement), intent(in)        :: e
      real(rp), intent(inout)                :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)
      
      real(rp)                   :: adiag
      integer(ip)                :: inode,ipoin,idofn
      integer(ip)                :: iffix_aux

      !Dirichlet boundary conditions
      do inode=1,e%pnode
         ipoin=e%lnods(inode)
         do idofn=1,a%ndofbc
            iffix_aux = a%kfl_fixno(idofn,ipoin)
            if(iffix_aux==1) then
               adiag=elmat(idofn,inode,idofn,inode)
               
               !Only if the Dirichlet columns are to be deleted
               if (a%kfl_DeleteDirichletColumns) then
                  elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                        - elmat(:,1:e%pnode,idofn,inode)*a%bvess(idofn,ipoin,1)
                  elmat(:,1:e%pnode,idofn,inode)=0.0_rp
               endif
               
               elrhs(idofn,inode)=adiag*a%bvess(idofn,ipoin,1)
               elmat(idofn,inode,:,1:e%pnode)=0.0_rp
               elmat(idofn,inode,idofn,inode)=adiag
               if(adiag<epsilon(0.0_rp)) then 
                  elrhs(idofn,inode)=a%bvess(idofn,ipoin,1)*e%detjm
                  elmat(idofn,inode,idofn,inode)= e%detjm
               end if
            end if
         end do
      end do

   end subroutine      
         
end module


