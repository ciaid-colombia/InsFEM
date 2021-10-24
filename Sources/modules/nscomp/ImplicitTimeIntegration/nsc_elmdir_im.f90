module Mod_nsc_elmdir_im

contains

   subroutine nsc_elmdir(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      use Mod_Mesh
      use Mod_NSCompressibleImplicit
      use Mod_NSCompressibleSubroutines
      implicit none
      class(NSCompressibleImplicitProblem)             :: a
      class(FiniteElement), intent(in)        :: e
      real(rp), intent(inout)                :: elmat(:,:,:,:), elrhs(:,:)
      
      integer(ip)                :: inode,ipoin,idofn
      integer(ip)                :: iffix_den,iffix_tem,iffix_vel
      real(rp)                   :: adiag
      real(rp)                   :: energ,acvis,actco,accph,accvh


      call a%GetPhysicalParameters(acvis,actco,accph,accvh)

      !Dirichlet boundary conditions
      do inode=1,e%pnode
         ipoin=e%lnods(inode)
         iffix_den = a%kfl_fixno(1,ipoin)
         if(iffix_den==1) then
            !density value is multiplied in the dirichlet values for momentum and energy
            a%densf(ipoin,1) = a%bvess(1,ipoin,1)
            adiag=elmat(1,inode,1,inode)
            !Only if the Dirichlet columns are to be deleted
            if (a%kfl_DeleteDirichletColumns) then
                elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                      - elmat(:,1:e%pnode,1,inode)*a%bvess(1,ipoin,1)
                elmat(:,1:e%pnode,1,inode)=0.0_rp
            endif
            elrhs(1,inode)=adiag*a%bvess(1,ipoin,1)
            elmat(1,inode,:,1:e%pnode)=0.0_rp
            elmat(1,inode,1,inode)=adiag
            if(adiag<epsilon(0.0_rp)) then 
               elrhs(1,inode)=a%bvess(1,ipoin,1)*e%detjm
               elmat(1,inode,1,inode)= e%detjm
            end if
         end if
         do idofn=1,a%ndofbc-2
            iffix_vel = a%kfl_fixno(idofn+1,ipoin)
            if(iffix_vel==1) then
               !velocity value is multiplied in the dirichlet value for energy
               a%veloc(idofn,ipoin,1) = a%bvess(idofn+1,ipoin,1)
               adiag=elmat(idofn+1,inode,idofn+1,inode)
               !Only if the Dirichlet columns are to be deleted
               if (a%kfl_DeleteDirichletColumns) then
                  elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                        - elmat(:,1:e%pnode,idofn+1,inode)*a%densf(ipoin,1)*a%bvess(idofn+1,ipoin,1)
                  elmat(:,1:e%pnode,idofn+1,inode)=0.0_rp
               endif
               
               elrhs(idofn+1,inode)=adiag*a%densf(ipoin,1)*a%bvess(idofn+1,ipoin,1)
               elmat(idofn+1,inode,:,1:e%pnode)=0.0_rp
               elmat(idofn+1,inode,idofn+1,inode)=adiag
            end if
         end do
         iffix_tem = a%kfl_fixno(a%ndofbc,ipoin)
         if(iffix_tem==1) then
            a%tempe(ipoin,1)= a%bvess(a%ndofbc,ipoin,1)
            call nsc_ComputeEnergy(e%ndime,accvh,a%tempe(ipoin,1),a%veloc(:,ipoin,1),energ)

            adiag=elmat(a%ndofbc,inode,a%ndofbc,inode)
               
            !Only if the Dirichlet columns are to be deleted
            if (a%kfl_DeleteDirichletColumns) then
               elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                     - elmat(:,1:e%pnode,a%ndofbc,inode)*a%densf(ipoin,1)*energ
               elmat(:,1:e%pnode,a%ndofbc,inode)=0.0_rp
            endif
            
            elrhs(a%ndofbc,inode)=adiag*a%densf(ipoin,1)*energ
            elmat(a%ndofbc,inode,:,1:e%pnode)=0.0_rp
            elmat(a%ndofbc,inode,a%ndofbc,inode)=adiag
         end if
      end do

   end subroutine      
         
end module


