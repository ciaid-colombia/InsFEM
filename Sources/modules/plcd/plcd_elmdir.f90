module Mod_plcd_elmdir
contains

   subroutine plcd_elmdir(a,e,elmat,elrhs)
      !------------------------------------------------------------------------
      ! This routine modifies the element stiffness to impose the correct 
      ! boundary conditions
      !------------------------------------------------------------------------
      use typre
      use Mod_PLCD
      use Mod_Element
      use Mod_plcd_Stages
      use Mod_plcdExacso
      implicit none
      class(PLCDProblem), target                 :: a
      class(FiniteElement), intent(in)        :: e
      real(rp), intent(inout)                :: elmat(:,:,:,:),elrhs(:,:)
     
      real(rp)                   :: adiag, pressval, pressgrad(3)
      integer(ip)                :: inode,ipoin,idofn,iffix_aux
      type(plcdExacso) :: Exacso
      real(rp), pointer :: coord(:)
     

      !Dirichlet boundary conditions
      do inode=1,e%pnode
         ipoin=e%lnods(inode)
         do idofn=1,a%ndofbc
            iffix_aux = a%cs%kfl_fixno(idofn,ipoin)
            if(iffix_aux==1) then
               adiag=elmat(idofn,inode,idofn,inode)
               
               !Only if the Dirichlet columns are to be deleted
               if (a%kfl_DeleteDirichletColumns) then
                  elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                        - elmat(:,1:e%pnode,idofn,inode)*(a%cs%bvess(idofn,ipoin)*a%css%CurrentLoadFactor - a%Displacement(idofn,ipoin,1))
                  elmat(:,1:e%pnode,idofn,inode)=0.0_rp
               endif

               elrhs(idofn,inode)=adiag*(a%cs%bvess(idofn,ipoin)*a%css%CurrentLoadFactor - a%Displacement(idofn,ipoin,1))
               elmat(idofn,inode,:,1:e%pnode)=0.0_rp
               elmat(idofn,inode,idofn,inode)=adiag
            end if
         end do
      end do
      
      !FixPressure if required
      if (a%UseUPFormulation .and. (a%kfl_confi == 1)) then
         do inode = 1,e%pnode
            ipoin=e%lnods(inode)
            if (ipoin == a%nodpr) then
               adiag=elmat(a%ndofbc+1,inode,a%ndofbc+1,inode)
               if (a%kfl_DeleteDirichletColumns) then
                  elmat(:,1:e%pnode,a%ndofbc+1,inode)=0.0_rp
               endif
               elmat(a%ndofbc+1,inode,:,1:e%pnode)=0.0_rp
               elrhs(a%ndofbc+1,inode) = 0.0_rp
               elmat(a%ndofbc+1,inode,a%ndofbc+1,inode) = adiag
               if (a%kfl_exacs > 0) then
                  call a%Mesh%GetPointCoord(ipoin,coord)
                  call Exacso%ComputeSolution(e%ndime,coord,a)
                  call Exacso%GetPressure(e%ndime,pressval,pressgrad)
                  elrhs(a%ndofbc+1,inode) = adiag*(pressval*a%css%CurrentLoadFactor - a%Pressure(ipoin,1))
               endif
            endif
         enddo
      endif
   end subroutine plcd_elmdir
end module
