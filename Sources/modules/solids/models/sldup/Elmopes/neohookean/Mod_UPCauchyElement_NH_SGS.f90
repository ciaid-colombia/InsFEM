submodule(Mod_CauchyElement) UPCauchyElement_NH_SGS

   implicit none

contains

   module subroutine up_elm_usgs_VS_NH_dynsgs(e,nd,tn,dvol,tautau,divstr,Finv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental matrix 
      !            (v_h,i*(1-tau_k^-1*tau_t), (1/J)F^-1_Kl(du_l/dX_K)div(S_ij),i)
      !            remember: div(S_ij)=(mu/J)*divstr
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip),intent(in)     :: nd,tn
      real(rp),   intent(in)     :: dvol,tautau,Finv(nd,nd),divstr(nd)
      real(rp),   intent(inout)  :: elmat(nd,e%pnode,nd,e%pnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: B(nd,nd)

      do inode=1,e%pnode
          do jnode=1,e%pnode

              call sld_calculateVTauTau_JS(e,nd,tn,inode,jnode,tautau,divstr,Finv,B)
              
              elmat(1:nd,inode,1:nd,jnode) = elmat(1:nd,inode,1:nd,jnode) + B(:,:)*dvol

          end do
      end do

   end subroutine up_elm_usgs_VS_NH_dynsgs

end submodule UPCauchyElement_NH_SGS
