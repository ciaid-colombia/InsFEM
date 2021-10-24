module Mod_NSCompressibleElement
   use typre
   use Mod_Element
   implicit none

   real(rp), parameter :: zeroc = epsilon(0.0_rp)

contains

   subroutine nsc_ComputeElementVelocity(e,elden,elmom,elvel)
      implicit none
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: elden(e%pnode),elmom(e%ndime,e%pnode)
      real(rp),    intent(out)   :: elvel(e%ndime,e%pnode)
      
      integer(ip)                :: idime,inode

      do idime = 1,e%ndime
         do inode = 1,e%pnode
            elvel(idime,inode) = elmom(idime,inode) / elden(inode)
         end do
      end do
      
   end subroutine     

   subroutine nsc_ComputeTau(e,gpspd,kvisc,ktcnd,vnorm,staco,chale,timom)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco(4)
      real(rp),    intent(in)    :: kvisc,ktcnd,gpspd
      real(rp),    intent(in)    :: vnorm,chale(2)
      real(rp),    intent(out)   :: timom(3)  

      real(rp)                   :: fretd,fretm,frete,freqc,freqm,freqe
      integer(ip)                :: idime

      ! Characteristic velocity                                 
      freqc = staco(2)*(vnorm+gpspd)*e%npol/chale(1)      ! Convective term (using h= h_stream)
      freqm = staco(1)*kvisc*e%npol4/(chale(2)*chale(2))  ! Momentum diffusive term (using h=vol^(1/3))
      freqe = staco(1)*ktcnd*e%npol4/(chale(2)*chale(2))  ! Energy diffusive term (using h=vol^(1/3))
      fretd = freqc                                       ! Total frequency for mass
      fretm = freqc + freqm                               ! Total frequency for momentum
      frete = freqc + freqe                               ! Total frequency for energy

      timom = 0.0_rp

      if(fretd>zeroc) then
         timom(1) = 1.0_rp/fretd
      else
         timom(1) = 1.0e12_rp
      end if
      if(fretm>zeroc) then
         timom(2) = 1.0_rp/fretm
      else
         timom(2) = 1.0e12_rp
      end if
      if(frete>zeroc) then
         timom(3) = 1.0_rp/frete
      else
         timom(3) = 1.0e12_rp
      end if

   end subroutine

   subroutine nsc_ComputeSourcesTau(e,gpspd,kvisc,ktcnd,vnorm,noext,noexh,staco,chale,timom)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco(3),chale(2)
      real(rp),    intent(in)    :: kvisc,ktcnd,gpspd
      real(rp),    intent(in)    :: vnorm,noext,noexh
      real(rp),    intent(out)   :: timom(3)  

      real(rp)                   :: fretd,fretm,frete,freqc,freqm,freqe,freqs
      real(rp)                   :: sour1,souroot


      ! Characteristic velocity                                 
      souroot = sqrt(noexh*noexh*noexh*noexh+4.0_rp*gpspd*gpspd*noext*noext*noexh*noexh) 
      sour1 = (souroot + noexh*noexh+2.0_rp*gpspd*gpspd*noext*noext)/(2.0_rp*gpspd*gpspd*gpspd*gpspd)

      freqc = staco(2)*(vnorm+gpspd)*e%npol/chale(1)      ! Convective term (using h= h_stream)
      freqm = staco(1)*kvisc*e%npol4/(chale(2)*chale(2))  ! Momentum diffusive term (using h=vol^(1/3))
      freqe = staco(1)*ktcnd*e%npol4/(chale(2)*chale(2))  ! Energy diffusive term (using h=vol^(1/3))
      freqs = staco(3)*sqrt(sour1)                        ! Sources term

      fretd = freqc + freqs                                     ! Total frequency for mass
      fretm = freqc + freqm + freqs                             ! Total frequency for momentum
      frete = freqc + freqe + freqs                             ! Total frequency for energy

      timom = 0.0_rp

      if(fretd>zeroc) then
         timom(1) = 1.0_rp/fretd
      else
         timom(2) = 1.0e12_rp
      end if
      if(fretm>zeroc) then
         timom(2) = 1.0_rp/fretm
      else
         timom(2) = 1.0e12_rp
      end if
      if(frete>zeroc) then
         timom(3) = 1.0_rp/frete
      else
         timom(3) = 1.0e12_rp
      end if

   end subroutine

   subroutine nsc_vgvar(e,elvar,vgrad,vgvar)
      !-----------------------------------------------------------------------
      !
      ! This routine computes vel·grad var
      !
      !-----------------------------------------------------------------------
     implicit none
     class(FiniteElement) :: e
     
     real(rp),    intent(in)    :: elvar(e%mnode)
     real(rp),    intent(in)    :: vgrad(e%mnode)
     real(rp),    intent(inout) :: vgvar
 
     vgvar = dot_product(elvar,vgrad)

   end subroutine nsc_vgvar

   subroutine nsc_vgvec(e,elvec,vgrad,vgvec)
      !-----------------------------------------------------------------------
      !
      ! This routine computes vel·grad vec
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: elvec(e%ndime,e%mnode)
      real(rp),    intent(in)    :: vgrad(e%mnode)
      real(rp),    intent(inout) :: vgvec(e%ndime)
   
      integer(ip)                :: idime
   
   
      do idime=1,e%ndime
   
         vgvec(idime) = dot_product(elvec(idime,:),vgrad)
     
      end do

   end subroutine nsc_vgvec

   subroutine nsc_TimeIntegrationToElTemp(e,ndime,Integrator,dtinv,gpvar,eltemp)
      use typre
      use Mod_Element
      use Mod_TimeIntegrator
      implicit none
      class(FiniteElement)        :: e
      integer(ip)                :: ndime
      type(TimeIntegratorDt1)    :: Integrator
      real(rp)                   :: gpvar(ndime,*),eltemp(ndime)
      real(rp)                   :: dtinv
      
      real(rp)                   :: gprhs(ndime)
      
      !Time integration
      call Integrator%GetRHS(ndime,gpvar(:,2),gprhs)
      eltemp = eltemp + dtinv*gprhs
      
      
   end subroutine

   subroutine nsc_ArtificialVis(e,shock,chale,elrem,grmom,arvis)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: shock,chale(2)
      real(rp),    intent(in)    :: elrem(e%ndime),grmom(e%ndime,e%ndime)
      real(rp),    intent(out)   :: arvis

      real(rp)                   :: aux,rnor,grnrm


      aux = 0.5_rp * shock * chale(2)

      call vecnor(elrem,e%ndime,rnor,2_ip)

      call matNormF(e%ndime,grmom,grnrm)

      if(grnrm>zeroc) then
         arvis = aux * rnor / grnrm
      else
         arvis = 0.0_rp
      end if
      

   end subroutine nsc_ArtificialVis

   subroutine nsc_ArtificialCnd(e,shock,chale,elree,grene,artco)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: shock,chale(2)
      real(rp),    intent(in)    :: elree,grene(e%ndime)
      real(rp),    intent(out)   :: artco

      real(rp)                   :: aux,gnor


      aux = 0.5_rp * shock * chale(2)

      call vecnor(grene,e%ndime,gnor,2_ip)

      if(gnor>zeroc) then
         artco = aux * abs(elree) / gnor
      else
         artco = 0.0_rp
      end if
      
   end subroutine nsc_ArtificialCnd

   subroutine nsc_GradOrthVis(e,shock,chale,gpvno,orthm,grmom,arvis)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: shock,chale(2),gpvno
      real(rp),    intent(in)    :: orthm(e%ndime,e%ndime),grmom(e%ndime,e%ndime)
      real(rp),    intent(out)   :: arvis

      real(rp)                   :: aux,ornrm,grnrm

      aux = 0.5_rp * shock * chale(2) * gpvno

      call matNormF(e%ndime,orthm,ornrm)

      call matNormF(e%ndime,grmom,grnrm)

      if(grnrm>zeroc) then
         arvis = aux * ornrm / grnrm
      else
         arvis = 0.0_rp
      end if
      

   end subroutine nsc_GradOrthVis

   subroutine nsc_GradOrthCnd(e,shock,chale,gpvno,orthe,grene,artco)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: shock,chale(2),gpvno
      real(rp),    intent(in)    :: orthe(e%ndime),grene(e%ndime)
      real(rp),    intent(out)   :: artco

      real(rp)                   :: aux,ornrm,gnor


      aux = 0.5_rp * shock * chale(2) * gpvno 

      call vecnor(orthe,e%ndime,ornrm,2_ip)

      call vecnor(grene,e%ndime,gnor,2_ip)

      if(gnor>zeroc) then
         artco = aux * ornrm / gnor
      else
         artco = 0.0_rp
      end if
      
   end subroutine nsc_GradOrthCnd

   subroutine nsc_StreamlineTensor(ndime,unitv,elstr)
      !-----------------------------------------------------------------------
      !
      ! This routine computes a streamline tensor with respect to unity vector U
      !    StreamlineTensor = UU
      !
      !-----------------------------------------------------------------------
      implicit none
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: unitv(ndime)
      real(rp),    intent(out)   :: elstr(ndime,ndime)

      integer(ip)                :: idime,jdime

      do idime = 1,ndime
         do jdime = 1,ndime
            elstr(idime,jdime) = unitv(idime)*unitv(jdime)
         enddo
      enddo

   end subroutine nsc_StreamlineTensor

   subroutine nsc_OrthogonalTensor(ndime,elstr,elort)
      !-----------------------------------------------------------------------
      !
      ! This routine computes an orthogonal tensor using the streamline tensor
      !    OrthogonalTensor = I - UU
      !
      !-----------------------------------------------------------------------
      implicit none
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: elstr(ndime,ndime)
      real(rp),    intent(out)   :: elort(ndime,ndime)

      integer(ip)                :: idime,jdime

      do idime = 1,ndime
         elort(idime,idime) = 1.0_rp
         do jdime = 1,ndime
            elort(idime,jdime) = elort(idime,jdime) - elstr(idime,jdime)
         enddo
      enddo

   end subroutine nsc_OrthogonalTensor

   subroutine nsc_FourthStreamOrthogonalTensor(ndime,untmv,forstr,forort)
      !-----------------------------------------------------------------------
      !
      ! This routine computes 4th order tensors
      !
      !-----------------------------------------------------------------------
      implicit none
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: untmv(ndime)
      real(rp),    intent(inout) :: forstr((ndime-1)*(ndime-1)+2,(ndime-1)*(ndime-1)+2)
      real(rp),    intent(inout) :: forort((ndime-1)*(ndime-1)+2,(ndime-1)*(ndime-1)+2)

      real(rp)                :: S_11, S_12, S_13
      real(rp)                :: S_22, S_23, S_33

      S_11 = untmv(1)*untmv(1)
      S_22 = untmv(2)*untmv(2)
      S_12 = untmv(1)*untmv(2)

      forstr(1,1) = S_11
      forstr(2,2) = S_22
      forstr(1,2) = S_12
      forstr(2,1) = forstr(1,2)

      forort(1,1) = 1.0_rp-S_11
      forort(2,2) = 1.0_rp-S_22
      forort(1,2) = -S_12
      forort(2,1) = forort(1,2)

      if(ndime.eq.2) then
         forstr(3,3) = S_12
         forort(3,3) = 1.0_rp-S_12
      elseif(ndime.eq.3) then
         S_33 = untmv(3)*untmv(3)
         S_13 = untmv(1)*untmv(3)
         S_23 = untmv(2)*untmv(3)
         forstr(3,3) = S_33
         forstr(1,3) = S_13
         forstr(2,3) = S_23
         forstr(3,1) = forstr(1,3)
         forstr(3,2) = forstr(2,3)
         forstr(4,4) = S_23
         forstr(5,5) = S_13
         forstr(6,6) = S_12

         forort(3,3) = 1.0_rp-S_33
         forort(1,3) = 1.0_rp-S_13
         forort(2,3) = 1.0_rp-S_23
         forort(3,1) = forort(1,3)
         forort(3,2) = forort(2,3)
         forstr(4,4) = 1.0_rp-S_23
         forstr(5,5) = 1.0_rp-S_13
         forstr(6,6) = 1.0_rp-S_12
      end if

   end subroutine nsc_FourthStreamOrthogonalTensor

   subroutine nsc_AnisotropicVisTen(ndime,acvis,arvis,arsvs,gpode,untmv,novst,elvst)
      !-----------------------------------------------------------------------
      !
      ! This routine computes anisotropic viscosity tensor using 4th order tensors
      !    visten = gpode*(acvis*I+arsvs*strea+arvis*orto)*visten
      !
      !-----------------------------------------------------------------------
      implicit none
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: acvis,arvis,arsvs
      real(rp),    intent(in)    :: gpode
      real(rp),    intent(in)    :: untmv(ndime)
      real(rp),    intent(in)    :: novst(ndime,ndime)
      real(rp),    intent(inout) :: elvst(ndime,ndime)

      real(rp)                :: strea(ndime,ndime)
      real(rp)                :: ortho(ndime,ndime)

      real(rp)                :: S_11, S_12, S_13
      real(rp)                :: S_22, S_23, S_33

      strea = 0.0_rp
      ortho = 0.0_rp

      S_11 = untmv(1)*untmv(1)
      S_22 = untmv(2)*untmv(2)
      S_12 = untmv(1)*untmv(2)

      strea(1,1) = S_11*novst(1,1)+S_12*novst(2,2)
      strea(2,2) = S_12*novst(1,1)+S_22*novst(2,2)
      strea(1,2) = S_12*novst(1,2)
      strea(2,1) = strea(1,2)

      ortho(1,1) = (1.0_rp-S_11)*novst(1,1)-S_12*novst(2,2)
      ortho(2,2) = -S_12*novst(1,1)+(1.0_rp-S_22)*novst(2,2)
      ortho(1,2) = (1.0_rp-S_12)*novst(1,2)
      ortho(2,1) = ortho(1,2)

      if(ndime.eq.3) then
         S_33 = untmv(3)*untmv(3)
         S_13 = untmv(1)*untmv(3)
         S_23 = untmv(2)*untmv(3)
         strea(1,1) = strea(1,1) + S_13*novst(3,3)
         strea(2,2) = strea(2,2) + S_23*novst(3,3)
         strea(3,3) = S_13*novst(1,1)+S_23*novst(2,2)+S_33*novst(3,3)
         strea(1,3) = S_13*novst(1,3)
         strea(2,3) = S_23*novst(2,3)
         strea(3,1) = strea(1,3)
         strea(3,2) = strea(2,3)
   
         ortho(1,1) = ortho(1,1) - S_13*novst(3,3)
         ortho(2,2) = ortho(2,2) - S_23*novst(3,3)
         ortho(3,3) = -S_13*novst(1,1)-S_23*novst(2,2)+(1.0_rp-S_33)*novst(3,3)
         ortho(1,3) = (1.0_rp-S_13)*novst(1,3)
         ortho(2,3) = (1.0_rp-S_23)*novst(2,3)
         ortho(3,1) = ortho(1,3)
         ortho(3,2) = ortho(2,3)
      end if

      elvst = gpode*(acvis*novst+arsvs*strea+arvis*ortho)

   end subroutine nsc_AnisotropicVisTen

   subroutine nsc_HessLap(e,hess,lapl)
    !-----------------------------------------------------------------------
    !
    ! This routine writes the Hessian into a expanded matrix
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(out) :: hess(e%ndime,e%ndime,e%pnode)
    real(rp),    intent(out) :: lapl(e%pnode)

    integer(ip) :: idime,itens

    if(e%ndime.eq.2) then
      hess(1,1,1:e%pnode)=e%hessi(1,1:e%pnode)
      hess(1,2,1:e%pnode)=e%hessi(3,1:e%pnode)
      hess(2,1,1:e%pnode)=e%hessi(3,1:e%pnode)
      hess(2,2,1:e%pnode)=e%hessi(2,1:e%pnode)
      lapl(1:e%pnode)=e%hessi(1,1:e%pnode)+e%hessi(2,1:e%pnode)
    else
      hess(1,1,1:e%pnode)=e%hessi(1,1:e%pnode)
      hess(1,2,1:e%pnode)=e%hessi(4,1:e%pnode)
      hess(2,1,1:e%pnode)=e%hessi(4,1:e%pnode)
      hess(2,2,1:e%pnode)=e%hessi(2,1:e%pnode)
      hess(1,3,1:e%pnode)=e%hessi(5,1:e%pnode)
      hess(3,1,1:e%pnode)=e%hessi(5,1:e%pnode)
      hess(2,3,1:e%pnode)=e%hessi(6,1:e%pnode)
      hess(3,2,1:e%pnode)=e%hessi(6,1:e%pnode)
      hess(3,3,1:e%pnode)=e%hessi(3,1:e%pnode)
      lapl(1:e%pnode)=e%hessi(1,1:e%pnode)+e%hessi(2,1:e%pnode)+e%hessi(3,1:e%pnode)
    end if

   end subroutine nsc_HessLap

end module
