module Mod_ConvectiveElement
   use typre
   use Mod_Element
   implicit none
   
   real(rp), parameter :: zeroc = epsilon(0.0_rp)
   
contains

   subroutine elmchl(e,kfl_advec,elvel,chale,kfl_hdifumin)
      !-----------------------------------------------------------------------
      !    This routine computes the characteristic element lengths. The first one
      !    is the length in the flow direction if there are convective terms and
      !    the minimum length otherwise. The second one is the square root in 2D 
      !    or cubic root in 3D of the volume.
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)   :: e
      integer(ip), intent(in)    :: kfl_advec
      real(rp),    intent(in)    :: elvel(e%ndime,e%mnode)
      real(rp),    intent(out)   :: chale(2)
      
      integer(ip), optional :: kfl_hdifumin
      
      integer(ip) :: hdifumin
      real(rp) :: chave(e%ndime,2)
      integer(ip)                :: idime,inode
      real(rp)                   :: elno1,elno2

      if (.not. present(kfl_hdifumin)) then
         hdifumin = 0
      else
         hdifumin = kfl_hdifumin
      endif
      
      !Advective element lenght
      !Initialization
      chale(1)=e%hleng(e%ndime)  
      ! Length in the flow direction
      if(kfl_advec>=1) then
         chave=0.0_rp                                ! Characteristic element velocity
         do idime=1,e%ndime
            do inode=1,e%pnode
               chave(idime,1)=chave(idime,1)+elvel(idime,inode)
            end do
            chave(idime,1)=chave(idime,1)/real(e%pnode)
         end do
         call mbvab1(chave(1,2),e%tragl,chave(1,1),e%ndime,e%ndime,elno2,elno1)
         if(elno2>zeroc.and.elno1>zeroc) &
            chale(1)=e%hnatu*elno1/elno2               ! Characteristic element length

      end if
      
      if (hdifumin == 0) then
         !The second one is square or cubic root of the volume
         chale(2)=e%hleng(1)
         do idime=2,e%ndime
            chale(2)=chale(2)*e%hleng(idime)
         end do
         chale(2)=chale(2)**(1.0_rp/real(e%ndime))
      else
         chale(2) = e%hleng(e%ndime)
      endif

   end subroutine elmchl
   
   subroutine ComputeAGradV(e,gpvel,vtemp)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: gpvel(e%ndime)
      real(rp) :: vtemp(e%pnode)
      integer(ip) :: inode
      
      do inode=1,e%pnode
         vtemp(inode) = dot_product(gpvel,e%cartd(:,inode))
      end do
   end subroutine 
   
   subroutine ComputeTau(e,denac,visac,vnorm,staco,chale,timom)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco(2)
      real(rp),    intent(in)    :: denac,visac
      real(rp),    intent(in)    :: vnorm,chale(2)
      real(rp),    intent(out)   :: timom
      real(rp)                   :: freto,freq1,freq2

      ! Characteristic velocity
      freq1 = staco(1)*visac*e%npol4/(chale(2)*chale(2)) ! Stokes term (using h=vol^(1/3))
      freq2 = staco(2)*denac*vnorm*e%npol/chale(1)       ! Convective term (using h= h_stream) 
      freto = freq1 + freq2                              ! Total frequency

      ! Compute stability paramters according to the case.
      timom = 0.0_rp
      if(freto>zeroc) then
         timom = 1.0_rp/freto
      else
         timom = 1.0e12_rp
      end if
      
   end subroutine
   
   subroutine ComputeTauCDR(e,denac,visac,react,vnorm,staco,chale,timom)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: staco(3)
      real(rp),    intent(in)    :: denac,visac,react
      real(rp),    intent(in)    :: vnorm,chale(2)
      real(rp),    intent(out)   :: timom
      real(rp)                   :: freto

      !Stabilization paramter
      freto =  staco(1)*visac*e%npol4/(chale(2)*chale(2)) &
            + staco(2)*denac*vnorm*e%npol/chale(1)       & 
            + staco(3)*react
 
      ! Compute stability paramters according to the case.
      timom = 0.0_rp
      if(freto>zeroc) then
         timom = 1.0_rp/freto
      else
         timom = 1.0e12_rp
      end if

   end subroutine 
   
   subroutine ComputeTransientTau(dtinv,acden,timom_static,timom)
      implicit none
      real(rp) :: dtinv, timom_static, timom,acden
      
      timom = 1/(1/timom_static + acden*dtinv)
   end subroutine
   
   subroutine ComputeTauBoundary(e,chale,visac,temom) !per element
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: visac,chale(2)
      real(rp),    intent(out)   :: temom
      real(rp)                   :: freq

      ! Characteristic velocity
      freq = visac*e%npol/chale(1)      ! n.grad u 

      temom = 0.0_rp
      if(freq>zeroc) then
         temom = 1.0_rp/freq
      else
         temom = 1.0e-12_rp
      end if
      
   end subroutine

   subroutine elmvis(e,dvolt0,visac,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for viscosity for block U,V 
      !     mu*(grad v, grad u)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: visac,dvolt0
      real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)

      integer(ip)                :: inode,jnode
      real(rp)                   :: tmp

      tmp = dvolt0*visac

      do inode=1,e%pnode
         do jnode=1,e%pnode
            elmat(inode,jnode) = tmp*dot_product(e%cartd(:,inode),e%cartd(:,jnode)) &
                  + elmat(inode,jnode)
         end do
      end do

  end subroutine elmvis
  
  subroutine elmmas(e,dvol,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the elemental mass matrix
    !   (v, u/dt)
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none
    class(FiniteElement)        :: e
    real(rp),    intent(in)    :: dvol
    real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
    integer(ip)                :: inode,jnode
    real(rp)                   :: aux1

    do jnode=1,e%pnode
       aux1 = e%shape(jnode,e%igaus)*dvol
       do inode=1,e%pnode
          elmat(inode,jnode) = e%shape(inode,e%igaus)*aux1 + elmat(inode,jnode)
       end do
    end do

  end subroutine elmmas
  
  subroutine elmmas_Lumped(e,dvol,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the elemental mass matrix
    !   (v, u/dt)
    !
    !-----------------------------------------------------------------------
    use typre
    implicit none
    class(FiniteElement)        :: e
    real(rp),    intent(in)    :: dvol
    real(rp),    intent(inout) :: elmat(e%mnode,e%mnode)
    integer(ip)                :: inode,jnode
    real(rp)                   :: aux1

    do jnode=1,e%pnode
       aux1 = e%shape(jnode,e%igaus)*dvol
       do inode=1,e%pnode
          elmat(inode,inode) = e%shape(inode,e%igaus)*aux1 + elmat(inode,inode)
       end do
    end do

  end subroutine elmmas_Lumped  
  
  subroutine grad_proj(e,ndime,dvolu,elgr,elrh)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the projection into the Finite element space
    !  (v, var)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    integer(ip), intent(in)    :: ndime
    real(rp),    intent(in)    :: elgr(ndime,e%ndime)
    real(rp),    intent(inout) :: elrh(ndime,e%ndime,e%mnode)

    integer(ip)                :: idime,jdime,inode

    do idime=1,ndime
         do jdime = 1,e%ndime
            do inode = 1,e%pnode
               elrh(idime,jdime,inode) = elrh(idime,jdime,inode) +&
                      elgr(idime,jdime)*e%shape(inode,e%igaus)*dvolu
            enddo
         enddo
    enddo

   end subroutine grad_proj

 end module  
