subroutine lmn_bouwal(a,e,gpvel,visac,denac,delta,&
                      wmatr,tract)

!-----------------------------------------------------------------------
!****f* Nstinc/lmn_bouwal
! NAME 
!    lmn_bouwal
! DESCRIPTION
!    This routine computes the surface traction for the NS equations at
!    a given integration point of a boundary IBOUN received by argument
!    due to the use of a turbulent wall law. The algorithm is:
!    - Compute the tangential velocity u at the integration point x.
!      In fact, u is not necessarily tangential at x. For example:
!      -> u    -> u      
!      o---x---o---x---o
!                      |\ u
!                      | 
!    - Given the distance to the wall y, compute (U*) (friction velocity)
!    - Compute y+=(y)(U*)/(nu)
!      if(y+>5) then
!        t=-(rho)*(U*)^2*(u/|u|)
!      else
!        u+=y+ => U*^2=u*nu/y so that
!        t=-mu*u/y
!      end if
! USES
!    vecnor
!    frivel
! USED BY
!    lmn_bouope
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_LowMach
   implicit none
   type(LowMachProblem) :: a
   class(FiniteElement)  :: e
   real(rp) :: gpvel(e%ndime),visac,denac,delta
   real(rp) :: wmatr(a%ndofn,e%mnode,a%ndofn,e%mnode),tract(e%ndime)

   integer(ip)                :: jlins
   integer(ip)                :: ldime,idime,ldimb

   real(rp)                   :: vikin,velfr,yplus    ! nu, U*, y+, y
   real(rp)                   :: tveno,tvelo(e%ndime),ovelo(e%ndime)             ! |u|, u
   real(rp)                   :: trano,xmuit,tcoef,dtrac
   integer(ip)                :: inodb,jnodb,jdime,kdime,inode,jnode
   real(rp)                   :: relax

   jlins = 1
   if(delta<zelmn) return

   !normal velocity, projection of a over b: a_b = (a . b) b
   ovelo = dot_product(e%baloc(:,e%ndime),gpvel)*e%baloc(:,e%ndime)
   
   tvelo = gpvel - ovelo
   call vecnor(tvelo,e%ndime,tveno,2)                       ! |u|
   vikin=visac/denac                                      ! nu
   if(tveno<=zelmn) then
      velfr=0.0_rp
      tract=0.0_rp
      return
   else
      call frivel(delta,tveno,vikin,velfr)                  ! U*
   end if

   ! Compute prescribed traction
   if(jlins==0) then                                       ! RHS
 
       yplus=delta*velfr/vikin
       if(yplus<5.0_rp) then
          tract = -denac*vikin*tvelo/delta                  ! t=-mu*u/y
       else
          tract = -denac*velfr*velfr*tvelo/tveno            ! t=-rho*U*^2*(u/|u|)
       end if

   else if(jlins==1) then                                  ! Assembly (the sign changes on the left)

      yplus=delta*velfr/vikin
      if(yplus<5.0_rp) then
         trano=denac*vikin/delta                              ! t=-mu*u/y
      else
         trano=denac*velfr*velfr/tveno                        ! t=-rho*U*^2*(u/|u|)
      end if
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do idime=1,e%ndime
            do jnodb=1,e%pnodb
               jnode = e%lboel(jnodb)
               xmuit=e%shapb(inodb,e%igaub)*e%shapb(jnodb,e%igaub)
               do jdime=1,e%ndime
                  do ldimb=1,e%ndime-1
                     wmatr(idime,inode,jdime,jnode) = &
                        wmatr(idime,inode,jdime,jnode) + &
                         xmuit*e%baloc(idime,ldimb)*e%baloc(jdime,ldimb)*trano
                  end do
               end do
            end do
         end do
      end do

   elseif(jlins==2) then                                   ! Newton on the traction

   else
      call runend('lmn_bouwal: jlins not defined')
   end if
 
end subroutine lmn_bouwal
