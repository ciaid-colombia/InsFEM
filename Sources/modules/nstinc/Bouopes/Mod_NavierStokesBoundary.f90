module Mod_NavierStokesBoundary
   use typre
   use Mod_Element
   
contains

   subroutine nsm_bopres(e,wmatr)
   !-----------------------------------------------------------------------
   !    This routine computes the contribution to WMATR for the Navier-
   !    Stokes equations due to integration on the boundary of v*p*n.
   !    It is the same idea as open boundaries but only the pressure
   !    part and it is done when there are walls.
   !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(inout) :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode)
      integer(ip) :: ldime,inodb,jnodb,inode,jnode
      real(rp)    :: prod1,xmuit
   
      do ldime=1,e%ndime
         prod1=e%baloc(ldime,e%ndime)
         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            xmuit=e%shapb(inodb,e%igaub)*prod1
            do jnodb=1,e%pnodb
               jnode = e%lboel(jnodb)
               wmatr(ldime,inode,e%ndime+1,jnode) = wmatr(ldime,inode,e%ndime+1,jnode)&
                  +xmuit*e%shapb(jnodb,e%igaub)
            end do
         end do
      end do
    
   end subroutine nsm_bopres

   subroutine nsm_bouopb(e,wmatr,acvis,fvins)
   !-----------------------------------------------------------------------
   !    This routine computes the contribution to WMATR for the Navier-
   !    Stokes equations due to open boundaries. In this case, the term
   !    S.n (S = -p I + 2 mu Sym(grad(u)) being the Cauchy stress tensor) 
   !    is considered unknown.
   !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(inout) :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode)
      real(rp), intent(in)    :: acvis,fvins
      real(rp), parameter     :: zensi = epsilon(1.0_rp)
      integer(ip) :: inodb,idime,jnode,kdime,inode,jnodb
      real(rp)    :: xmuit
      real(rp)    :: prod1
      
      ! Contribution from the viscous term: 2 mu Sym(grad(u).n
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         xmuit=e%shapb(inodb,e%igaub)*acvis
         do idime=1,e%ndime
            do jnode=1,e%pnode
               wmatr(idime,inode,idime,jnode) = wmatr(idime,inode,idime,jnode)  &
                -xmuit*dot_product(e%cartb(1:e%ndime,jnode),e%baloc(1:e%ndime,e%ndime))  
            end do
            if (fvins>zensi) then
               do jnode=1,e%pnode
                  do kdime=1,e%ndime
                     wmatr(idime,inode,kdime,jnode)=wmatr(idime,inode,kdime,jnode)&
                        -xmuit*e%cartb(idime,jnode)*e%baloc(kdime,e%ndime)
                  end do
               end do
            end if
         end do
      end do
      
      !Contribution from the pressure term. (v,n*p)_Gamma
      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb = 1,e%pnodb
            jnode = e%lboel(jnodb)
            do idime = 1,e%ndime
               wmatr(idime,inode,e%ndime+1,jnode) = wmatr(idime,inode,e%ndime+1,jnode) + &
                  e%shapb(inodb,e%igaub)*e%shapb(jnodb,e%igaub)*e%baloc(idime,e%ndime)
            enddo
         enddo
      enddo

   end subroutine nsm_bouopb
   
   subroutine nsm_bouopb_p(e,wmatr,acvis,fvins)
   !-----------------------------------------------------------------------
   !    This routine computes the contribution to WMATR for the Navier-
   !    Stokes equations due to open boundaries but prescribing the pressure 
   !    in a weak form in the boundary. In this case, only the term
   !    S.n (S =  2 mu Sym(grad(u)) being the Cauchy stress tensor) 
   !    is considered unknown.
   !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(inout) :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode)
      real(rp), intent(in)    :: fvins,acvis
      real(rp), parameter     :: zensi = epsilon(1.0_rp)
      integer(ip) :: inodb,ldime,jnode,kdime,inode,jnodb
      real(rp)    :: xmuit
      real(rp)    :: prod1

      ! Contribution from the viscous term: 2 mu Sym(grad(u).n
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         xmuit=e%shapb(inodb,e%igaub)*acvis
         do ldime=1,e%ndime
            do jnode=1,e%pnode
               wmatr(ldime,inode,ldime,jnode) = wmatr(ldime,inode,ldime,jnode)  &
                  -xmuit*dot_product(e%cartb(1:e%ndime,jnode),e%baloc(1:e%ndime,e%ndime))  
            end do
            if (fvins>zensi) then
               do jnode=1,e%pnode
                  do kdime=1,e%ndime
                     wmatr(ldime,inode,kdime,jnode)=wmatr(ldime,inode,kdime,jnode)&
                        -xmuit*e%cartb(ldime,jnode)*e%baloc(kdime,e%ndime)
                  end do
               end do
            end if
         end do
      end do   
     
   end subroutine nsm_bouopb_p

   subroutine nsm_bouwal(e,gpvel,visac,denac,delta,wmatr,tract,yplus)
   !-----------------------------------------------------------------------
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
   !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)    :: e
      real(rp), intent(in)    :: gpvel(e%ndime),visac,denac,delta
      real(rp), intent(inout) :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode),tract(e%ndime)
      real(rp),parameter      :: zensi = epsilon(1.0_rp)
      integer(ip) :: jlins 
      integer(ip) :: ldime,idime,ldimb
      real(rp)    :: vikin,velfr,yplus    ! nu, U*, y+, y
      real(rp)    :: tveno,tvelo(e%ndime),ovelo(e%ndime)             ! |u|, u
      real(rp)    :: trano,xmuit,tcoef,dtrac
      integer(ip) :: inodb,jnodb,jdime,kdime,inode,jnode
      real(rp)    :: relax
   
      jlins = 1
      if(delta<zensi) return
      !normal velocity, projection of a over b: a_b = (a . b) b
      ovelo = dot_product(e%baloc(:,e%ndime),gpvel)*e%baloc(:,e%ndime)
      tvelo = gpvel - ovelo
      call vecnor(tvelo,e%ndime,tveno,2)             ! |u|
      vikin=visac/denac                              ! nu
      if(tveno<=zensi) then
         velfr=0.0_rp
         tract=0.0_rp
         return
      else
         call frivel(delta,tveno,vikin,velfr)        ! U*
      end if
   
      ! Compute prescribed traction
      if(jlins==0) then                              ! RHS
          yplus=delta*velfr/vikin
          if(yplus<5.0_rp) then
             tract = -denac*vikin*tvelo/delta        ! t=-mu*u/y
          else
             tract = -denac*velfr*velfr*tvelo/tveno  ! t=-rho*U*^2*(u/|u|)
          end if
      else if(jlins==1) then                         ! Assembly (the sign changes on the left)
         yplus=delta*velfr/vikin
         if(yplus<5.0_rp) then
            trano=denac*vikin/delta                  ! t=-mu*u/y
         else
            trano=denac*velfr*velfr/tveno            ! t=-rho*U*^2*(u/|u|)
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
   
      elseif(jlins==2) then                          ! Newton on the traction
      else
         call runend('nsm_bouwal: jlins not defined')
      end if
       
   end subroutine nsm_bouwal


!----------------------------------------------------------------------------
!Boundary SGS

!Traction terms

   subroutine nsm_elmbuv_tr(e,fvins,dsurf,acvis,temom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,V 
      !  - tau1*(F*v, Fu)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: fvins,temom,dsurf,acvis
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb,idime,jdime,kdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            do idime=1,e%ndime
               elmat(idime,inode,idime,jnode) =  elmat(idime,inode,idime,jnode) - 0.25*acvis*acvis*temom*dsurf &
                  * (dot_product(e%cartb(1:e%ndime,jnodb),e%baloc(1:e%ndime,e%ndime)) &
                  * dot_product(e%cartb(1:e%ndime,inodb),e%baloc(1:e%ndime,e%ndime))) 
               if (fvins>zensi) then
                  do jdime=1,e%ndime
                     elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) &
                        -0.25*acvis*acvis*temom*(e%cartb(idime,jnodb)*e%baloc(jdime,e%ndime)* &
                        dot_product(e%cartb(1:e%ndime,inodb),e%baloc(1:e%ndime,e%ndime)) &
                        + e%cartb(idime,inodb)*e%baloc(jdime,e%ndime)* &
                        dot_product(e%cartb(1:e%ndime,jnodb),e%baloc(1:e%ndime,e%ndime)))
                     do kdime=1,e%ndime
                        elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) -0.25*acvis*acvis*temom* &
                           e%cartb(kdime,jnodb)*e%baloc(jdime,e%ndime)*e%cartb(kdime,inodb)*e%baloc(idime,e%ndime)
                     end do
                  end do
               end if
            end do
         end do
      end do
   end subroutine nsm_elmbuv_tr 

   subroutine nsm_elmbpv_tr(e,fvins,dsurf,acvis,temom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,V 
      !     + tau1*(F*v, F p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: fvins,dsurf,acvis,temom
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb,idime,jdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)

      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb = 1,e%pnodb
            jnode = e%lboel(jnodb)
            elmat(1:e%ndime,inode,1,jnode) = elmat(1:e%ndime,inode,1,jnode) + 0.5*dsurf*temom*acvis* &
               e%shapb(jnodb,e%igaub)* dot_product(e%cartb(1:e%ndime,inodb),e%baloc(1:e%ndime,e%ndime))
            if (fvins>zensi) then
               do idime = 1,e%ndime
                  do jdime = 1,e%ndime
                     elmat(idime,inode,1,jnode) = elmat(idime,inode,1,jnode) + 0.5*dsurf*temom*acvis* &
                        e%shapb(jnodb,e%igaub)*e%cartb(jdime,inodb)*e%baloc(idime,e%ndime)
                  enddo
               enddo
            end if
         enddo
      enddo
  end subroutine nsm_elmbpv_tr

  subroutine nsm_elmbuq_tr(e,fvins,dsurf,acvis,temom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q
      !     - tau*(F*q, Fu) = - tau*(q, grad(u)+gradT(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: fvins,dsurf,acvis,temom
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb,idime,jdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            elmat(1,inode,1:e%ndime,jnode) =  elmat(1,inode,1:e%ndime,jnode) - 0.5*temom*acvis*dsurf* &
               e%shapb(inodb,e%igaub)*dot_product(e%cartb(1:e%ndime,jnodb),e%baloc(1:e%ndime,e%ndime)) 
            if (fvins>zensi) then
               do idime = 1,e%ndime
                  do jdime = 1,e%ndime
                     elmat(1,inode,idime,jnode) = elmat(1,inode,idime,jnode) - 0.5*dsurf*temom*acvis* &
                        e%shapb(inodb,e%igaub)*e%cartb(jdime,jnodb)*e%baloc(idime,e%ndime)
                  enddo
               enddo
            end if
         end do
      end do
   end subroutine nsm_elmbuq_tr

   subroutine nsm_elmbpq_tr(e,dsurf,temom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,Q
      !    (q, p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dsurf,temom
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + temom*dsurf*&
               e%shapb(jnodb,e%igaub)*e%shapb(inodb,e%igaub)
         end do
      end do
   end subroutine nsm_elmbpq_tr

   subroutine nsm_elmrhu_tr(e,fvins,dsurf,acvis,temom,tract,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: tract(e%ndime)
      real(rp),    intent(in)    :: fvins,dsurf,temom,acvis
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,inodb,idime,jdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)

      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)
         elrhs(1:e%ndime,inode) = elrhs(1:e%ndime,inode) + 0.5*dsurf*temom*acvis* &
            tract(1:e%ndime)* dot_product(e%cartb(1:e%ndime,inodb),e%baloc(1:e%ndime,e%ndime))
         if (fvins>zensi) then
            do idime = 1,e%ndime
               do jdime = 1,e%ndime
                  elrhs(idime,inode) = elrhs(idime,inode) + 0.5*dsurf*temom*acvis &
                     *tract(jdime)*e%cartb(jdime,inodb)*e%baloc(idime,e%ndime)
               enddo
            enddo
         end if
      enddo
   end subroutine nsm_elmrhu_tr
   
   subroutine nsm_elmrhp_tr(e,dsurf,temom,tract,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    tau*(n.f,q)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: tract(e%ndime)
      real(rp),    intent(in)    :: temom,dsurf
      real(rp),    intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                :: inode,inodb

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         elrhs(1,inode) = elrhs(1,inode) + dsurf*temom* &
            e%shapb(inodb,e%igaub)*dot_product(tract,e%baloc(:,e%ndime))
      end do
   end subroutine nsm_elmrhp_tr

!Normal velocity terms (DG) (Nitsche)

   subroutine nsm_elmbuv_ni(e,dsurf,temom,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,V 
      !  - tau1*(F*v, Fu)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: temom,dsurf
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb,idime,jdime

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            do idime=1,e%ndime
               do jdime=1,e%ndime
                  elmat(idime,inode,jdime,jnode) =  elmat(idime,inode,jdime,jnode) - temom*dsurf &
                     *e%baloc(idime,e%ndime)*e%baloc(jdime,e%ndime)*e%shapb(inodb,e%igaub)*e%shapb(jnodb,e%igaub)
               end do
            end do
         end do
      end do
   end subroutine nsm_elmbuv_ni

   subroutine nsm_elmrhu_ni(e,dsurf,temom,gpvel,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: dsurf,temom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,inodb,idime,jdime

      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)
         elrhs(1:e%ndime,inode) = elrhs(1:e%ndime,inode) - dsurf*temom*e%shapb(inodb,e%igaub) &
            *e%baloc(1:e%ndime,e%ndime) * dot_product(gpvel(1:e%ndime),e%baloc(1:e%ndime,e%ndime))
      enddo
   end subroutine nsm_elmrhu_ni
   
!Consistency terms (DG)

   subroutine nsm_elmbuv_cn(e,fvins,dsurf,acvis,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,V 
      !  - tau1*(F*v, Fu)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: fvins,dsurf,acvis
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb,idime,jdime,kdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)
      real(rp)                   :: dx,dy,dz,n1,n2,n3
      real(rp)                   :: B(e%ndime,e%ndime)

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)

            B = 0.0_rp
            if (e%ndime .eq. 2) then

               dx=e%cartb(1,jnodb)
               dy=e%cartb(2,jnodb)
               n1=e%baloc(1,e%ndime)
               n2=e%baloc(2,e%ndime)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = e%shapb(inodb,e%igaub)*(n1*dx + n2*dy)
               B(2,2) = B(1,1)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = B(1,1) + e%shapb(inodb,e%igaub)*(n1*dx)
               B(1,2) = B(1,2) + e%shapb(inodb,e%igaub)*(n2*dx)
               B(2,1) = B(2,1) + e%shapb(inodb,e%igaub)*(n1*dy)
               B(2,2) = B(2,2) + e%shapb(inodb,e%igaub)*(n2*dy)

            elseif(e%ndime .eq. 3) then

               dx=e%cartb(1,jnodb)
               dy=e%cartb(2,jnodb)
               dz=e%cartb(3,jnodb)
               n1=e%baloc(1,e%ndime)
               n2=e%baloc(2,e%ndime)
               n3=e%baloc(3,e%ndime)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = e%shapb(inodb,e%igaub)*(n1*dx + n2*dy + n3*dz)
               B(2,2) = B(1,1)
               B(3,3) = B(1,1)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = B(1,1) + e%shapb(inodb,e%igaub)*(n1*dx)
               B(1,2) = B(1,2) + e%shapb(inodb,e%igaub)*(n2*dx)
               B(1,3) = B(1,3) + e%shapb(inodb,e%igaub)*(n3*dx)
               B(2,1) = B(2,1) + e%shapb(inodb,e%igaub)*(n1*dy)
               B(2,2) = B(2,2) + e%shapb(inodb,e%igaub)*(n2*dy)
               B(2,3) = B(2,3) + e%shapb(inodb,e%igaub)*(n3*dy)
               B(3,1) = B(3,1) + e%shapb(inodb,e%igaub)*(n1*dz)
               B(3,2) = B(3,2) + e%shapb(inodb,e%igaub)*(n2*dz)
               B(3,3) = B(3,3) + e%shapb(inodb,e%igaub)*(n3*dz)
            endif

            elmat(1:e%ndime,inode,1:e%ndime,jnode) = elmat(1:e%ndime,inode,1:e%ndime,jnode) - 0.25*acvis*dsurf*B(1:e%ndime,1:e%ndime)

         end do
      end do

      !do inodb=1,e%pnodb
      !   inode = e%lboel(inodb)
      !   do jnodb=1,e%pnodb
      !      jnode = e%lboel(jnodb)
      !      do idime=1,e%ndime
      !         elmat(idime,inode,idime,jnode) =  elmat(idime,inode,idime,jnode) - 0.5*acvis*dsurf* &
      !            (dot_product(e%cartb(1:e%ndime,jnodb),e%baloc(1:e%ndime,e%ndime)) * e%shapb(inodb,e%igaub)) !&
      !            !+e%shapb(jnodb,e%igaub)* dot_product(e%cartb(1:e%ndime,inodb),e%baloc(1:e%ndime,e%ndime))!)
      !         if (fvins>zensi) then
      !            do jdime=1,e%ndime
      !               elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) -0.5*acvis* &
      !                  *(e%shapb(inodb,e%igaub)*e%cartb(idime,jnodb)*e%baloc(jdime,e%ndime)) !&
      !         !         + e%shapb(jnodb,e%igaub)*e%cartb(jdime,inodb)*e%baloc(idime,e%ndime))
      !            end do
      !         end if
      !      end do
      !   end do
      !end do
   end subroutine nsm_elmbuv_cn

   subroutine nsm_elmbpv_cn(e,dsurf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block P,V 
      !     + tau1*(F*v, F p)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dsurf
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb

      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb = 1,e%pnodb
            jnode = e%lboel(jnodb)
            elmat(1:e%ndime,inode,1,jnode) = elmat(1:e%ndime,inode,1,jnode) + 0.5*dsurf* &
               e%shapb(jnodb,e%igaub)*e%shapb(inodb,e%igaub)*e%baloc(1:e%ndime,e%ndime)
         enddo
      enddo
   end subroutine nsm_elmbpv_cn

   subroutine nsm_elmbuv_acn(e,fvins,dsurf,acvis,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for OSS for block U,V 
      !  - tau1*(F*v, Fu)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: fvins,dsurf,acvis
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb,idime,jdime,kdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)
      real(rp)                   :: dx,dy,dz,n1,n2,n3
      real(rp)                   :: B(e%ndime,e%ndime)

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)

            B = 0.0_rp

            if(e%ndime .eq. 2) then

               dx=e%cartb(1,inodb)
               dy=e%cartb(2,inodb)
               n1=e%baloc(1,e%ndime)
               n2=e%baloc(2,e%ndime)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = e%shapb(jnodb,e%igaub)*(n1*dx + n2*dy)
               B(2,2) = B(1,1)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = B(1,1) + e%shapb(jnodb,e%igaub)*(n1*dx)
               B(2,1) = B(2,1) + e%shapb(jnodb,e%igaub)*(n2*dx)
               B(1,2) = B(1,2) + e%shapb(jnodb,e%igaub)*(n1*dy)
               B(2,2) = B(2,2) + e%shapb(jnodb,e%igaub)*(n2*dy)

            elseif(e%ndime .eq. 3) then

               dx=e%cartb(1,inodb)
               dy=e%cartb(2,inodb)
               dz=e%cartb(3,inodb)
               n1=e%baloc(1,e%ndime)
               n2=e%baloc(2,e%ndime)
               n3=e%baloc(3,e%ndime)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = e%shapb(jnodb,e%igaub)*(n1*dx + n2*dy + n3*dz)
               B(2,2) = B(1,1)
               B(3,3) = B(1,1)

               !Symmetric part 1/2(nj du_i/dxj,v_i)
               B(1,1) = B(1,1) + e%shapb(jnodb,e%igaub)*(n1*dx)
               B(2,1) = B(2,1) + e%shapb(jnodb,e%igaub)*(n2*dx)
               B(3,1) = B(3,1) + e%shapb(jnodb,e%igaub)*(n3*dx)
               B(1,2) = B(1,2) + e%shapb(jnodb,e%igaub)*(n1*dy)
               B(2,2) = B(2,2) + e%shapb(jnodb,e%igaub)*(n2*dy)
               B(3,2) = B(3,2) + e%shapb(jnodb,e%igaub)*(n3*dy)
               B(1,3) = B(1,3) + e%shapb(jnodb,e%igaub)*(n1*dz)
               B(2,3) = B(2,3) + e%shapb(jnodb,e%igaub)*(n2*dz)
               B(3,3) = B(3,3) + e%shapb(jnodb,e%igaub)*(n3*dz)
            endif

            elmat(1:e%ndime,inode,1:e%ndime,jnode) = elmat(1:e%ndime,inode,1:e%ndime,jnode) - 0.25*acvis*dsurf*B(1:e%ndime,1:e%ndime)

         end do
      end do
   end subroutine nsm_elmbuv_acn

  subroutine nsm_elmbuq_acn(e,dsurf,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,Q
      !     - tau*(F*q, Fu) = - tau*(q, grad(u)+gradT(u))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement)       :: e
      real(rp),    intent(in)    :: dsurf
      real(rp),    intent(inout) :: elmat(1,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,inodb,jnodb

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            elmat(1,inode,1:e%ndime,jnode) =  elmat(1,inode,1:e%ndime,jnode) - 0.5*dsurf* &
               e%shapb(inodb,e%igaub)*e%shapb(jnodb,e%igaub)*e%baloc(1:e%ndime,e%ndime)
         end do
      end do
   end subroutine nsm_elmbuq_acn

   subroutine nsm_elmrhu_cn(e,dsurf,trac,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: trac(e%ndime)
      real(rp),    intent(in)    :: dsurf
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,inodb,idime,jdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)

      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)
         do idime = 1,e%ndime
            elrhs(idime,inode) = elrhs(idime,inode) - 0.5*dsurf* &
               trac(idime)*e%shapb(inodb,e%igaub)
         end do
      enddo
   end subroutine nsm_elmrhu_cn
  
   subroutine nsm_elmrhu_acn(e,fvins,dsurf,acvis,gpvel,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    (v, f) + tau1*(L*v, f) + (v, u_n/dt) + tau1*(L*v, u_n/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: fvins,dsurf,acvis
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,inodb,idime,jdime
      real(rp),    parameter     :: zensi = epsilon(1.0_rp)
      real(rp)                   :: dx,dy,dz,n1,n2,n3
      real(rp)                   :: B(e%ndime)

      do inodb = 1,e%pnodb
         inode = e%lboel(inodb)

         B = 0.0_rp

         if(e%ndime .eq. 2) then

            dx=e%cartb(1,inodb)
            dy=e%cartb(2,inodb)
            n1=-1.0_rp*e%baloc(1,e%ndime)
            n2=-1.0_rp*e%baloc(2,e%ndime)

            !Symmetric part 1/2(nj du_i/dxj,v_i)
            B(1) = (n1*dx + n2*dy)*gpvel(1)
            B(2) = (n1*dx + n2*dy)*gpvel(2)

            !Symmetric part 1/2(nj du_i/dxj,v_i)
            B(1) = B(1) + (n1*dx + n1*dy)*gpvel(1)
            B(2) = B(2) + (n2*dx + n2*dy)*gpvel(2)

         elseif(e%ndime .eq. 3) then

            dx=e%cartb(1,inodb)
            dy=e%cartb(2,inodb)
            dz=e%cartb(3,inodb)
            n1=-1.0_rp*e%baloc(1,e%ndime)
            n2=-1.0_rp*e%baloc(2,e%ndime)
            n3=-1.0_rp*e%baloc(3,e%ndime)

            !Symmetric part 1/2(nj du_i/dxj,v_i)
            B(1) = (n1*dx + n2*dy + n3*dz)*gpvel(1)
            B(2) = (n1*dx + n2*dy + n3*dz)*gpvel(2)
            B(3) = (n1*dx + n2*dy + n3*dz)*gpvel(3)

            !Symmetric part 1/2(nj du_i/dxj,v_i)
            B(1) = B(1) + (n1*dx + n1*dy + n1*dz)*gpvel(1)
            B(2) = B(2) + (n2*dx + n2*dy + n2*dz)*gpvel(2)
            B(3) = B(3) + (n3*dx + n3*dy + n3*dz)*gpvel(3)
         endif
          
         elrhs(:,inode) = elrhs(:,inode) + 0.25*dsurf*acvis*B(:)

         !do idime = 1,e%ndime
         !   elrhs(idime,inode) = elrhs(idime,inode) + 0.25*dsurf* &
         !      (acvis*grvel(idime)*e%shapb(inodb,e%igaub)) !- gpvel(idime)* &
         !   !   dot_product(e%cartb(1:e%ndime,inodb),e%baloc(1:e%ndime,e%ndime)))
         !   !if (fvins>zensi) then
         !   !   do jdime = 1,e%ndime
         !   !      elrhs(idime,inode) = elrhs(idime,inode) - 0.5*dsurf &
         !   !         *gpvel(jdime)*e%cartb(jdime,inodb)*e%baloc(idime,e%ndime)
         !   !   enddo
         !   !end if
         !end do
      enddo
   end subroutine nsm_elmrhu_acn
  
   subroutine nsm_elmrhp_acn(e,dsurf,gpvel,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for ASGS & OSS
      !    tau*(n.f,q)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: dsurf
      real(rp),    intent(inout) :: elrhs(1,e%pnode)
      integer(ip)                :: inode,inodb

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         elrhs(1,inode) = elrhs(1,inode) + 0.5*dsurf*e%shapb(inodb,e%igaub)* &
            dot_product(gpvel,-1.0_rp*e%baloc(:,e%ndime))
      end do
   end subroutine nsm_elmrhp_acn

end module
