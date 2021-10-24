subroutine frivel(delta,tveno,vikin,velfr)
!------------------------------------------------------------------------
!****f* Mathru/frivel
! NAME 
!    frivel
! DESCRIPTION
!    Computes the friction velocity by the Newton-Raphson method. This 
!    friction velocity is computed from U_p/U_*=(1/X)log(U_* D/nu)+C, 
!    where X=0.41 (VONKA) is the von Karman constant, U_p is the tangent 
!    wall velocity at the current integration point, C is a constant 
!    which value is 5.5 and D = DELTA is a user defined parameter 
!    (DELTA, its physical meaning being the boundary layer width. 
!    The result is VELFR.
! USES
! USED BY
!    nsm_bouwal
!    tur_updnbc
!***
!------------------------------------------------------------------------
  use      typre
  implicit none
  real(rp), intent(in)  :: delta,tveno,vikin
  real(rp), intent(out) :: velfr
  integer(ip)           :: itera,ndelt,j
  real(rp)              :: xmuit,fdvfr,devfr
  real(rp)              :: vkinv,diffd,parco,dplus
  real(rp)              :: ypele,expye,expyt,oneoe,firsl
!
! Wall law taking into account the buffer zone and log layer:
! u+ =1/k*ln(1+0.4*y+)+7.8*[1-exp(-y+/11-y+/11*exp(-0.33*y+)]   
!           
  if (delta>0.0_rp) then
     velfr=sqrt(tveno*vikin/delta)
     if(velfr*delta/vikin>5.0_rp) then
        vkinv=1.0_rp/0.41_rp
        xmuit=delta/vikin                                               ! D/nu
        itera=0                                                         ! i=0
        velfr=tveno                                                     ! U_*^0
        parco=1.0_rp
        oneoe=1.0_rp/11.0_rp
        do while((parco>=1.0e-6).and.itera<100)
           itera=itera+1                                                ! i=i+1
           if(velfr<0.0_rp) velfr=0.0_rp
           dplus=velfr*xmuit
           ypele=dplus*oneoe
           expye=exp(-ypele)
           expyt=exp(-dplus*0.33_rp) 
           firsl=vkinv*log(1.0_rp+0.4_rp*dplus)
           fdvfr=velfr*(firsl+7.8_rp*(1.0_rp-expye-ypele*expyt))-tveno  ! F(U_*^i)
           diffd=firsl+vkinv*0.4_rp*dplus/(1.0_rp+0.4_rp*dplus)&        ! DF(U_*^i)
                +7.8_rp*(1.0_rp-expye*(1.0_rp+ypele)&
                -ypele*expyt*(2.0_rp-dplus*0.33_rp))               
           devfr=-fdvfr/diffd                                           ! delta(U_*^0)
           parco=abs(devfr/tveno)                                       ! U_*^i=U_*^i-1
           velfr=velfr+devfr                                            ! +delta(U_*)
        end do
     end if
!
! If DELTA = 0, the friction velocity is set to zero.
!
  else
     velfr=0.0_rp
  end if
  
end subroutine frivel
