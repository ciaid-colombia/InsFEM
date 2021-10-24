function funcre(param,npara,ifuge,timev)
!------------------------------------------------------------------------
!
! This function yields a parabolic or periodic evolution
!
!------------------------------------------------------------------------
  use      typre
  use      def_parame
  implicit none
  real(rp)                 :: funcre
  integer(ip), intent(in)  :: npara,ifuge
  real(rp),    intent(in)  :: param(npara),timev
  integer(ip)              :: ipara
  real(rp)                 :: timea,timeb,funca,funcb,zerom,timec
  real(rp)                 :: timei,timef,Tgp,T0,Tp,Tn
  
  funcre=0.0_rp
  zerom=epsilon(1.0_rp)

  if(ifuge==0) then
!
! No time dependence 
!
     funcre=1.0_rp

  else if(ifuge==1) then
!
! Parabolic evolution
!
     if(param(1)-zerom<=timev.and.timev<=param(2)+zerom) then 
        funcre=param(3)*timev*timev+param(4)*timev+param(5)
     else if (timev>param(2)+zerom) then
        timea=param(2)
        funcre=param(3)*timea*timea+param(4)*timea+param(5)
     else if (timev<param(1)-zerom) then
        timea=param(1)
        funcre=param(3)*timea*timea+param(4)*timea+param(5)
     end if
  else if(ifuge==2) then

     !Periodic evolution
     !param(1) : tini
     !param(2) : tfin
     !param(3) : amplitude
     !param(4) : frequency
     !param(5) : phase
     !param(6) : mean value
     
     if(param(1)-zerom<=timev.and.timev<=param(2)+zerom) then 
        funcre=param(3)*sin(param(4)*timev+param(5))+param(6)
     else if (timev>param(2)+zerom) then
        timea=param(2)
        funcre=param(3)*sin(param(4)*timea+param(5))+param(6)
     else if (timev<param(1)-zerom) then
        timea=param(1)
        funcre=param(3)*sin(param(4)*timea+param(5))+param(6)
     end if
  else if(ifuge==3) then
!
! Discrete evolution
!
     timei=param(1)
     timef=param((npara/2-1)*2+1)

     if(timev<=timei) then
        funcre=param(2)
     else
        if(timev>=timef) then            ! Look for the time inside the period
           timec=timev
           do while(timec>timef)
              timec=timec-(timef-timei)
           end do
        else
           timec=timev
        end if
        ipara=0
        do while(ipara<npara/2)
           ipara=ipara+1
           if(timec<param((ipara-1)*2+1)) then
              timea=param((ipara-1)*2+1)
              funca=param((ipara-1)*2+2)
              timeb=param((ipara-2)*2+1)
              funcb=param((ipara-2)*2+2)
              funcre=(funcb-funca)/(timeb-timea)*(timec-timea)+funca
              ipara=npara/2
           end if
        end do
     end if

  else if(ifuge==4) then
!
! Special function to change boundary values
!
     funcre=param(2)

  else if(ifuge==5) then
!
! Parabolic discontinuous
!
     if(param(1)-zerom<=timev.and.timev<=param(2)+zerom) then 
        funcre=param(5)*timev*timev+param(4)*timev+param(3)
     else if (timev>param(2)+zerom) then
        funcre=param(6)
     else if (timev<param(1)-zerom) then
        funcre=param(7)
     end if
  else if(ifuge==6) then
! Time-Stretched pulse, as appeared in takemoto et al Acoustic analysis of the vocal tract during vowel production by finite-difference time-domain method
    Tgp = 0.646/param(1)
    funcre = exp(-(((timev-Tgp)/(0.29*Tgp))**2))
    
  else if(ifuge==7) then
! Rosenberg glottal model
  ! param(1) t0 initial time
  ! param(2) tf final time
  ! param(3) A  amplitude factor
  ! param(4) F0 fundamental frequency
  !  fs sampling frequency computed from time step
    if(param(1)-zerom<=timev.and.timev<=param(2)+zerom) then 
        timea = timev
    else if (timev>param(2)+zerom) then
        timea = param(2)
    else if (timev<param(1)-zerom) then
        timea = param(1)
    end if
    T0 = 1.0/param(4)
    Tp = 0.4*T0
    Tn = 0.16*T0
    timea = mod(timea,T0)
    if (timea<=Tp) then
        funcre = param(3)*0.5*(1.0-cos(pi*timea/Tp))
    else if(timea<=Tp+Tn) then
        funcre = param(3)*cos(pi*(timea-Tp)/(2.0*Tn))
    else
        funcre = 0.0
    end if

    else if(ifuge==8) then
    ! Constant velocity
      funcre = param(1)*timev + param(2)

    else if(ifuge==9) then
    ! Velocity ramp
      if(timev<=param(1)+zerom) then
         funcre = param(2)*(timev**2.0_rp) + param(3)*timev + param(4)
      else
         !funcre = (2*param(2)*param(1) + param(3))*(timev-param(1)) + (param(2)*(param(1)**2.0_rp) + param(3)*param(1) + param(4))
         funcre = param(2)*(param(1)**2.0_rp) + param(3)*param(1) + param(4)
      endif
      
    else if(ifuge==10) then
      !Hyperbolic
      funcre = 0.5_rp*(1.0_rp+tanh(8*(timev-param(1))))
    
    else if(ifuge==11) then
      !Sinusoidal function (Taylor 1998)
      funcre = 1.0_rp + sin(2.0_rp*4.0_rp*atan(1.0_rp)*timev/param(1))
    else if(ifuge==12) then
      !Sinusoidal function (Wall2003)
      funcre = 1.0_rp - cos(2.0_rp*4.0_rp*atan(1.0_rp)*timev/param(1))
      !funcre = cos(2.0_rp*4.0_rp*atan(1.0_rp)*timev/param(1))
    
  end if
  
end function funcre
