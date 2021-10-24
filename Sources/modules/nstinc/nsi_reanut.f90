subroutine nsi_reanut(a,itask)
!-----------------------------------------------------------------------
!****f* Nstinc/nsi_reanut
! NAME 
!    nsi_reanut
! DESCRIPTION
!    This routine reads the numerical treatment for the incompressible NS
!    equations.
! USES
!    a%Listener%listen
! USED BY
!    nsi_turnon
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes
   implicit none
   
   class (NavierStokesProblem) :: a
   integer(ip) :: itask
   integer(ip) :: istab,ipart,rpalo

   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   

   call a%Listener%getarrs(words,param,nnpar,nnwor)
  
   if (itask == 0) then         !  Initializations (defaults)

      a%staco(1)  = 4.0_rp                               ! Viscous part of tau1
      a%staco(2)  = 2.0_rp                               ! Convective part of tau1
      a%staco(3)  = 1.0_rp                               ! Rotation part of tau1
      a%staco(4)  = 1.0_rp                               ! Coefficient of tau2
      a%kfl_repro = 0                                    ! Res. projection not used
!       a%kfl_penalty = 0                                ! Penalty form
      a%kfl_repro_SkipFE = 0                             ! Default is do not skip the Finite Element Part
      a%kfl_adapsgs = 0                                  ! 1: Scaled L_2 norm, 2: Entropy function. Default is L_2 norm
      a%kfl_shock = 0                                    ! Shock capturing off
      a%shock     = 0.0_rp                               ! SC parameter
      a%kfl_wtemp = 1                                    ! dT/dt is in the residual
      a%kfl_wlapl = 1                                    ! nu Lapl(u) is in the 
      a%kfl_hdifumin = 0                                 ! Use minimum h for the diffusive stabilization term

      a%kfl_trasg = 0                                    ! Tracking of subscales
      a%mtrit     = 1                                    ! Tracking iterations
      a%kfl_tacsg = 0                                    ! Time accuracy of subscales
      a%kfl_nolsg = 0                                    ! Non-linear subscales

      a%kfl_bousg = 0                                    ! Boundary subscales
      a%bstco(1)  = 10.0_rp                              ! Coefficient for boundary SGS n.gradu
      a%bstco(2)  = 10.0_rp                              ! Coefficient for Nitsche

      a%tosgs = 1.0_rp                                   ! Tracking tolerance
      a%relsg = 1.0_rp                                   ! Tracking relaxation
      a%kfl_nolsgNewtonRaphson = 0                       ! Non-linear subscales, newton-Raphson scheme
      a%kfl_stabm = 1                                    ! Default is Variational Multiscale
      a%kfl_tausm = 0                                    ! Default is do not smooth tau
      a%subrelax = 1                                     ! Relaxation parameter
      a%kfl_doAitken= .false.                            ! Aitken
      a%relax = 1.0_rp                           ! Aitken parameter domain decomposition
      a%relax_max = 1.0_rp                       ! Max aitken parameter
      a%kfl_PressTempSubscale = 0                        !Temporal residual pressure term

   
   elseif (itask == 1) then     !Inside the Numerical Treatment Block   
      if(words(1)=='STABI') then
        do istab = 1,4
           a%staco(istab) = param(istab)
        end do
      
      else if(words(1).eq.'HDIFU') then
         if(a%Listener%exists('MINIM')) then    
            a%kfl_hdifumin = 1
         endif
        
      else if(words(1).eq.'TAUSM') then
         if(a%Listener%exists('ON   ')) then
             a%kfl_tausm = max(1,a%Listener%getint('ON   ',1,'#Number of tau smoothing cycles'))
         end if
        
     else if(words(1).eq.'PENAL') then
         if(a%Listener%exists('ON   ')) a%kfl_penal = 1
         if(a%Listener%exists('CLASS')) a%kfl_penal = 1
         if(a%Listener%exists('ITERA')) a%kfl_penal = 2
         if(a%kfl_penal>=1)& 
              a%penal=param(2)

     else if(words(1)=='TYPEO') then
        if(words(2)=='RESID') then
            a%kfl_repro = 1 
            if(words(3) == 'SKIPF') then
               a%kfl_repro_SkipFE = 1
            elseif(words(3) == 'NOSKI') then
               a%kfl_repro_SkipFE = 0
            endif
        elseif (words(2)=='SPLIT') then
            a%kfl_repro = 2
            if(a%Listener%exists('OSS  ')) then 
               a%kfl_repro = 2 
            elseif (a%Listener%exists('OSSFO')) then
               a%kfl_repro = 3
            else
               !call runend('Split option not implemented')
            endif
        endif

     else if(words(1)=='TESTF') then
        if(words(2)=='VMS  ') then 
          a%kfl_stabm = 1
        elseif(words(2)=='SUPG ') then  
          a%kfl_stabm = 0
        elseif(words(2) == 'GLS  ') then
          a%kfl_stabm = -1
        endif
      
     else if(words(1)=='SHOCK') then
        if(a%Listener%exists('ISOTR').or.a%Listener%exists('ON   ')) then
           a%kfl_shock = 1
           a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#Shock capturing parameter')
        else if(a%Listener%exists('ANISO')) then
           a%kfl_shock = 2
           a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#Shock capturing parameter')
        end if

     else if(words(1)=='TEMPO') then
        if(a%Listener%exists('OFF  ')) a%kfl_wtemp = 0

     else if(words(1)=='VISCO') then
        if(a%Listener%exists('OFF  ')) a%kfl_wlapl = 0

     
     else if(words(1)=='TRACK') then
        if (words(2) == 'ON   ') then
            a%kfl_trasg = 1
        else
           a%kfl_trasg = 0
        end if
     else if(words(1)=='BOUND') then
        if (words(2) == 'ON   ') then
           a%kfl_bousg = 1
           do istab = 2,3
              a%bstco(istab-1) = param(istab)
           end do
        else
           a%kfl_bousg = 0
        end if
     elseif (words(1) == 'ERRSU') then
        if (a%Listener%exists('SCALE')) then
           a%kfl_adapsgs = 1
        endif
     else if(words(1)=='NONLI') then
        if(a%Listener%exists('LINEA')) then
           a%kfl_nolsg=0
        else
           a%kfl_nolsg=1
           a%mtrit = a%Listener%getint('ITERA',20,'#Tracking iterations')
           a%tosgs = a%Listener%getrea('TOLER',1.0e-8_rp,'#Tracking tolerance')
           if (a%Listener%exists('NEWTO')) a%kfl_nolsgNewtonRaphson = 1
        end if
     else if(words(1)=='RELAX') then
        a%relsg = a%Listener%getrea('RELAX',1.0_rp,'#Tracking relaxation')
     else if(words(1)=='ACCUR') then
        if(a%Listener%exists('PRESS')) then
         a%kfl_tacsg=1
         a%kfl_PressTempSubscale = 1
        end if  
        if(a%Listener%exists('DYNAM')) a%kfl_tacsg=1
        if(a%Listener%exists('QUASI')) a%kfl_tacsg=0

     else if(words(1) == 'SUBRE')then
         if(words(2) == 'ON ') a%subrelax = param(2)
         
     elseif(words(1) == 'FSSTA') then
        if (words(2) == 'ON   ') a%kfl_StabilizeFreeSurface = 1
        a%StabilizeFreeSurface_Param = param(2)

     elseif(words(1)=='AITKE') then
        if (a%Listener%exists('ON   ')) then 
           a%kfl_doAitken= .true.
        end if
        a%relax = param(2)
        a%relax_max = a%relax
     elseif(words(1) == 'ROBIN') then
        if (words(2) == 'ON   ') then 
            a%doRobin= .true.
           !            =    (1/dt)*rho_s*H_s       +  dt*beta_s see [Badia2008]
            a%alfa_robin= a%dtinv*param(2)*param(3) + (1.0/a%dtinv)*param(4)
        endif 
     endif
     
   elseif (itask == 100) then   !Finalize reading operations  
      
      if (a%kfl_tacsg > 0) a%kfl_trasg = 1
      if (a%kfl_nolsg > 0) a%kfl_trasg = 1
      if (a%kfl_adapsgs == 1) a%kfl_trasg = 1
      if (a%kfl_doAitken) a%kfl_docoupconv = .true.
   endif
      
end subroutine nsi_reanut
