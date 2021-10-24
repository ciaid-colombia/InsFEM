subroutine sup_reanut(a,itask)
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
   use Mod_Listen
   use Mod_ThreeField
   implicit none
   
   class (ThreeFieldNSProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: istab,ipart,rpalo,auxit

   real(rp)  :: auxre,auxre2
   !For Listener
   real(rp), pointer     :: param(:)
   character(5), pointer :: words(:)
   integer(ip), pointer  :: nnpar,nnwor
   

   call a%Listener%getarrs(words,param,nnpar,nnwor)
  
   if (itask == 0) then                                  ! Initializations (defaults)

      a%staco(1)  = 4.0_rp                               ! Viscous part of tau1
      a%staco(2)  = 2.0_rp                               ! Convective part of tau1
      a%staco(3)  = 0.25_rp                              ! Rotation part of tau1
      a%staco(4)  = 0.25_rp                              ! Coefficient of tau2 
      a%kfl_repro = 0                                    ! Res. projection not used
      a%kfl_reproBoundzero=1
      a%kfl_penal = 0                                    ! Penalty form
      a%kfl_repro_SkipFE = 0                             ! Default is do not skip the Finite Element Part
      a%kfl_shock = 0                                    ! Shock capturing off
      a%shock     = 0.0_rp                               ! SC parameter
      a%kfl_wtemp = 1                                    ! dT/dt is in the residual
      a%kfl_wlapl = 1                                    ! nu Lapl(u) is in the residual
      a%kfl_tiacc = 1                                    ! First order time integ.
      
      a%kfl_hdifumin=0
      a%kfl_tausm=0                                      ! Smooth taus
      a%kfl_ntausmooth=1                                 ! Number of consecutive taus smooths
      
      a%kfl_trasg = 0                                    ! Tracking of subscales
      a%mtrit     = 1                                    ! Tracking iterations
      a%kfl_tacsg = 0                                    ! Time accuracy of subscales
      a%kfl_nolsg = 0                                    ! Non-linear subscales
      a%tosgs = 1.0_rp                                   ! Tracking tolerance
      a%relsg = 1.0_rp                                   ! Tracking relaxation
      a%kfl_nolsgNewtonRaphson = 0                       ! Non-linear subscales, newton-Raphson scheme
      a%kfl_stabm = 1                                    ! Default is Variational Multiscale
      a%subrelax  = 1.0_rp                               ! Relaxation Parameter
      a%incremental = 60_ip                              ! Incremental parameter 
      a%penal = 6
      a%kfl_splitOSSMomentum = 1                         ! Split-OSS in momentum equation. 1=Off, 0=On
      a%kfl_splitOSSConstitutive = 1                     ! Split-OSS in constitutive equation. 1=Off, 0=On
      a%kfl_linearConstitutiveTerms = 1                  ! Linearization terms in constitutive equation. 0=Picard, 1=Newton
      a%kfl_linearConvectiveTerm = 1                     ! Linearization terms of convective term.  0=Picard, 1=Newton
      a%kfl_reproTemporalTermZero=1                      ! Temporal term in residual. 1=Off, 0=On
   
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
             a%kfl_tausm = 1
             a%kfl_ntausmooth=param(2)
         end if
         
      else if(words(1) == 'PENAL') then
         if(words(2)=='CLASS') a%kfl_penal = 1
            a%penal=param(2)

      else if(words(1)=='TYPEO') then
         if(words(2)=='RESID') then
            a%kfl_repro = 1 !OSS
            if(words(3) == 'SKIPF') then
               a%kfl_repro_SkipFE = 1
            elseif(words(3) == 'NOSKI') then
               a%kfl_repro_SkipFE = 1
            endif
         endif
        
         if(words(2)=='SPLIT'.and. words(3) == 'OSS')then 
            a%kfl_repro = 2 !Split OSS
            a%kfl_splitOSSMomentum = 0
            !a%kfl_splitOSSConstitutive=0
            
            if (a%MatProp(1)%lawvi>=0) then
               a%kfl_repro=4 
               a%kfl_splitOSSConstitutive=0
            end if   
         end if
        
         if(words(2)=='SPLIT' .and. words(3) == 'TOTAL')then
            a%kfl_repro = 3 !Split ASGS
            a%kfl_splitOSSMomentum = 0
         end if

      else if(words(1)=='BOUND' .and. words(2)=='OFF') then
         a%kfl_reproBoundzero=0
          write(*,*) 'Projection on the boundaries is available now.'

      else if(words(1)=='TESTF') then
         if(words(2)=='VMS  ') then 
            a%kfl_stabm = 1
         elseif(words(2)=='SUPG ') then  
            a%kfl_stabm = 0
         elseif(words(2) == 'GLS  ') then
            a%kfl_stabm = -1
         endif
      
      else if(words(1)=='SHOCK') then
         if(words(2) == 'ON ') then
            a%kfl_shock = 1
            a%shock(1:2) = param(2:3)            
         else if(words(2) == 'GRADU  ') then
            a%kfl_shock = 2
            a%shock(1:2) = param(2:3)                
         else if(words(2) == 'DEVSS  ') then
            a%kfl_shock = 3
            a%shock(1:2) = param(2:3)                
         else if(words(2) == 'ALL') then
            a%kfl_shock = 4
            a%shock(1:2) = param(2:3)              
         end if
     
      else if(words(1) == 'SUBRE')then
         if(words(2) == 'ON ') a%subrelax = param(2)
     
         
      else if(words(1) == 'INCRE')then
         if(words(2) == 'ON ') a%incremental = param(2)
         
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
         if(a%Listener%exists('DYNAM')) then
            a%kfl_tacsg=1
            if (a%kfl_repro==1 .or. a%kfl_repro==2) a%kfl_reproTemporalTermZero=0 !Switch off temporary term at the residual.
         end if 
         
         if(a%Listener%exists('QUASI')) a%kfl_tacsg=0
        
      else if(words(1)=='LINEA') then
         if(words(2)=='NEWTO') then
            a%kfl_linearConvectiveTerm=1
            a%kfl_linearConstitutiveTerms=1
         else
            a%kfl_linearConvectiveTerm=1
            a%kfl_linearConstitutiveTerms=0 
         end if
         
       else if(words(1)=='LINCO') then
         if(words(2)=='NEWTO') then
            a%kfl_linearConstitutiveTerms=1
         else
            a%kfl_linearConstitutiveTerms=0 
         end if   
      endif 

      
   elseif (itask == 100) then   !Finalize reading operations   
      if (a%kfl_tacsg > 0) a%kfl_trasg = 1
      if (a%kfl_nolsg > 0) a%kfl_trasg = 1
      
      if (a%kfl_advec==0) a%kfl_linearConvectiveTerm=0 
      
   endif 
   
   
end subroutine sup_reanut
