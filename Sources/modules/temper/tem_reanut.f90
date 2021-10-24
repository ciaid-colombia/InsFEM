subroutine tem_reanut(a,itask)
!    This routine reads the numerical treatment for the incompressible NS
!    equations.
   use typre
   use Mod_Listen
   use Mod_Temperature
   implicit none
   
   class (TemperatureProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   integer(ip) :: istab
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   
   !Initializations
   if (itask == 0) then
      a%kfl_repro = 0                                ! Res. projection not used 
      a%kfl_adapsgs = 0                              ! 1: Scaled L_2 norm, 2: Entropy function. Default is L_2 norm
      a%kfl_shock = 0                                ! a%shock capturing off
      a%kfl_wtemp = 1                                ! dT/dt is in the residual
      a%kfl_wlapl = 1                                ! -Lap(T) is in the residual
      a%kfl_trasg = 0                                ! Tracking of subscales
      a%kfl_tacsg = 0                                ! Time accuracy of subscales
      a%kfl_nolsg = 0                                ! Non-linearity of the subscales
      a%kfl_stabm = 1                                ! Default is Variational Multiscale
      a%staco(1)  = 4.0_rp                           ! Diffusive term
      a%staco(2)  = 2.0_rp                           ! Convective term
      a%staco(3)  = 1.0_rp                           ! Reactive term
      a%shock     = 0.0_rp                           ! SC parameter
      a%relax     = 1.0_rp                           ! Relaxation factor
   
   
   elseif (itask == 1) then
   
      if(words(1)=='STABI') then
         do istab = 1,3
            a%staco(istab) = param(istab)
         end do
      else if(words(1)=='TYPEO') then
         if(words(2)=='RESID') a%kfl_repro = 1
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
            !Residual based istropic
            a%kfl_shock = 1
            a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#a%shock capturing parameter')
         else if(a%Listener%exists('ANISO')) then
            !Residusl based anisotropic
            a%kfl_shock = 2
            a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#a%shock capturing parameter')
         elseif(a%Listener%exists('GISOT').or.a%Listener%exists('ON   ')) then
            !Gradient based isotropic
            a%kfl_shock = 3
            a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#a%shock capturing parameter')
         elseif(a%Listener%exists('GANIS')) then
            !Gradient based anisotropic
            a%kfl_shock = 4
            a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#a%shock capturing parameter')   
         endif
      else if(words(1)=='TEMPO') then
         if(a%Listener%exists('OFF  ')) a%kfl_wtemp = 0
      else if(words(1)=='DIFUS') then
         if(a%Listener%exists('OFF  ')) a%kfl_wlapl = 0
      else if(words(1)=='TRACK') then
         if(a%Listener%exists('ON   ')) then
            a%kfl_trasg=1
         end if
         if(a%kfl_trasg/=0) then
            a%relsg = a%Listener%getrea('RELAX',1.0_rp,'#Tracking relaxation')
            a%kfl_tacsg = a%Listener%getint('TIMEA',1,'#Tracking accuracy')
            if(a%Listener%exists('NONLI')) a%kfl_nolsg = 1
            if (a%kfl_tacsg > 1) call runend('Temperature: Dynamic subscales: only BDF1 implemented for the subscales')
            if (a%relsg /= 1.0_rp) call runend('Temperature: Dynamic subscales: relaxation not implemented')
         end if
      elseif (words(1) == 'ERRSU') then
         if (a%Listener%exists('SCALE')) then
            a%kfl_adapsgs = 1
         endif
      endif
      
   elseif (itask == 100) then !Final operations
   
      if (a%kfl_tacsg /= 0) then
         if (a%kfl_tacsg /= 1) then
            call runend('Temperature: only bdf1 for dynamics subgrid scaels')
         endif
         
         !We need to track the subscales
         a%kfl_trasg = 1
      endif
      if (a%kfl_adapsgs == 1) a%kfl_trasg = 1
   
   endif

   
end subroutine   
