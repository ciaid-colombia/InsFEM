subroutine lmn_reanut(a,itask)
   use typre
   use Mod_LowMach
   implicit none
   class (LowMachProblem) :: a
   integer(ip) :: itask
   integer(ip) :: istab
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()

   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !  Initializations (defaults)

      a%staco(1)  = 4.0_rp      ! Viscous part of tau1
      a%staco(2)  = 2.0_rp      ! Convective part of tau1
      a%staco(3)  = 1.0_rp      ! Rotation part of tau1
      a%staco(4)  = 1.0_rp      ! Coefficient of tau2
      a%kfl_repro = 0           ! Res. projection not used
      a%kfl_tiacc = 1           ! First order time integ.
      a%kfl_adapsgs = 0         ! 1: Scaled L_2 norm, 2: Entropy function. Default is L_2 norm

      a%kfl_trasg = 0           ! Tracking of subscales
      a%mtrit     = 1           ! Tracking iterations
      a%kfl_tacsg = 0           ! Time accuracy of subscales
      a%kfl_nolsg = 0           ! Non-linear subscales
      a%kfl_nolsgScheme = 0     ! Linearization model subscales
      a%tosgs = 1.0_rp          ! Tracking tolerance
      a%relsg = 1.0_rp          ! Tracking relaxation
      a%kfl_stabm = 1           ! Default is Variational Multiscale
      a%epspe = 0.0001_rp       ! Penalization in confined flow 
      

   
   elseif (itask == 1) then     !Inside the Numerical Treatment Block   
      if(words(1)=='STABI') then
        do istab = 1,4
           a%staco(istab) = param(istab)
        end do

      else if(words(1)=='TYPEO') then
         if(words(2)=='RESID') then
             a%kfl_repro = 1 
         endif

      else if(words(1)=='TESTF') then
         if(words(2)=='VMS  ') then 
            a%kfl_stabm = 1
         elseif(words(2)=='SUPG ') then  
            a%kfl_stabm = 0
         elseif(words(2) == 'GLS  ') then
            a%kfl_stabm = -1
         endif
      
      else if(words(1)=='EPRES') then
         a%epspe = a%Listener%getrea('EPRES',0.0001_rp,'Pressure penalization')
      else if(words(1)=='TRACK') then
         if (words(2) == 'ON   ') then
            a%kfl_trasg = 1
         else
            a%kfl_trasg = 0
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
            if(a%Listener%exists('PICAR')) a%kfl_nolsgScheme = 1
            if(a%Listener%exists('GAUSS')) a%kfl_nolsgScheme = 2
            if(a%Listener%exists('NEWTO')) a%kfl_nolsgScheme = 0
            if(a%Listener%exists('SEGRE')) a%kfl_nolsgScheme = 3
         end if
      else if(words(1)=='RELAX') then
         a%relsg = a%Listener%getrea('RELAX',1.0_rp,'#Tracking relaxation')
      else if(words(1)=='ACCUR') then
         if(a%Listener%exists('DYNAM')) a%kfl_tacsg=1
         if(a%Listener%exists('QUASI')) a%kfl_tacsg=0
    
      endif 
     
   elseif (itask == 100) then   !Finalize reading operations  
      
      if (a%kfl_tacsg > 0) a%kfl_trasg = 1
      if (a%kfl_nolsg > 0 .and. a%kfl_tacsg == 0) then
         a%kfl_trasg = 1
         a%kfl_tacsg = 1
         write(*,*) 'WARNING: Quasi-static nonlinear SGS not implemented, Dynamic SGS set'
      end if
      if (a%kfl_stabm == 0) then
         a%kfl_trasg = 0
         a%kfl_tacsg = 0
         a%kfl_nolsg = 0
      end if
      if (a%kfl_stabm == -1) then
         a%kfl_trasg = 0
         a%kfl_tacsg = 0
         a%kfl_nolsg = 0
      end if
      if (a%kfl_adapsgs == 1) a%kfl_trasg = 1
   endif
      
end subroutine lmn_reanut

