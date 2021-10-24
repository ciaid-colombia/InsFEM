subroutine nsc_reanut(a,itask)
   use typre
   use Mod_Listen
   use Mod_NSCompressible
   implicit none
   
   class (NSCompressibleProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: istab

   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   

   call a%Listener%getarrs(words,param,nnpar,nnwor)
  
   if (itask == 0) then         !  Initializations (defaults)

      a%staco(1)  = 4.0_rp                               ! Viscous part of tau
      a%staco(2)  = 2.0_rp                               ! Convective part of tau
      a%staco(3)  = 1.0_rp                               ! Reaction part of tau
      a%staco(4)  = 1.0_rp                               ! Coefficient of tau continuity eq
      a%kfl_repro = 0                                    ! Res. projection not used
      a%kfl_shock = 0                                    ! Shock capturing off
      a%kfl_sctyp = 0                                    ! Shock capturing isotropic type
      a%shock     = 0.0_rp                               ! SC parameter
      a%kfl_trasg = 0                                    ! Not tracking subscales
      a%kfl_tacsg = 0                                    ! Quasistatic subscales
      a%kfl_nolsg = 0                                    ! Linear subscales
      a%kfl_stabm = 0                                    ! Default is SUPG
      a%kfl_jacgr = 0                                    ! Default is A·grad V
      a%epspe = 0.001_rp                                 !Penalization in confined flow 
      a%ErrorEstimatorTypeOfSubscales = 0                !0: Scaled L_2 norm, 1: Entropy function. Default is L_2 norm
      
   call a%SpecificNSCompReanut(0)
   
   elseif (itask == 1) then     !Inside the Numerical Treatment Block   
      if(words(1)=='STABI') then
        do istab = 1,4
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
     elseif(words(1)=='JACOB') then 
        if (words(2) == 'ON   ') a%kfl_jacgr = 1

     else if(words(1)=='SHOCK') then
        if(a%Listener%exists('TOTAL')) then
           a%kfl_shock = 1
           if(a%Listener%exists('ISOTR')) then
              a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           elseif(a%Listener%exists('ANISO')) then
              a%kfl_sctyp = 1 
              a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if
        else if(a%Listener%exists('GRADI')) then
           a%kfl_shock = 2
           if(a%Listener%exists('ISOTR')) then
              a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           elseif(a%Listener%exists('ANISO')) then
              a%kfl_sctyp = 1 
              a%shock     = a%Listener%getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if
        end if

     else if(words(1)=='TRACK') then
        if (words(2) == 'ON   ') then
            a%kfl_trasg = 1
        else
           a%kfl_trasg=0
        end if
     else if(words(1)=='NONLI') then
        if(a%Listener%exists('LINEA')) then
           a%kfl_nolsg=0
        else
           a%kfl_nolsg=1
        end if
     else if(words(1)=='ACCUR') then
        if(a%Listener%exists('DYNAM')) a%kfl_tacsg=1
        if(a%Listener%exists('QUASI')) a%kfl_tacsg=0
    
     else if(words(1)=='EPRES') then
         a%epspe = a%Listener%getrea('EPRES',0.0001_rp,'Pressure penalization')
     elseif (words(1) == 'ERRSU') then
         if (a%Listener%exists('SCALE')) then
            a%ErrorEstimatorTypeOfSubscales = 1
         elseif (a%Listener%exists('ENTROP')) then
            a%ErrorEstimatorTypeOfSubscales = 2
         endif
     endif 
     
      call a%SpecificNSCompReanut(1)

   elseif (itask == 100) then   !Finalize reading operations  
      
      if (a%kfl_tacsg > 0) a%kfl_trasg = 1
      if (a%kfl_nolsg > 0) a%kfl_trasg = 1
      if (a%kfl_repro > 0) a%kfl_trasg = 1
      if (a%RefinerErrorEstimator == 'SUBSC') a%kfl_trasg = 1
      
      call a%SpecificNSCompReanut(100)     !Final operations

   endif
      

end subroutine nsc_reanut
