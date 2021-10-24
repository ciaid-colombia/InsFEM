subroutine lmn_reaphy(a,itask)
  use typre
  use Mod_LowMach
  implicit none
  integer(ip) :: itask
  class(LowMachProblem) :: a
  real(rp)    :: dummr
  !For Listener
  real(rp), pointer     :: param(:) => NULL()
  character(5), pointer :: words(:) => NULL()
  integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
  
  call a%Listener%getarrs(words,param,nnpar,nnwor)
  
  !Initializations
  if (itask == 0) then
      a%kfl_eqnst = 0                                   ! Equation of state
      a%kfl_advec = 0                                    ! Convection is off
      a%kfl_confi = 0                                    ! Flow is not confined
      a%kfl_visco = 1                                    ! Viscous term on
!      a%kfl_cotur = 0
!      a%fcons     = 0.0_rp                               ! Non conservative form default
      a%fvins     = 1.0_rp                               ! Divergence form default
      a%kfl_sourc = 0                                    ! heat sources are off
      a%kfl_pther = .true.                                    ! heat sources are off
!      a%kfl_joule = 0                                    ! Joule effect is off

      a%itpre     = 0.0_rp                               ! Initial thermodynamic pressure
      a%visco     = 1.0_rp                               ! Dynamic viscosity (mu)
      a%cphea     = 1.0_rp                               ! Specific heat constant pressure (cp)
      a%tcond     = 1.0_rp                               ! Thermal conductivity (K)
      a%texpc     = 1.0_rp                               ! Thermal expansion (alpha)
      a%grnor     = 0.0_rp                               ! Gravity norm
      a%gravi     = 0.0_rp                               ! Gravity vector
      a%sgasc     = 0.0_rp                               ! Specific gas constant (R) 
      a%react     = 0.0_rp                               ! a%reaction term
!      a%prtur     = 0.0_rp                               ! a%turbulent Prandtl number Prt = 0
      a%sourc     = 0.0_rp                               ! a%source term
!      a%turbu     = 0.0_rp 
!      a%lawsp     = 0
!      a%lawtc     = 0
!      a%lawvi     = 0
      a%epspe     = 0.0_rp                               ! Penalization in pressure for closed flows 

   !Problem data   
   elseif (itask == 1) then

      if(words(1)=='CONVE') then               ! Convective term
         if(words(2)/='OFF  ') a%kfl_advec = 1
!         if(words(2)=='NONCO') then
!            a%fcons=0.0_rp
!         else if(words(2)=='CONSE') then
!            a%fcons=1.0_rp                 
!         else if(words(2)=='SKEWS') then                 
!            a%fcons=0.5_rp
!         end if                 
!         if(words(3)=='VELOC') then
!            if(words(4)=='NAVIE') then
!               a%kfl_advec = 1
!            else
!               a%kfl_advec = a%Listener%getint('VELOC',0,'#Velocity function')
!            end if
!         end if
!      else if(words(1)=='VISCO') then               ! Viscous term
!         if(words(2)=='OFF  ') a%kfl_visco=0
!         if(words(2)=='DIVER') a%fvins=1.0_rp
!         if(words(2)=='ON   ') a%fvins=1.0_rp     ! Divergence form default
!         if(words(2)=='LAPLA') a%fvins=0.0_rp     ! Laplacian form
!         if(words(2)=='COMPL') a%fvins=2.0_rp     ! Complete form
      else if(words(1).eq.'GRAVI') then
         a%grnor    = a%Listener%getrea('NORM ',0.0_rp,'#Gravity norm')
         a%gravi(1) = a%Listener%getrea('GX   ',0.0_rp,'#x-component of g')
         a%gravi(2) = a%Listener%getrea('GY   ',0.0_rp,'#y-component of g')
         a%gravi(3) = a%Listener%getrea('GZ   ',0.0_rp,'#z-component of g')
         call vecuni(3,a%gravi,dummr)
!      else if(words(1).eq.'TURBU') then
!         if(words(2)=='LESMO') then
!            if(a%Listener%exists('SMAGO')) then                ! Smagorinsky LES model
!               a%kfl_cotur=-1
!               a%turbu(1)=a%Listener%getrea('PARAM',0.0_rp,'#Coefficient c')  
!            endif
!         endif
      else if(words(1)=='SOURC') then               ! a%source term
         if(a%Listener%exists('CONST')) then
            a%kfl_sourc = 1
         else if(a%Listener%exists('VARIA')) then
            a%kfl_sourc = 2
         end if
         if(a%kfl_sourc/=0) then
            a%sourc = a%Listener%getrea('VALUE',0.0_rp,'#SOURCE term')
         end if
      end if
   
   !Properties
     
  elseif(itask == 2) then

      if(words(1)=='EQNST') then                   ! Equation of state
         if(a%Listener%exists('CONST')) then
            a%kfl_eqnst=0
         else if(a%Listener%exists('IDEAL')) then
            a%kfl_eqnst=1
         end if
      elseif(words(1)=='DENSI') then               ! Initial density
         a%densi=param(1) 
      elseif(words(1)=='VISCO') then               ! Viscosity (mu)
         a%visco=param(1) 
      elseif(words(1)=='SPECI') then               ! Specific heat cons pressure(Cp)
         a%cphea=param(1)      
      elseif(words(1)=='TCOND') then               ! Thermal conductivity (k)
         a%tcond=param(1)      
      elseif(words(1)=='TEXPA') then               ! Thermal expansion (alpha)
         a%texpc=param(1)      
      elseif(words(1)=='GASCO') then               ! Specific gas constant
         a%sgasc=param(1)
      elseif(words(1)=='TPRES') then               ! Thermodynamic pressure
         a%itpre=param(1)
      elseif(words(1)=='PCONS') then               ! Thermodynamic pressure
         a%kfl_pther = .false.
!       else if(words(1)=='LAWSP') then               ! Law Specific heat (Cp)
!         a%lawsp=0
!      else if(words(1)=='LAWTH') then               ! Law Thermal conductivity (k)
!         a%lawtc=0
!         if(words(2)=='VARIA') a%lawtc=1
!      else if(words(1)=='LAWVI') then               ! Law Thermal a%viscosity
!         a%lawvi=0
      else if(words(1)=='REACT') then               ! a%reaction term (s)
         a%react = a%Listener%getrea('REACT',0.0_rp,'#a%reaction parameter')
!      else if(words(1)=='TURBU') then               ! a%turbulent Prandtl number
!         a%prtur = a%Listener%getrea('TURBU',0.0_rp,'#a%turbulent Prandtl number')
      end if            
       
   !Final operations      
   elseif(itask==100) then      
   
   endif
   

end subroutine lmn_reaphy
