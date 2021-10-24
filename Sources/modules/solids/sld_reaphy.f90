subroutine sld_reaphy(a,itask)
  use typre
  use Mod_Solids
  implicit none
  class(SolidsProblem) :: a
  integer(ip) :: itask
  real(rp)    :: dummr,gnorm
  
  !For Listener
  real(rp), pointer     :: param(:) => NULL()
  character(5), pointer :: words(:) => NULL()
  integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
  
  call a%Listener%getarrs(words,param,nnpar,nnwor)

  !Initializations
  if (itask == 0) then

      a%grnor     = 0.0_rp                                ! Gravity norm
      a%gravi     = 0.0_rp                                ! Gravity vector
      a%densi     = 0.0_rp                                ! Gravity vector
      !Linear elasticity
      a%young     = 0.0_rp                                ! youngs modulus
      a%poisson     = 0.0_rp                              ! poisson modulus
      !Non-Linear elasticity
      a%lambda = 0.0_rp                                   ! lambda (lame parm)
      a%mu = 0.0_rp                                       ! mu (lame param) 
      a%kfl_traction = 0                                  ! tractions are off
      !a%rpsou     = 0.0_rp                               ! a%source term parameters
      a%traction     = 0.0_rp                              ! a%source term
      a%sld_model= 'LINEA'

  elseif (itask == 1) then

      if(words(1).eq.'GRAVI') then
          a%grnor    = a%Listener%getrea('NORM ',0.0_rp,'#Gravity norm')
          a%gravi(1) = a%Listener%getrea('GX   ',0.0_rp,'#x-component of g')
          a%gravi(2) = a%Listener%getrea('GY   ',0.0_rp,'#y-component of g')
          a%gravi(3) = a%Listener%getrea('GZ   ',0.0_rp,'#z-component of g')
          call vecuni(3,a%gravi,dummr)
      else if(words(1)=='BODYF') then               ! a%tractione term
            a%kfl_traction = 1
         if(a%kfl_traction/=0) then
            a%traction(1) = a%Listener%getrea('FX   ',0.0_rp,'#tractionX term')
            a%traction(2) = a%Listener%getrea('FY   ',0.0_rp,'#tractionY term')
            a%traction(3) = a%Listener%getrea('FZ   ',0.0_rp,'#tractionZ term')
         end if
      end if

  elseif (itask == 2) then

      if(words(1).eq.'SLDMO') then
          a%sld_model = words(2)
      endif

      if(words(1).eq.'DENSI') then
          a%densi    = a%Listener%getrea('DENSI',0.0_rp,'#Density')
      endif

      if(words(1).eq.'YOUNG') then
          a%young    = a%Listener%getrea('YOUNG',0.0_rp,'#Young')
      endif
      if(words(1).eq.'POISS') then
          a%poisson  = a%Listener%getrea('POISS',0.0_rp,'#Poisson')
      endif

  !Final operations      
   elseif(itask==100) then      

   endif

   !write(*,*),'neohookean solid: lambda= ',a%lambda,' , mu= ',a%mu

end subroutine sld_reaphy
