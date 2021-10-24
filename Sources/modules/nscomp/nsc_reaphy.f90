subroutine nsc_reaphy(a,itask)
  use typre
  use Mod_Listen
  use Mod_NSCompressible
  implicit none
  
  integer(ip) :: itask
  class(NSCompressibleProblem) :: a
  
  integer(ip) :: ipara,idamp,idime,ndime
  real(rp)    :: dummr
  
  !For Listener
  real(rp), pointer     :: param(:) => NULL()
  character(5), pointer :: words(:) => NULL()
  integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
  
  character(150) :: outstr
  
  call a%Listener%getarrs(words,param,nnpar,nnwor)
  outstr = adjustl(trim(a%exmod))//'_REAPHY'
  
  call a%Mesh%GetNdime(ndime)
  
  !Initializations
  if (itask == 0) then
      a%kfl_advec = 1                                    ! Convection is on
      a%kfl_visco = 1                                    ! Viscous term on
      a%kfl_react = 0                                    ! No reaction term
      a%kfl_sourc = 0                                    ! No sources
      a%kfl_confi = 0                                    ! Flow is not confined
      a%cphea     = 1.0_rp                                ! [J/kg.K]=[m2/(s2.K)]			
      a%cvhea     = 1.0_rp                                ! [J/kg.K]=[m2/(s2.K)]
      a%visco     = 1.0_rp                                ! [kg/m.s]
      a%tcond     = 1.0_rp                                ! [W/m.K]=[Kg.m/(s3.K)]  
      a%lawde     = 0                                    ! Equation of state
      a%lawdep    = 0.0_rp                               ! State law eqn Param
      a%relpre    = 0.0_rp                               ! Relative pressure
      a%reltem    = 0.0_rp                               ! Relative temperature
      a%grnor     = 0.0_rp                               ! Gravity norm
      a%gravi     = 0.0_rp                               ! Gravity vector
      a%srce      = 0.0_rp                               ! Heat sources

      a%ndamp = 0                               ! No damping zones
      a%dampco = 0.0_rp                         ! Coeficients of damping
      a%dampxo = 0.0_rp                         ! Coordinates of initial vertex rectangular region
      a%dampxf = 0.0_rp                         ! Coordinates of final vextex rectangular region
      a%dampff = 0.0_rp                         ! Unperturbed values of fields

      a%rdamp   = 0                             ! No Radial damping
      a%rdampco = 0.0_rp                        ! Coeficients of radial damping
      a%rdampxo = 0.0_rp                        ! Coordinates of origin
      a%rdampro = 0.0_rp                        ! Initial radius of damping region
      a%rdampex = 0                             ! Exponential of damping
      a%rdampff = 0.0_rp                        ! Unperturbed values of fields

      ! Restart from incompressible flow
      a%speed = 0.0_rp
      a%molmass = 0.0_rp
      a%densi = 0.0_rp
      

   call a%SpecificNSCompReaphy(0)

   !Problem data   
   elseif (itask == 1) then

      if(words(1).eq.'GRAVI') then
         a%grnor    = a%Listener%getrea('NORM ',0.0_rp,'#Gravity norm')
         a%gravi(1) = a%Listener%getrea('GX   ',0.0_rp,'#x-component of g')
         a%gravi(2) = a%Listener%getrea('GY   ',0.0_rp,'#y-component of g')
         a%gravi(3) = a%Listener%getrea('GZ   ',0.0_rp,'#z-component of g')
         call vecuni(3,a%gravi,dummr)
      else if(words(1)=='CONVE') then               ! Convective term
         if(words(2)=='OFF  ') a%kfl_advec = 0
      else if(words(1)=='VISCO') then               ! Viscous term
         if(words(2)=='OFF  ') a%kfl_visco=0
      else if(words(1).eq.'SRCES') then
         a%srce = a%Listener%getrea('SRCE ',0.0_rp,'#Heat sources')
      ! Damping zones.
      else if(words(1)=='DAMPI') then
         ipara=1
         a%ndamp=int(param(ipara))
         do idamp = 1,a%ndamp
            do idime = 1,ndime
               a%dampco(idime,idamp) = (param(ipara+idime))         
               a%dampxo(idime,idamp) = (param(ipara+idime+3))         
               a%dampxf(idime,idamp) = (param(ipara+idime+6))         
               a%dampff(idime+1,idamp) = (param(ipara+idime+10))         
            end do
            ipara=ipara+10
            a%dampff(1,idamp) = (param(ipara))         
            ipara=ipara+4
            a%dampff(5,idamp) = (param(ipara))         
         end do
      ! Radial Damping
      else if(words(1)=='RADIA') then
         a%rdamp= 1 
         a%rdampco = (param(1))         
         do idime = 1,ndime
         a%rdampxo(idime) = (param(idime+1))         
         a%rdampff(idime+1) = (param(idime+7))         
         end do
         a%rdampro = (param(5))         
         a%rdampex = int(param(6))         
         a%rdampff(1) = (param(7))         
         a%rdampff(5) = (param(11))         
      else if(words(1)=='INCOM') then               ! Viscous term
         if(words(2)=='ON   ') then
            a%kfl_nsirstar=1
            call a%Listener%listen(outstr)
            a%speed = param(1)
            call a%Listener%listen(outstr)
            a%molmass = param(1)
            call a%Listener%listen(outstr)
            a%densi = param(1)   
         endif
      end if
   
      call a%SpecificNSCompReaphy(1)

   !Properties
     
  elseif(itask == 2) then
  
      if(words(1)=='EQNST') then               ! Equation of state
         if(words(2)=='IDEAL') a%lawde = 1 ! Ideal gas eqn. of state.
      end if                 
      if(words(1)=='CPHEA') then                   ! Specific heat cons pressure(Cp)
         a%cphea=param(1)      
      elseif(words(1)=='CVHEA') then               ! Specific heat cons volume (Cv)
         a%cvhea=param(1)      
      elseif(words(1)=='VISCO') then               ! Viscosity (mu)
         a%visco=param(1)      
      elseif(words(1)=='THERM') then               ! Thermal conductivity (k)
         a%tcond=param(1)      
      elseif(words(1)=='RELAT') then
         if(words(2) == 'ON ') then 
            a%relpre = param(2)
            a%reltem = param(3)
         end if
     end if    
      
     call a%SpecificNSCompReaphy(2)
   
   !Final operations      
   elseif(itask==100) then 
      
     call a%SpecificNSCompReaphy(100)
   
   endif

   

end subroutine nsc_reaphy
