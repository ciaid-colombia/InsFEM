subroutine opt_reaphy(a,itask)
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_MPIObject
   use Mod_Optics
   implicit none
   
   integer(ip) :: itask
   class(OpticsProblem) :: a
	
	integer(ip) :: imate,istat,npara,ifunp,nfunp
	!For a%Listener%listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
 
   call a%Listener%getarrs(words,param,nnpar,nnwor)
      
   if (itask == 0) then
   !Initializations
   a%kfl_velocity = 1                                    !Velocity is external
   a%kfl_pressure = 1                                    !Pressure is external
   a%kfl_temperature = 1                                 !Temperature is external
   
   
   a%densi = 1.0_rp
   a%sphea = 1.0_rp
   a%tcond = 1.0_rp
   a%visco = 1.0_rp
   a%prtur = 1.0_rp
   a%units = 0
   
   !Problem data
   elseif (itask == 1) then 
      if(words(1)=='FLOWF') then               ! Convective term
         if(a%Listener%exists('EXTER')) then
            a%kfl_velocity = 1
            a%kfl_pressure = 1
         else
            call runend('opt_reaphy: Velocity field only implemented for external')
         endif 
      else if(words(1)=='TEMPE') then               ! Temperature term
         if(a%Listener%exists('EXTER')) then
            a%kfl_temperature = 1
         else
            call runend('opt_reaphy: Temperature field only implemented for external')
         endif    
      endif
   
   !Properties
   elseif(itask == 2) then      
      if(words(1)=='DENSI') then                    ! a%density (rho)
         a%densi = param(1)
      else if(words(1)=='SPECI') then               ! Specific heat (Cp)
         a%sphea=param(1)
      else if(words(1)=='THERM') then               ! Thermal conductivity (k)
         a%tcond=param(1)
      else if(words(1)=='VISCO') then               ! a%viscosity (mu)
         a%visco=param(1)
      else if(words(1)=='TURBU') then               ! a%turbulent Prandtl number
         a%prtur = a%Listener%getrea('TURBU',0.0_rp,'#a%turbulent Prandtl number')
      else if(words(1)=='UNITS') then               ! a%turbulent Prandtl number
         if (a%Listener%exists('METER')) a%units = 0
      endif
   
   !Other initializations
   elseif (itask == 100) then
   
   
   endif

end subroutine opt_reaphy

