subroutine tem_reaphy(a,itask)
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_MPIObject
   use Mod_Temperature
   implicit none
   
   integer(ip) :: itask
   class(TemperatureProblem) :: a
   integer(ip) :: imate,istat,npara,ifunp,nfunp,ndime
   !For a%Listener%listener

   integer(ip) :: idnum
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
 
   call a%Listener%getarrs(words,param,nnpar,nnwor)
      
   if (itask == 0) then
   !Initializations
   a%kfl_advec = 0                                    ! Convection is off
   a%kfl_sourc = 0                                    ! a%sources are off
   a%kfl_cotur = 0                                    ! a%turbulence coupling
   a%kfl_joule = 0                                    ! Joule effect is off
   a%kfl_CouplingThreeField = 0                       ! Three Field Coupling
   
   a%ipsou     = 0
   
   a%react     = 0.0_rp                               ! a%reaction term
   a%prtur     = 0.0_rp                               ! a%turbulent Prandtl number Prt = 0
   a%rpsou     = 0.0_rp                               ! a%source term parameters
   a%sourc     = 0.0_rp                               ! a%source term
   a%turbu     = 0.0_rp 
   
   a%visco = 1.0_rp
   a%react = 0.0_rp


   !Problem data
   elseif (itask == 1) then 
         if(words(1)=='TEMPO') then                    ! Temporal evolution
         if(a%Listener%exists('ON   ')) a%kfl_timei = 1 
      else if(words(1)=='CONVE') then               ! Convective term
         if(a%Listener%exists('ON   ')) then
            a%kfl_advec = 1
            if(words(3)=='VELOC') then
               if(words(4)=='EXTER') then
                  a%kfl_advec = 1
               else
                  a%kfl_advec = a%Listener%getint('VELOC',0,'#Velocity function')
               end if
            end if
         end if
      else if(words(1).eq.'TURBU') then
         if(words(2)=='LESMO') then
            if(a%Listener%exists('SMAGO')) then
               a%kfl_cotur=-1
               !Turbu = Cs^2
               a%turbu=a%Listener%getrea('PARAM',0.0_rp,'#Coefficient c')
               
            elseif(a%Listener%exists('WALE ')) then
               a%kfl_cotur=-2
               !Turbu = Cw
               a%turbu=a%Listener%getrea('PARAM',0.0_rp,'#Coefficient c')
            end if   
          end if  
       else if(words(1).eq.'THREE') then      
            if (words(2)=='WLF' .or. words(2)=='ON') then
               a%kfl_CouplingThreeField=1  
            end if
      else if(words(1)=='SOURC') then               ! a%source term
         if(a%Listener%exists('CONST')) then
            a%kfl_sourc = 1
         else if(a%Listener%exists('VARIA')) then
            a%kfl_sourc = 2
         else if(a%Listener%exists('USERD')) then
               a%kfl_sourc = 3
         end if
         if(a%kfl_sourc/=0) then
            a%sourc(2) = a%Listener%getrea('VALUE',0.0_rp,'#SOURCE term')
            a%sourc(1) = a%sourc(2)
            if(a%Listener%exists('JOULE')) a%kfl_joule = 1
         end if
         if(a%kfl_sourc==3) then
            call a%Listener%listen('tem_reaphy')
            do while(words(1)/='ENDSO')
               if(words(1)=='SPATI') then
                  a%rpsou(1:10)=param(1:10)
               else if(words(1)=='TIMEF') then
                  if(words(2)=='PARAB') then
                     a%ipsou(1)=1
                     a%ipsou(2)=6
                  else if(words(2)=='PERIO') then
                     a%ipsou(1)=2
                     a%ipsou(2)=6
                  else if(words(2)=='DISCO') then
                     a%ipsou(1)=5
                     a%ipsou(2)=6
                  else if(words(2)=='DISCR') then
                     a%ipsou(1)=3
                     a%ipsou(2)=2*a%Listener%getint('NUMBE',1,'#Number of data')
                  else
                     a%ipsou(1)=0
                  end if
                  if(a%ipsou(1)/=0) then
                     call a%Memor%alloc(a%ipsou(2),a%tfsou,'tfsou','tem_reaphy')
                     if(a%ipsou(1)==3) then
                        ifunp=0
                        nfunp=a%ipsou(2)/2
                        if(nfunp<1) call runend('TEMPER: WRONG DISCRETE FUNCTION PARAMETER')
                        call a%Listener%listen('tem_reaphy')
                        do while(words(1)/='ENDFU')
                           ifunp=ifunp+1
                           if(ifunp>nfunp) call runend('TEMPER: WRONG DISCRETE FUNCTION DATA')
                           a%tfsou((ifunp-1)*2+1)=param(1)
                           a%tfsou((ifunp-1)*2+2)=param(2)
                           call a%Listener%listen('tem_reabcs')
                        end do
                        ! Order the function field
                        call ordena(nfunp,a%tfsou)
                     else
                        a%tfsou(1:6)=param(3:8)
                     end if
                  end if
               end if
               call a%Listener%listen('tem_reaphy')
            end do
         end if
      else if(words(1)=='DUSTT') then               ! Dust Transport
         if (a%Listener%exists('ON   ')) then
            a%DT_kfl_DustTransport = 1
            call a%Listener%listen('tem_reaphy')
            do while(words(1)/='ENDDU')
               if ((words(1) == 'PARTI') .and. (words(2) == 'SIZE ')) then
                  a%DT_ParticleSize = param(2)
               elseif ((words(1) == 'MEDIU') .and. (words(2) == 'DENSI')) then
                  a%DT_MediumDensity = param(2)
               elseif ((words(1) == 'MEDIU') .and. (words(2) == 'DYNAM')) then
                  a%DT_MediumDynamicViscosity = param(3)   
               elseif ((words(1) == 'GRAVI') ) then
                  a%DT_GravityForce = param(2:4)   
               endif 
               call a%Listener%listen('tem_reaphy')
            enddo
         endif
      
      end if
   
   !Properties
   elseif(itask == 2) then  
      !Old compatibility
      if(words(1)=='DENSI') then                    ! a%density (rho)
         a%NumberOfMaterials = 1
         allocate(a%Materials(a%NumberOfMaterials))
         call a%Memor%allocObj(0,'Materials','tempe_memall',1*a%NumberOfMaterials)
         call a%Materials(1)%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
         call a%Materials(1)%SetDensity(param(1))
      elseif(words(1)=='SPECI') then               ! Specific heat (Cp)
         call a%Materials(1)%SetSpecificHeat(param(1))
      elseif(words(1)=='THERM') then               ! Thermal conductivity (k)
         call a%Materials(1)%SetThermalConductivity(param(1))
      endif      

      !New reading
      if(a%Listener%words(1) == 'NUMBE' .and. a%Listener%words(3) == 'MATER') then
         a%NumberOfMaterials = a%Listener%param(3)
         allocate(a%Materials(a%NumberOfMaterials))
         call a%Memor%allocObj(0,'Materials','tempe_memall',1*a%NumberOfMaterials)
      elseif (a%Listener%words(1) == 'MATER') then
         do while(a%Listener%words(1)/='ENDMA')
            call a%Listener%Listen('tempe_reaphy')
            if (a%Listener%words(1) == 'BEGIN') then
               idnum = -1
               do while(a%Listener%words(1)/='ENDMA')
                  call a%Listener%Listen('tempe_reaphy')
                  if(a%Listener%words(1) == 'IDNUM') then
                     idnum = a%Listener%param(1)
                     if(idnum > a%NumberOfMaterials) call runend('The ID material number can not be greater than Number of Materials')
                     call a%Materials(idnum)%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
                     call a%Materials(idnum)%ReadData(a%Listener)
                  endif   
                  
               enddo              
            endif
 
         enddo
      
      else if(words(1)=='VISCO') then               ! a%viscosity (mu)
         a%visco=param(1)
      else if(words(1)=='REACT') then               ! a%reaction term (s)
         a%react = a%Listener%getrea('REACT',0.0_rp,'#a%reaction parameter')
      else if(words(1)=='TURBU') then               ! a%turbulent Prandtl number
         a%prtur = a%Listener%getrea('TURBU',0.0_rp,'#a%turbulent Prandtl number')
      endif
   
   !Other initializations
   elseif (itask == 100) then

   
   endif

end subroutine tem_reaphy

