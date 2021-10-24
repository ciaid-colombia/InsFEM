subroutine nsi_reaphy(a,itask)
  use typre
  use Mod_NavierStokes
  implicit none
  
  integer(ip) :: itask
  class(NavierStokesProblem) :: a
  
  integer(ip) :: istat,npara,ipara,aux,idime,ndime
  real(rp)    :: dummr,gnorm
  
  !For Listener
  real(rp), pointer     :: param(:) => NULL()
  character(5), pointer :: words(:) => NULL()
  integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
  integer(ip)           :: imat  =1_ip,irest
  character(150) :: outstr
  
  call a%Listener%getarrs(words,param,nnpar,nnwor)
  outstr = adjustl(trim(a%exmod))//'_REAPHY'  
  
  !Initializations
  if (itask == 0) then
      a%kfl_advec = 0                                    ! Convection is off
      a%kfl_confi = 0                                    ! Flow is not confined
      a%kfl_visco = 1                                    ! Viscous term on
      a%fcons     = 0.0_rp                               ! Non conservative form default
      a%fvins     = 1.0_rp                               ! Divergence form default
      a%MatProp(:)%visco=1.0_rp
      a%MatProp(:)%densi=1.0_rp
      a%MatProp(:)%lawvi=0
      
      a%kfl_advco = 0                                    ! Constant convection is off
      a%advco     = 0.0_rp                               ! Constant convection velocity vector
      a%grnor     = 0.0_rp                               ! Gravity norm
      a%gravi     = 0.0_rp                               ! Gravity vector

      a%kfl_cotem = 0                                    ! No coupling with temperature
      a%boube     = 0.0_rp                               ! Beta
      a%boutr     = 0.0_rp                               ! Tr
      a%bougr     = 0.0_rp                               ! Gravity for Boussinesq
      a%kfl_cotur = 0
      a%kfl_SwitchOff = 0_ip

      a%kfl_colev = 0_ip                                  !No coupling with levelset
      a%kfl_fsurf = 0_ip                                  !Free surface flag
      a%kfl_fsurfLapla = 0_ip                             !Free surface with laplacian problem outside
      a%kfl_fsurfDirichlet = 0_ip                         !DirichletBoundary Conditions on the free surface
      a%kfl_EnrichElem = 0_ip                             !Enriched element flag used in cut elements
      a%kfl_EnrichPressure=0_ip                           !Pressure enriched flag
      a%kfl_EnrichVelocity=0_ip                           !Velocity enriched flag
      a%nmat=1                                            !Number of Materials

      a%kfl_FORAcceleration = 0
      a%kfl_FORFunty = 0
      a%FORParam = 0
      a%FORAcceleration = 0.0_rp
      a%kfl_FORAxesRotation = 0_ip                        !Non-inertial rotating frame of reference
      a%FORAxesAngularVeloc = 0.0_rp

      a%kfl_TurbulentInletBoundaryConditions = 0_ip
      a%kfl_PorousMedia = 0_ip                            ! Porous Media 
      a%kfl_Plasma = 0_ip
      
      a%kfl_RVE = 0_ip
      
   !Problem data   
   elseif (itask == 1) then
      if(words(1)=='CONVE') then               ! Convective term
         if(words(2)/='OFF  ') a%kfl_advec = 1
         if(words(2)=='CONST') then
             a%kfl_advco = 1
             a%advco(1) = a%Listener%getrea('VX   ',0.0_rp,'#x-component of a')
             a%advco(2) = a%Listener%getrea('VY   ',0.0_rp,'#y-component of a')
             a%advco(3) = a%Listener%getrea('VZ   ',0.0_rp,'#z-component of a')
         end if
         if(words(2)=='NONCO') then
            a%fcons=0.0_rp
         else if(words(2)=='CONSE') then
            a%fcons=1.0_rp                 
         else if(words(2)=='SKEWS') then                 
            a%fcons=0.5_rp
         end if                 
         if(words(3)=='VELOC') then
            if(words(4)=='NAVIE') then
               a%kfl_advec = 1
            else
               a%kfl_advec = a%Listener%getint('VELOC',0,'#Velocity function')
            end if
         end if
      else if(words(1)=='VISCO') then               ! Viscous term
         if(words(2)=='OFF  ') a%kfl_visco=0
         if(words(2)=='DIVER') a%fvins=1.0_rp
         if(words(2)=='ON   ') a%fvins=1.0_rp     ! Divergence form default
         if(words(2)=='LAPLA') a%fvins=0.0_rp     ! Laplacian form
         if(words(2)=='COMPL') a%fvins=2.0_rp     ! Complete form
      else if(words(1).eq.'GRAVI') then
         a%grnor    = a%Listener%getrea('NORM ',0.0_rp,'#Gravity norm')
         a%gravi(1) = a%Listener%getrea('GX   ',0.0_rp,'#x-component of g')
         a%gravi(2) = a%Listener%getrea('GY   ',0.0_rp,'#y-component of g')
         a%gravi(3) = a%Listener%getrea('GZ   ',0.0_rp,'#z-component of g')
         call vecuni(3,a%gravi,dummr)
      else if((words(1)=='FRAME') .and. a%Listener%exists('ON   ')) then
         a%kfl_FORAcceleration=1
         a%kfl_FORFunty = a%Listener%getint('FUNTY ',0_ip,'#Function type')
         a%FORParam(1:6)=a%Listener%param(4:9)
         a%FORAcceleration = a%Listener%param(10:12)
      else if((words(1)=='CORIO') .and. a%Listener%exists('ON   ')) then
         a%kfl_FORAxesRotation=1
         a%FORAxesAngularVeloc(1:3)=a%Listener%param(2:4)
     else if((words(1)=='TEMPE'.and.words(2)=='BOUSS').or.&
            &  (words(1)=='BOUSS'.and.a%Listener%exists('ON   '))) then
         a%kfl_cotem=1
         a%boube = a%Listener%getrea('BETA ',0.0_rp,'#Beta coefficient')
         a%boutr = a%Listener%getrea('TA   ',0.0_rp,'#Reference temperature')
         a%bougr(1) = a%Listener%getrea('GX   ',0.0_rp,'#x-component of g')
         a%bougr(2) = a%Listener%getrea('GY   ',0.0_rp,'#y-component of g')
         a%bougr(3) = a%Listener%getrea('GZ   ',0.0_rp,'#z-component of g')
      else if(words(1).eq.'TURBU') then
         if(words(2)=='LESMO') then
            if(a%Listener%exists('SMAGO')) then                !Smagorinsky LES model
               a%kfl_cotur=-1
               !turbu = Cs^2
               a%turbu(1)=a%Listener%getrea('PARAM',0.0_rp,'#Coefficient c')  
            
            elseif (a%Listener%exists('WALE ')) then
               a%kfl_cotur=-2
               !turbu = Cw
               a%turbu(1)=a%Listener%getrea('PARAM',0.0_rp,'#Coefficient c')  
            
            endif
         endif
      else if(words(1).eq.'SWITC') then
         if(words(2)=='ON') then
               a%kfl_SwitchOff = 1_ip
         endif
      else if(words(1) == 'LEVEL'.and. words(2)=='ON') then
         a%kfl_colev=1
      else if(words(1) == 'FREES')then
         if(words(2) == 'ON')then
            a%kfl_fsurf=1
            if(words(3) == 'LAPLA') a%kfl_fsurfLapla=1
            if(words(3) == 'STOKE') a%kfl_fsurfLapla=2
            if(words(3) == 'DINAV') a%kfl_fsurfLapla=3
            if (words(4) == 'DIRIC') then
               a%kfl_fsurfDirichlet = 1
            elseif (words(4) == 'DIRI0') then
               a%kfl_fsurfDirichlet = 2
            endif
            if (words(5) == 'NITSC') then
               a%kfl_fsurfDirichletMethod = 1
            elseif (words(5) == 'STRON') then
               a%kfl_fsurfDirichletMethod = 2   
            elseif (words(5) == 'LLM  ') then
               a%kfl_fsurfDirichletMethod = 3   
            endif
               
         elseif(words(2) == 'Off')then
            a%kfl_fsurf=0
         end if      
      else if(words(1) == 'NUMBE')then
         a%nmat=param(1)
         if (a%nmat > 2) call runend('NavierStokes not ready for more than 2 materials')
      else if(words(1) == 'ENRIC')then
         if(words(2) == 'ON')then
            a%kfl_EnrichElem=1
            if(words(3) == 'PRESS') a%kfl_EnrichPressure=1
            if(words(3) == 'VELOC') a%kfl_EnrichVelocity=1
            if(words(3) == 'BOTH')then
               a%kfl_EnrichPressure=1
               a%kfl_EnrichVelocity=1
            end if
         elseif(words(2) == 'Off')then
            a%kfl_EnrichElem=0
         end if    
      !Elastic boundary Dirichlet boundary conditions
      elseif(words(1) == 'ELAST') then
         a%EB_Stiffness = a%Listener%getrea('STIFF',100.0_rp,'Stiffness')
         a%EB_Damping = a%Listener%getrea('DAMPI',0.0_rp,'Damping')
         a%EB_Density = a%Listener%getrea('DENSI',1.0_rp,'Damping')  
          a%EB_ShearStiffness = a%Listener%getrea('SHEAR',100.0_rp,'ShearStiffness')
         a%EB_NDelaySteps = a%Listener%getint('DELAY',0,'DelaySteps')
      
      else if(words(1) == 'RESTR')then
         if(words(2)=='ON') then
            a%kfl_restrictions = 1_ip
            call a%Listener%listen(outstr)
            if(words(1)=='NREST') then
               a%nrest = param(1)
               do irest=1,a%nrest
                  call a%Listener%listen(outstr)
                  a%keyrest(irest)=words(2)
               enddo
            endif
         endif 

      else if(words(1)=='INLET') then   
         if (a%Listener%exists('ON   ')) then
            a%kfl_TurbulentInletBoundaryConditions = 1
            call a%TIBC%ReadData(a%Listener)
         endif 
      else if(words(1)=='EXITB') then   
         if (a%Listener%exists('ON   ')) then
            a%kfl_ExitBodyForces = 1
            call a%EBF%ReadData(a%Listener)
         endif    
      
      else if(words(1)=='TURBO') then   
         if (a%Listener%exists('ON   ')) then
            a%kfl_TurbulentBodyForces = 1
            call a%TBF%ReadData(a%Listener)
         endif 
         
      elseif(words(1)=='EXTRA' .and. words(2) == 'VISCO') then
         !Extra Initial Viscosity
         if (a%Listener%exists('ON   ')) then
            a%kfl_ExtraInitialViscosity = 1
            a%EIV_visco = param(3)
            a%EIV_nsteps = param(4)
         endif
         
      elseif(words(1)=='PORUS') then             ! Porous Media
         if (a%Listener%exists('ON   ')) then
            a%kfl_PorousMedia = 1
            call a%PME%ReadData(a%Listener)
         endif
         
      elseif(words(1)=='PLASM') then             ! Plasma
         if (a%Listener%exists('ON   ')) then
            a%kfl_Plasma = 1
            call a%PAC%ReadData(a%Listener)
         endif
         
      elseif(words(1)=='RVE  ') then             ! Plasma
         if (a%Listener%exists('ON   ')) then
            a%kfl_RVE = 1
         endif

      elseif(words(1)=='CORIO')then   
         if (a%Listener%exists('ON   ')) then
            a%kfl_CoriolisForce = 1
            a%CoriolisW(1:3) = a%Listener%param(2:4)
         endif
      endif
   
   !Properties
     
   elseif(itask == 2) then
      if(words(1)=='MATER') then 
         imat = param(1)
         if (imat > a%nmat) call runend('wrong material number')
      elseif(words(1)=='DENSI') then                 ! Density (rho)
         a%MatProp(imat)%densi=param(1)         
      elseif(words(1)=='VISCO') then             ! Viscosity (mu)
         a%MatProp(imat)%visco=param(1)
      elseif(words(1)=='LAWVI')then
         if(words(2)=='CONST')then
            a%MatProp(imat)%lawvi=0
         elseif(words(2)=='POWER')then           ! Power Law model
            a%MatProp(imat)%lawvi=1           
            a%MatProp(imat)%LawViParam(1:2) = param(2:3)            
         elseif(words(2)=='CARRE')then           ! Careau-Yasuda model 
            a%MatProp(imat)%lawvi=2  
            a%MatProp(imat)%LawViParam(1:5) = param(2:6)            
         end if
      end if  
      
      
         
      !Final operations      
      elseif(itask==100) then 
      
         !If Free surface, do NOT allow different materials
         !Important for stabilization, which requires the same material in all elements
         if (a%kfl_fsurf == 1) then
            a%MatProp(2:size(a%MatProp)) = a%MatProp(1) 
         endif
      
   endif
   
end subroutine nsi_reaphy
