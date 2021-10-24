subroutine sup_reaphy(a,itask)
   use typre
   use Mod_Listen
   use Mod_ThreeField
   implicit none
   
   integer(ip) :: itask
   class(ThreeFieldNSProblem) :: a
   
   integer(ip) :: istat,npara,ipara,aux
   real(rp)    :: dummr,gnorm
   
   !For Listener
   real(rp), pointer     :: param(:)
   character(5), pointer :: words(:)
   integer(ip), pointer  :: nnpar,nnwor
   integer(ip)           :: aux_BC_number
   integer(ip)           :: imat
   character(150)        :: outstr   
   
   interface
      subroutine nsi_reaphy(a,itask)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
         
      end subroutine
   end interface
   
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   outstr = adjustl(trim(a%exmod))//'_REAPHY'    
   
   !Call the default one
   if(itask==0 .or. itask==1) call nsi_reaphy(a,itask)
   
   
  
   !Initializations
   if (itask == 0) then
      a%kfl_bc_number=0
      
      !Constitutive models
      a%PTT_model=1
      a%Giesekus_model=0
      
      !LCR reformulation
      a%LogFormulation=0 
      
      !Temperature coupling
      a%kfl_cotem=0
      a%kfl_cotem_Boussinesq=0
      a%kfl_cotem_WLF=0 
      a%kfl_cotem_Arrhenius=0 
      a%ReferenceTemp=0.0_rp
      a%c1_WLF=0.0_rp
      a%c2_WLF=0.0_rp
      a%alpha_Arrhenius=0.0_rp
      
   !Problem data   
   elseif (itask == 1) then
   
      if(words(1)=='BOUND')then
         a%kfl_bc_number= param(1)
      elseif(words(1)=='TEMPE') then 
         if (words(2)=='BOUSS') then !Boussinesq model
            a%kfl_cotem = 1 
            a%kfl_cotem_Boussinesq = 1
            !The parameters are setted in nsi_reaphy
         else if (words(2)=='WLF') then !WLF Model
            a%kfl_cotem = 1 
            a%kfl_cotem_WLF = 1
            a%ReferenceTemp= a%Listener%getrea('T0   ',0.0_rp,'#Reference temperature')
            a%c1_WLF= a%Listener%getrea('C1   ',0.0_rp,'#Parameter c1')
            a%c2_WLF= a%Listener%getrea('C2   ',0.0_rp,'#Parameter c2')
         else if (words(2)=='ARRHE') then
            a%kfl_cotem = 1 
            a%kfl_cotem_Arrhenius = 1
            a%ReferenceTemp= a%Listener%getrea('T0   ',0.0_rp,'#Reference temperature')
            a%alpha_Arrhenius= a%Listener%getrea('ALPHA',0.0_rp,'#Parameter alpha')
         end if
      end if   
      
      if(words(1)=='LOGAR') then
         if(words(2)=='ON') then
            a%LogFormulation=1
            a%nu0_LCR=param(2) !Factor value.
            if (a%nu0_LCR==0.0_rp) a%nu0_LCR=0.01_rp    
         end if   
      end if
      

  !Properties
  elseif(itask == 2) then
      if(a%nmat==1)then         
         if(words(1)=='DENSI') then                 ! Density (rho)
            a%MatProp(1)%densi=param(1)         
         elseif(words(1)=='VISCO') then             ! Viscosity (mu)
            a%MatProp(1)%visco=param(1)
         elseif(words(1)=='LAWVI')then
            if(words(2)=='CONST')then
               a%MatProp(1)%lawvi=0
            elseif(words(2)=='POWER')then           ! Power Law model
               a%MatProp(1)%lawvi=1           
               a%MatProp(1)%LawViParam(1:2) = param(2:3)            
            elseif(words(2)=='CARRE')then           ! Careau-Yasuda model 
               a%MatProp(1)%lawvi=2  
               a%MatProp(1)%LawViParam(1:5) = param(2:6) 
            elseif(words(2)=='OLDRO')then               ! Oldroyd-B model  
               a%MatProp(1)%lawvi=-1  
               a%MatProp(1)%LawViParam(1:3) = param(2:4) 
               a%MatProp(1)%LawViParam(4) = 0.0_rp  
               a%MatProp(1)%LawViParam(5)= a%nu0_LCR
            elseif(words(2)=='GIESE')then               ! Giesekus  
               a%MatProp(1)%lawvi=-2  
               a%MatProp(1)%LawViParam(1:4) = param(2:5)
               a%MatProp(1)%LawViParam(5)= a%nu0_LCR
               a%Giesekus_model=1
            elseif(words(2)=='PTT')then          ! Phan-Thien-Tanner model            
               a%MatProp(1)%lawvi=-3  
               a%MatProp(1)%LawViParam(1:4) = param(2:5)  
               a%PTT_model=0
            end if 
         end if   
     
      elseif(a%nmat==2)then
         do imat=1,a%nmat
          if(words(1)=='MATER' .and. param(1)==imat)then
            call a%Listener%listen(outstr)
            do while(words(1)/='ENDMA')
               if(words(1)=='DENSI') then                 ! Density (rho)
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
                  elseif(words(2)=='OLDR')then             ! Oldroyd-B model  
                     a%MatProp(imat)%lawvi=-1  
                     a%MatProp(imat)%LawViParam(1:3) = param(2:4) 
                     a%MatProp(imat)%LawViParam(4) = 0.0_rp  
                  elseif(words(2)=='GIESE')then            ! Giesekus  
                     a%MatProp(imat)%lawvi=-2  
                     a%MatProp(imat)%LawViParam(1:4) = param(2:5) 
                     a%Giesekus_model=1
                  elseif(words(2)=='PTT')then          ! Phan-Thien-Tanner model            
                     a%MatProp(imat)%lawvi=-3  
                     a%MatProp(imat)%LawViParam(1:4) = param(2:5)
                     a%PTT_model=0
                  end if
                end if  
               call a%Listener%listen(outstr)
            end do
          end if
         end do
      elseif(a%nmat>2)then
         call runend('nsi_reaphy: first change the problem type to get more than 2 fluids')              
      end if  
         
   !Final operations      
   elseif(itask==100) then 
      
      if(a%PTT_model==1 .and. a%kfl_cotem_WLF==1)&
         call runend('WLF model only can be used with PTT viscoelastic model.')  
      

   endif
   
   aux_BC_number=a%kfl_bc_number
   

end subroutine sup_reaphy
