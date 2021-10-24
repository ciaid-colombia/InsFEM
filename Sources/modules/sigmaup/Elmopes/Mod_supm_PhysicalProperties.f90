module Mod_supm_PhysicalProperties
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_TemporalDerivatives
   use Mod_supm_InterpolateGradients
   use Mod_nsm_Viscosity
   implicit none
   private
   public SetPointersPhysicalPropertiesSUP
   integer(ip), allocatable :: kfl_IsSet
   real(rp) :: lambdaT
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersPhysicalPropertiesSUP(itask,task)
      integer(ip) :: itask
      character(6) :: task
      !procedure() :: NULL()
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            !Viscosity law
            if (task .eq. 'Elmope') then
               if(a%MatProp(imat)%lawvi > 0)then      
                  call ConcatenateProcedures(ProcHook_PhysicalProp,ViscosityLawElmope)
               endif
               if (a%MatProp(imat)%lawvi < 0) then !viscoelastic fluid
                  ProcPointer%TwoFieldTerms   => NULLSUB       
                  call SetPointersTemporalDerivatives(1)
                  call ConcatenateProcedures(ProcHook_Physicalprop,ViscoelasticParameters)    
                  if(a%LogFormulation==1) call ConcatenateProcedures(ProcHook_Physicalprop,LogarithmicReformulationParameters)

                  ProcPointer%ConstitutiveComponents => elementalVE       
               elseif(a%MatProp(imat)%lawvi>=0)then !only Viscous fluid
                  beta=0.0_rp
                  lambda=0.0_rp
                  auxG=0.0_rp
                  auxVE=0.0_rp      
               endif    
      
            elseif (task .eq. 'Endite') then
               if(a%MatProp(imat)%lawvi > 0)then
                  call ConcatenateProcedures(ProcHook_PhysicalProp,ViscosityLawEnditeElmope)
               endif
               
               if (a%MatProp(imat)%lawvi < 0) then !viscoelastic fluid 
                  call ViscoelasticParameters 
                  if (a%LogFormulation==1) call LogarithmicReformulationParameters
               elseif (a%MatProp(imat)%lawvi>=0)then !only Viscous fluid
                  beta=0.0_rp
                  lambda=0.0_rp
                  auxG=0.0_rp
                  auxVE=lambda/(2.0_rp*acvis)      
               endif             
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
      
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !-------------------------------------------------------------------
   !Physical Properties
   subroutine ViscosityLawElmope 
      implicit none   
      call nsi_vislaw(e,grvel,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis)  
      if ((a%kfl_repro == 1).or.(a%npp_stepi(5)==1)) then
         a%viscarray(ielem)%a(e%igaus) = acvis
      endif
   end subroutine
   
   subroutine ViscoelasticParameters
      implicit none 
       
      if (a%kfl_cotem_WLF==0 .and. a%kfl_cotem_Arrhenius==0)  then
         acvis  = a%MatProp(imat)%LawViParam(1) 
         !incremental scheme for Weissenberg number
         call a%IncrementalLambda(imat,a%MatProp(imat)%LawViParam(3), lambda)  
         !a%lambdarray(ielem)%a(e%igaus) = lambda
         auxVE  = lambda/(2.0_rp*acvis)
      endif
      beta    = a%MatProp(imat)%LawViParam(2)
      auxG    = a%MatProp(imat)%LawviParam(4)/((1.0_rp-beta) + 0.00001_rp)
   end subroutine
   
   subroutine LogarithmicReformulationParameters
      implicit none
      !Logarithmic Reformulation not coupled with tempe
      if (a%kfl_cotem_WLF==0 .and. a%kfl_cotem_Arrhenius==0) then
         !Definition of lambda0
         if (a%MatProp(imat)%LawviParam(3) == 0 .or. a%MatProp(imat)%LawviParam(5) == 0) then
            lambda0=0.01_rp
         else   
            if (a%kfl_exacs/=0) then
               lambda0=a%MatProp(imat)%LawviParam(5) !SETLAMBDA0
            else
               if(a%kfl_timei==0) then
                  lambda0=lambda*a%MatProp(1)%LawviParam(5)
               else  
                  lambda0=a%MatProp(imat)%LawViParam(3)*a%MatProp(1)%LawviParam(5)
               end if  
            end if   
         end if   
         auxL = ((1.0_rp-beta)*acvis)/lambda0
         auxG   = a%MatProp(imat)%LawviParam(4)*lambda/(2.0_rp*(lambda0**2_ip) + 0.00001_rp) !GIESEKUSLCR 
      end if
      
   end subroutine   
   
   subroutine ViscosityLawEnditeElmope 
      implicit none              
      acvis = a%viscarray(ielem)%a(igaus) 
   end subroutine 
   
   subroutine elementalVE 
      !Galerking terms      
      call  ProcPointer%ViscoelasticGalerkin 
      !Momentum and continuity stabilization terms 
      call  ProcPointer%ViscoelasticEstab1  
      !Constitutive stabilization terms
      call  ProcPointer%ViscoelasticEstab2 
   end subroutine
  
end module 
