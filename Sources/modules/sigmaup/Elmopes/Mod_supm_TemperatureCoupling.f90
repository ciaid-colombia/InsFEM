module Mod_supm_TemperatureCoupling
   use Mod_supm_BaseElmope
   use Mod_supm_ExternalForces
   use Mod_nsi_BoussinesqForces
   use typre
   implicit none
   private
   public :: SetPointersTemperatureCoupling
   real(rp), allocatable :: eltem(:)
   real(rp) :: gptem(1), lambdaT
   integer(ip), allocatable :: kfl_IsSet 

contains

   subroutine SetPointersTemperatureCoupling(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            !Pointers are already set
            kfl_IsSet = 1

            if (a%kfl_cotem == 1) then
               call SetPointersExternalForces(1)
               call ConcatenateProcedures(ProcHook_Initializations,AllocTem)
               call ConcatenateProcedures(ProcHook_Gathers,GatherTem)
               call ConcatenateProcedures(ProcHook_Interpolates,InterpolateTem)
               if (a%kfl_ExternalTemperatureSGS == 1) then
                  call ConcatenateProcedures(ProcHook_Interpolates,InterpolateTemSGS)
               endif
               
               if (a%kfl_cotem_Boussinesq == 1) call ConcatenateProcedures(ProcPointer%ExternalForces_sup,BoussinesqForces)
               if (a%kfl_cotem_WLF==1)          call ConcatenateProcedures(ProcHook_Physicalprop,ViscoelasticTempParametersWLFmodel)
               if (a%kfl_cotem_Arrhenius==1)    call ConcatenateProcedures(ProcHook_Physicalprop,ViscoelasticTempParameters_ArrheniusModel)
               
               if (a%kfl_cotem_WLF==1 .or. a%kfl_cotem_Arrhenius==1) then 
                  !Update physical parameters
                  call ConcatenateProcedures(ProcHook_Physicalprop,UpdateViscoelasticAuxiliarParameter)
                  if (a%LogFormulation==1)call ConcatenateProcedures(ProcHook_Physicalprop,UpdateLogarithmicAuxiliarParameters)
                  !Print viscosity and lambda
                  if (a%npp_stepi(5)>=1)  call ConcatenateProcedures(ProcHook_Physicalprop,PlotViscosity)
                  if (a%npp_stepi(23)>=1) call ConcatenateProcedures(ProcHook_Physicalprop,PlotLambda)
               end if   

               
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocTem)
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-------------------------------------------------------------------
   !Computation Subroutines
   subroutine AllocTem
      implicit none
      call a%Memor%alloc(e%mnode,eltem,'eltem','nsm_elmope')
   end subroutine
   
   subroutine DeallocTem
      implicit none
      call a%Memor%dealloc(e%mnode,eltem,'eltem','nsm_elmope')
   end subroutine
   
   subroutine GatherTem
      implicit none
      call e%gather(1,eltem,a%tempe)
   end subroutine
   
   subroutine InterpolateTem
      implicit none
      call e%interpg(1,eltem,gptem(1))
  
   end subroutine   
   
   subroutine InterpolateTemSGS
      implicit none
      gptem = gptem + a%tesgs(ielem)%a(1,e%igaus)
   end subroutine
   
   subroutine BoussinesqForces
      implicit none
      call nsi_BoussinesqForces(e,acden,a%bougr,a%boube,a%boutr,gptem(1),elext)
   end subroutine  
   
   
   subroutine ViscoelasticTempParametersWLFmodel
      implicit none  
      call sup_templaw_WLF(e, gptem, a%ReferenceTemp, a%c1_WLF, a%c2_WLF, a%MatProp(imat)%LawViParam, acvis, lambdaT)
   end subroutine
   
   subroutine ViscoelasticTempParameters_ArrheniusModel
      implicit none
      call sup_templaw_Arrhenius(e, gptem, a%ReferenceTemp, a%alpha_Arrhenius, a%MatProp(imat)%LawViParam, acvis, lambdaT)
   end subroutine
   
   subroutine UpdateViscoelasticAuxiliarParameter
      implicit none
      call a%IncrementalLambda(imat,lambdaT,lambda)
      auxVE  = lambda/(2.0_rp*acvis)
   end subroutine   
   
   subroutine UpdateLogarithmicAuxiliarParameters
      implicit none
      if(a%kfl_timei==0) then
        lambda0=lambda*a%MatProp(1)%LawviParam(5)
      else  
        lambda0=a%MatProp(imat)%LawViParam(3)*a%MatProp(1)%LawviParam(5)
      end if  
            
      auxL = ((1.0_rp-beta)*acvis)/lambda0
      auxG   = a%MatProp(imat)%LawviParam(4)*lambda/(2.0_rp*(lambda0**2_ip) + 0.00001_rp) !GIESEKUSLCR    
   end subroutine
   
   subroutine PlotViscosity
      implicit none
      a%viscarray(ielem)%a(e%igaus) = acvis  
   end subroutine   
    
    subroutine PlotLambda
      implicit none
      a%lambdarray(ielem)%a(e%igaus) = lambda
    end subroutine  

end module