 module Mod_supm_ComputeTaus
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_ComputeTauSmoothing
   implicit none
   private
   public SetPointersComputeTaus, timom_static
   
   !Tau Smoothing
   real(rp), allocatable :: elSmoothedTau(:,:)
   
   !TransientSubgridScales
   real(rp) :: timom_static
   
   integer(ip), allocatable :: kfl_IsSet
   
   !For ALE
   logical :: isALE
   real(rp) :: gpvelnor
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersComputeTaus(itask)
      implicit none
      integer(ip) :: itask, step
      select case (itask)   
      
      case(0) 
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            !call a%Mesh%GetALE(isALE)
            if (a%kfl_timei==0) then
               step=a%itera
            else
               step=a%istep
            end if   
          
            !Taus are computed as usual
            if (a%kfl_Tausm == 0 .or. (step == 1)) then
               !Initialization
               tisig_static=0.0_rp
               
               if(a%MatProp(imat)%lawvi < 0)then
                  ProcPointer%TauConstitutive => TauElastic
                  if (a%LogFormulation==0) then
                     call ConcatenateProcedures(ProcPointer%TauConstitutive,TisigElastic_STD)
                  else if (a%LogFormulation==1) then
                     call ConcatenateProcedures(ProcPointer%TauConstitutive,TisigElastic_LCR)
                  end if   
               else
                  ProcPointer%TauConstitutive => TauNoElastic
               end if   
               
               if (a%npp_stepi(26)>=1) call ConcatenateProcedures(ProcPointer%TauConstitutive,PlotTaus)
            
            !Taus are interpolated from an elementary array
            elseif (a%kfl_Tausm == 1) then
            
               call ConcatenateProcedures(ProcHook_Initializations,AllocSmoothedTau)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocSmoothedTau)
               call ConcatenateProcedures(ProcHook_Gathers,GatherTau)
               call ConcatenateProcedures(ProcPointer%TauConstitutive,InterpolateTau)
               
            endif

            if (a%kfl_tacsg == 1) then
               call ConcatenateProcedures(ProcPointer%TauConstitutive,TransientTaus)
              if(a%MatProp(imat)%lawvi < 0) call ConcatenateProcedures(ProcPointer%TauConstitutive,TransientTaus_VE)
            endif
            
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !Compute Tau values
   subroutine TauElastic
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none

      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)      
      call sup_ComputeTidivVE(e,timom,a%staco,chale,tidiv)  

   end subroutine 
   
    subroutine TisigElastic_LCR
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none
      real(rp) :: auxG, w2,w3

      !GIESEKUSLCR, SET AUXG
      auxG   = (a%MatProp(imat)%LawviParam(4)*lambda)/(2.0_rp*(lambda0**2_ip) + 0.00001_rp)
      call sup_LCR_ComputeTisig(e,lambda,lambda0,auxG,beta,acvis,gpvno_VE,grvel,ExpGpPsi_Matrix,a%staco,chale,a%PTT_model,gpsig,w2,w3,tisig) 
      
   end subroutine 
   
   subroutine TisigElastic_STD
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none
      real(rp) :: auxG, w2,w3

      auxG   = a%MatProp(imat)%LawviParam(4)/((1.0_rp-beta) + 0.00001_rp)
      !GIESEKUSLCR, SET AUXG
      call sup_ComputeTiSig(e,lambda,acvis,gpvno_VE,grvel,a%staco,chale,auxG,a%PTT_model,gpsig,tisig) 
      
   end subroutine 
   
   subroutine PlotTaus
      implicit none
      a%tau_mom(ielem)%a(e%igaus)=timom
      a%tau_sig(ielem)%a(e%igaus)=tisig
   end subroutine   
   
   
   subroutine TauNoElastic
      use Mod_sup_ComputeTidivVE
      implicit none
      !Compute the stability parameters  
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)          
      if(a%kfl_colev==1)then
         call sup_ComputeTidivVE(e,timom,a%staco,chale,tidiv)
      else
         call sup_ComputeTidiv(e,acden,acvis,gpvno,a%staco(4),chale,tidiv)       
      end if
      !Stabilization Constituive parameter 
      call sup_ComputeTidiv(e,acden,acvis,gpvno,a%staco(4),chale,tisig)
      
   end subroutine     
   
   !To interpolate the tau values from smoothed array Tausmo
   subroutine AllocSmoothedTau
      implicit none
      
      call a%Memor%alloc(3,e%mnode,elSmoothedTau,'elSmoothedTau','nsm_EndElmope')
   end subroutine
   
   subroutine DeallocSmoothedTau
      implicit none
      
      call a%Memor%dealloc(3,e%mnode,elSmoothedTau,'elSmoothedTau','nsm_EndElmope')
   end subroutine
   
   subroutine GatherTau
      implicit none
      
      call e%gather(3,elSmoothedTau,a%Tausmo)
   end subroutine
   
   subroutine InterpolateTau
      implicit none
      real(rp) :: gptau(3), ig
      
      call e%interpg(3,elSmoothedTau,gptau)
     
      timom = gptau(1)
      tidiv = gptau(2)
      tisig = gptau(3)
   end subroutine
   
   !DynamicSS
   !Computes the transient stabilization parameter
   subroutine TransientTaus
      !timom is the transient one, the static one goes to timom_static
      timom_static = timom
      call ComputeTransientTau(a%dtinv,acden,timom_static,timom)  
   end subroutine
   
   subroutine TransientTaus_VE
      tisig_static=tisig
      call ComputeTransientTisig(a%dtinv,auxVE,tisig_static,tisig) 
   end subroutine
   
   
    subroutine ComputeTransientTisig(dtinv,auxVE,tisig_static,tisig)
      implicit none
      real(rp) :: auxVE, auxVEpol
      real(rp) :: dtinv, tisig_static, tisig

      auxVEpol=auxVE/(1.0_rp-beta*a%LogFormulation)
      tisig = 1.0_rp/(1.0_rp/tisig_static + auxVEpol*dtinv)
   end subroutine
   
end module

 
