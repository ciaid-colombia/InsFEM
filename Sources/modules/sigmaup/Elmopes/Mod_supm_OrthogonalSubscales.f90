module Mod_supm_OrthogonalSubscales
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateGradients
   use Mod_supm_ComputeAdvectionVelocity
   use Mod_supm_InterpolateResidualProjection
   implicit none
   private
   public SetPointersOrthogonalSubscales
   
   integer(ip), allocatable :: kfl_IsSet

contains

 subroutine SetPointersOrthogonalSubscales(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1

            if (a%kfl_repro >= 1 .and. a%kfl_repro/=4) then
               call SetPointersInterpolateResidualProjection(1)
               call ConcatenateProcedures(ProcHook_Initializations,AllocGrpre)

               
               if (a%LogFormulation==0) then
                  if (a%kfl_repro==1) then   
                     call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepOSS)
                  elseif (a%kfl_repro==2) then
                     call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepSplit)
                  end if
                  
               elseif (a%LogFormulation==1) then
                  if (a%kfl_repro==1) then   
                     call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepOSS_LCR)
                  elseif (a%kfl_repro==2) then 
                     call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepSplitOSS_LCR)
                  elseif (a%kfl_repro==3) then 
                     call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepSplitASGS_LCR)
                  end if
               endif
               
              call ConcatenateProcedures(ProcHook_Finalizations,DeallocGpre)   
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select
      
   end subroutine  
   
   subroutine AllocGrpre
      implicit none
      call a%Memor%alloc(1_ip     ,e%ndime,grpre,'grpre','supm_elmope')
   end subroutine
   
   subroutine DeallocGpre
      implicit none
      call a%Memor%dealloc(1_ip     ,e%ndime,grpre,'grpre','supm_elmope')
   end subroutine

   !For ResidualProjection
   subroutine InGaussElmatsRepSplit
      implicit none     
      !Compute contribution to RHS in the blocks S y U
      call  ResProRHS
   end subroutine
   
   subroutine InGaussElmatsRepOSS
      implicit none     
      !Compute contribution to RHS in the blocks S y U
      call  ResProRHS
      !Compute contributions to RHS : Block P 
      !Only if OSS this term works
      call supm_elmrhp_oss(e,timom,dvol,gprep,auxtens,elrhp)
   end subroutine
   
   subroutine InGaussElmatsRepSplitASGS_LCR
      implicit none
      call supm_LCR_elmrhu_splitoss(e,acden,tidiv,dvol,gprep,auxtens,a%kfl_splitOSSMomentum,elrhu)
   end subroutine   
   
   
   subroutine InGaussElmatsRepSplitOSS_LCR
      implicit none
      call ResProRHS_LCR
   end subroutine   
   
   
   subroutine InGaussElmatsRepOSS_LCR
      implicit none
      call ResProRHS_LCR
      call supm_elmrhp_oss(e,timom,dvol,gprep,auxtens,elrhp)
   end subroutine   
   
   
   subroutine ResProRHS
      implicit none
      !If kfl_repro=2 (S-OSS) only tisig((1/2eta)*t+grad_sym(v),gprep) + tidiv*(div(v),gprep) 
      !all the terms are considered.
      !Compute contributions to RHS : Block C
      call supm_elmrhc_oss(e,timom,tisig,dvol,gprep,acvis,beta,auxtens,a%kfl_splitOSSMomentum,elrhc) 
      !Compute contributions to RHS : Block U
      call supm_elmrhu_oss(e,acden,tidiv,tisig,timom,dvol,AgradV,gprep,auxtens,a%kfl_splitOSSMomentum,elrhu) 
      !Viscoelastic case, remainder terms
      call supm_elmrhcVES_oss(e,auxVE,auxG,a%PTT_model,tisig,acvis,dvol,auxtens,AgradV_VE,grvel,gpsig,gprep,elrhc) 
   end subroutine  


   subroutine ResProRHS_LCR
      implicit none
      call supm_elmrhc_oss(e,timom,tisig,dvol,gprep,acvis*(1.0_rp-beta),0.0_rp,auxtens,a%kfl_splitOSSMomentum,elrhc)
      call supm_elmrhu_oss(e,acden,tidiv,tisig,timom,dvol,AgradV,gprep,auxtens,a%kfl_splitOSSMomentum,elrhu)
      call supm_LCR_elmrhc_oss(e,lambda0,lambda,beta,acvis,tisig,grvel,AgradV_VE,dvol,auxtens,gprep,elrhc)  
   end subroutine  

end module
