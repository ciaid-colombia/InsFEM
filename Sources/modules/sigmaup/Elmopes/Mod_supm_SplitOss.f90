 module Mod_supm_OrthogonalSplitTerms
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateSplitTermsProjection
   implicit none
   private
   public SetPointersOrthogonalSplitTherms
   integer(ip) :: kfl_nonlinear
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersOrthogonalSplitTherms(itask)
      implicit none
      integer(ip) :: itask
      !procedure() :: NULL()
      select case (itask)   
 
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1

            if (a%kfl_repro >= 2) then
               call SetPointersInterpolateSplitTermsProjection(1)
               if (a%LogFormulation==0) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsOssSplit_STD)
               if (a%LogFormulation==1) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsOssSplit_LCR)

               
               if (a%kfl_repro==4) call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsOssSplitConstitutive)
               if(kfl_nonlinear==1)then
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsNonLinearOssSplit)
               end if
            endif
               
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select
      
   end subroutine  

 
   
!--------------------------------------------------------------------------------------------------
!SUBROUTINES

   subroutine InGaussElmatsOssSplit_STD
      implicit none
         call supm_elmrhc_splitoss(e,auxtens,beta,timom,dvol,gpdivs,a%kfl_splitOSSMomentum,elrhc) 
         call InGaussElmatsOssSplitRHUandRHP
         
   end subroutine   
   
   
   subroutine InGaussElmatsOssSplit_LCR
      implicit none

         call supm_LCR_elmrhc_splitoss(e,auxtens,beta,auxL,timom,dvol,gpdivs,elrhc) 
         call InGaussElmatsOssSplitRHUandRHP
         
   end subroutine   
   
   subroutine InGaussElmatsOssSplitRHUandRHP
      implicit none
      
      call supm_elmrhu_splitoss(e,acden,timom,dvol,gpconv,AGradV,a%kfl_splitOSSMomentum,elrhu)
      call supm_elmrhp_splitoss(e,timom,dvol,gpgrap,a%kfl_splitOSSMomentum,elrhp)

   end subroutine
   
   subroutine InGaussElmatsOssSplitConstitutive
      implicit none

      call supm_elmrhu_splitoss_tidiv(e,tidiv,dvol,gpdivu,elrhu) 
   
      if (a%LogFormulation==0) then
      
         call supm_elmrhu_splitoss_tisig(e,tisig,beta,dvol,auxtens,gpgrau,elrhu) !beta=0 if three field model only
         
         if (a%MatProp(imat)%lawvi<0) then
            call supm_elmrhcVES_splitoss(e,auxVE,tisig,dvol,auxtens,AgradV_VE,grvel,gpconvsigma,gpdeform,elrhc)
         end if 
            
      else if (a%LogFormulation==1) then
      
         call supm_elmrhu_splitoss_tisig(e,tisig,0.0_rp,dvol,auxtens,gpgrau,elrhu)
         
         if (a%MatProp(imat)%lawvi<0) then
            call supm_elmrhcLCR_splitoss(e,auxVE,beta,lambda,lambda0,tisig,dvol,auxtens,AgradV_VE,grvel,gpconvsigma,gpdeform,elrhc)
         end if 
         
      end if

      
   end subroutine   
      
!----------------------------------------------------------------------------------------------------------
   !Laplacian term      
   subroutine InGaussElmatsNonLinearOssSplit
      implicit none
      
      call supm_elmrhu_splitossNonLinear(e,acvis,beta,timom,dvol,gplapl,a%kfl_splitOSSMomentum,elrhu)       
   end subroutine 
      
end module
