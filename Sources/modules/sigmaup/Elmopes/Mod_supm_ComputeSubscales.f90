module Mod_supm_ComputeSubscales
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_SubgridSpaceResidual
   use Mod_supm_SubgridSpaceSplitTerms
   use Mod_supm_ComputeTaus
   implicit none
   private
   public SetPointersComputeSubscales_sup

   
   integer(ip), allocatable :: kfl_IsSet, kfl_IsSetGetSubscales
   real(rp)    ::  auxVEnew, auxlamb !For Logarithmic formulation
 
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeSubscales_sup(itask)
      implicit none
      integer(ip) :: itask
      integer(ip) :: kfl_nonlinear
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
      
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            !We need to compute the residual and project it 
            !to the subgrid scale space
            
            if ( a%kfl_tacsg == 1) then
               !We need to modify the momentum tau so that it includes 
               !the time step, but we keep the old timom in timom_static
               call SetPointersComputeTaus(1) 
               
               !We compute the transient subscales
               if(a%kfl_repro == 0 .or. a%kfl_repro == 1)then
                  call SetPointersComputeSubgridSpaceResidual(1) 
                  call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS)
                  if (a%MatProp(imat)%lawvi<0) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_VE)
                    
               elseif (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
                  call SetPointersComputeSubgridSpaceResidual(1) 
                  call SetPointersComputeSubgridSpaceSplitTerms(1)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split_mixed)  
                  if (a%LogFormulation==0) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split_mixed_vesgs3_STD) 
                  if (a%LogFormulation==1) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split_mixed_vesgs3_LCR)
                  if (a%MatProp(imat)%lawvi<0) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split_mixed_VE)

 
               elseif(a%kfl_repro == 4)then
                  call SetPointersComputeSubgridSpaceSplitTerms(1)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split)  
                  if (a%MatProp(imat)%lawvi<0) then 

                  if(a%LogFormulation==0) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split_VESTD) 
                  if(a%LogFormulation==1) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split_VELCR) 

                  end if

               end if
            endif
         endif
      
      case(100)
      
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   

   !-------------------------------------------------------------------
   !Compute 

   !Tracking of Subscales
   !Dynamic subscales
   subroutine ComputeSubgridScaleDSS
      implicit none
      auxtens=(e%ndime-1)*(e%ndime-1)+2    
      !Uses a Backward Euler scheme
      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = timom*(-gpSubscaleSpaceResidual(auxtens+1:auxtens+e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)

   end subroutine
   
   subroutine ComputeSubgridScaleDSS_VE
      implicit none
      auxVEnew=auxVE/(1.0_rp-beta*a%LogFormulation)
      a%sisgs(ielem)%a(1:auxtens,1,e%igaus) = tisig*(-gpSubscaleSpaceResidual(1:auxtens) + auxVEnew*a%sisgs(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
   end subroutine   
   
   
   !Dynamic subscales for Split OSS
   subroutine ComputeSubgridScaleDSS_Split_mixed
      implicit none
      !Uses a Backward Euler scheme
      !First subscale (convective component)
      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)  = timom*(-acden*gpSubscaleSpaceSplitV1(1:e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      !Second subscale (pressure gradient component)
      a%vesgs2(ielem)%a(1:e%ndime,1,e%igaus) = timom*(-gpSubscaleSpaceSplitV2(1:e%ndime)       + acden*a%vesgs2(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      
      if (a%LogFormulation==0) then
         !Third subscale (divergence of sigma)
         a%vesgs3(ielem)%a(1:e%ndime,1,e%igaus) = timom*(gpSubscaleSpaceSplitV3(1:e%ndime)       + acden*a%vesgs3(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      else if (a%LogFormulation==1) then
         a%vesgs3(ielem)%a(1:e%ndime,1,e%igaus) = timom*(auxL*gpSubscaleSpaceSplitV3(1:e%ndime)  + acden*a%vesgs3(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      end if
      
   end subroutine   
   
   subroutine ComputeSubgridScaleDSS_Split_mixed_vesgs3_STD
      implicit none
       !Third subscale (divergence of sigma)
       a%vesgs3(ielem)%a(1:e%ndime,1,e%igaus) = timom*(gpSubscaleSpaceSplitV3(1:e%ndime)       + acden*a%vesgs3(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
   end subroutine   
   
   subroutine ComputeSubgridScaleDSS_Split_mixed_vesgs3_LCR
      implicit none
       a%vesgs3(ielem)%a(1:e%ndime,1,e%igaus) = timom*(auxL*gpSubscaleSpaceSplitV3(1:e%ndime)  + acden*a%vesgs3(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
   end subroutine   
   
   
   subroutine ComputeSubgridScaleDSS_Split_mixed_VE
      implicit none
      auxVEnew=auxVE/(1.0_rp-beta*a%LogFormulation)
      a%sisgs(ielem)%a(1:auxtens,1,e%igaus) = tisig*(-gpSubscaleSpaceResidual(1:auxtens) + auxVEnew*a%sisgs(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
   end subroutine   
   
   
   
   !Dynamic subscales for Split OSS
   subroutine ComputeSubgridScaleDSS_Split
      implicit none
      !Uses a Backward Euler scheme
      !First subscale (convective component)
      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)  = timom*(-acden*gpSubscaleSpaceSplitV1(1:e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      !Second subscale (pressure gradient component)
      a%vesgs2(ielem)%a(1:e%ndime,1,e%igaus) = timom*(-gpSubscaleSpaceSplitV2(1:e%ndime)       + acden*a%vesgs2(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      !Third subscale (divergence of sigma)
      a%vesgs3(ielem)%a(1:e%ndime,1,e%igaus) = timom*(gpSubscaleSpaceSplitV3(1:e%ndime)       + acden*a%vesgs3(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      
   end subroutine   
   
   subroutine ComputeSubgridScaleDSS_Split_VESTD
      implicit none

      !First subscale (velocity gradient)
      a%sisgs(ielem)%a(1:auxtens,1,e%igaus)  = tisig*((1.0_rp-beta)*gpSubscaleSpaceSplitS1(1:auxtens) + auxVE*a%sisgs(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
      !Second subscale (convective term)
      a%sisgs2(ielem)%a(1:auxtens,1,e%igaus) = tisig*(-auxVE*gpSubscaleSpaceSplitS2(1:auxtens)        + auxVE*a%sisgs2(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
      !Third subscale (deformation terms)
      a%sisgs3(ielem)%a(1:auxtens,1,e%igaus) = tisig*(auxVE*gpSubscaleSpaceSplitS3(1:auxtens)        + auxVE*a%sisgs3(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
            
      
   end subroutine   
   
   subroutine ComputeSubgridScaleDSS_Split_VELCR
      implicit none
      auxVEnew=auxVE/(1.0_rp-beta)
      auxlamb=lambda/(2.0_rp*lambda0)
      
      !First subscale (velocity gradient)
      a%sisgs(ielem)%a(1:auxtens,1,e%igaus)  = tisig*(gpSubscaleSpaceSplitS1(1:auxtens)                 + auxVEnew*a%sisgs(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
      !Second subscale (convective term)
      a%sisgs2(ielem)%a(1:auxtens,1,e%igaus) = tisig*(-auxlamb*gpSubscaleSpaceSplitS2(1:auxtens)        + auxVEnew*a%sisgs2(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)
      !Third subscale (deformation terms)
      a%sisgs3(ielem)%a(1:auxtens,1,e%igaus) = tisig*(auxlamb*gpSubscaleSpaceSplitS3(1:auxtens)        + auxVEnew*a%sisgs3(ielem)%a(1:auxtens,2,e%igaus)*ReferenceDtinv)

   end subroutine   
   
    
end module 
