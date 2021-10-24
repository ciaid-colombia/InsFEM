module Mod_supm_InterpolateSplitTermsProjection
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersInterpolateSplitTermsProjection

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersInterpolateSplitTermsProjection(itask)
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
         
            if (a%kfl_repro >=2) then
               call ConcatenateProcedures(ProcHook_Initializations,InitOssSplit)
               call ConcatenateProcedures(ProcHook_Finalizations,FinOssSplit)
               call ConcatenateProcedures(ProcHook_Gathers,GatherOssSplit)
               call ConcatenateProcedures(ProcHook_Interpolates,InterpolateOssSplit)

               
               if (a%kfl_repro==4)  then !Constitutive equation term-by-term
                  call ConcatenateProcedures(ProcHook_Initializations,InitOssSplitConstitutive)             
                  call ConcatenateProcedures(ProcHook_Finalizations,FinOssSplitConstitutive)
                  call ConcatenateProcedures(ProcHook_Gathers,GatherOssSplitConstitutive)
                  call ConcatenateProcedures(ProcHook_Interpolates,InterpolateOssSplitConstitutive) 
                  
                  if (a%LogFormulation==1) then
                     call ConcatenateProcedures(ProcHook_Initializations,InitOssSplitConstitutive_LCR)  
                     call ConcatenateProcedures(ProcHook_Finalizations,FinOssSplitConstitutive_LCR) 
                     call ConcatenateProcedures(ProcHook_Gathers,GatherOssSplitConstitutive_LCR)
                     call ConcatenateProcedures(ProcHook_Interpolates,InterpolateOssSplitConstitutive_LCR) 

                  end if
                  
                  if (a%MatProp(imat)%lawvi<0) then
                     call ConcatenateProcedures(ProcHook_Initializations,InitOssSplitConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_Finalizations,FinOssSplitConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_Gathers,GatherOssSplitConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_Interpolates,InterpolateOssSplitConstitutive_VE)
                  end if   
                  
               end if   
            end if
         end if
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
  
   subroutine InitOssSplit
      implicit none
      call a%Memor%alloc(e%ndime,gpdivs, 'gpdivs','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrepdivs, 'elrepdivs','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpconv, 'gpconv','supm_EniteElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrepconv, 'elrepconv','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpgrap, 'gpgrap','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elrepgrap, 'elrepgrap','supm_EnditeElmope')

   end subroutine
   
   
   subroutine InitOssSplitConstitutive
      implicit none
      
         call a%Memor%alloc(1,gpdivu, 'gpdivu','supm_EnditeElmope')
         call a%Memor%alloc(1,e%mnode,elrepdivu, 'elrepdivu','supm_EnditeElmope')
         call a%Memor%alloc(auxtens,gpgrau, 'gpgrau','supm_EnditeElmope')
         call a%Memor%alloc(auxtens,e%mnode,elrepgrau, 'elrepgrau','supm_EnditeElmope')

      
   end subroutine   
   
   subroutine InitOssSplitConstitutive_LCR
      implicit none

         call a%Memor%alloc(auxtens,gpexp, 'gpgexp','supm_elmope')
         call a%Memor%alloc(auxtens,e%mnode,elrepexp, 'elrepexp','supm_elmope')
      
   end subroutine 
   
   subroutine InitOssSplitConstitutive_VE
      implicit none
         call a%Memor%alloc(auxtens,gpconvsigma, 'gpconvsigma','supm_elmope')
         call a%Memor%alloc(auxtens,e%mnode,elrepconvsigma, 'elrepconvsigma','supm_elmope')
         call a%Memor%alloc(auxtens,gpdeform, 'gpdeform','supm_elmope')
         call a%Memor%alloc(auxtens,e%mnode,elrepdeform, 'elrepdeform','supm_elmope')  
   end subroutine   
   
   
    subroutine FinOssSplit
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(e%ndime,gpdivs, 'gpdivs','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepdivs, 'elrepdivs','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpconv, 'gpconv','supm_EniteElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepconv, 'elrepconv','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpgrap, 'gpgrap','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elrepgrap, 'elrepgrap','supm_EnditeElmope') 

   end subroutine
   
   subroutine FinOssSplitConstitutive
   
      call a%Memor%dealloc(1,gpdivu, 'gpdivu','supm_elmope')
      call a%Memor%dealloc(1,e%mnode,elrepdivu, 'elrepdivu','supm_elmope')
      call a%Memor%dealloc(auxtens,gpgrau, 'gpgrau','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elrepgrau, 'elrepgrau','supm_elmope')
      
   end subroutine  
   
   subroutine FinOssSplitConstitutive_LCR
      call a%Memor%dealloc(auxtens,gpexp, 'gpgexp','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elrepexp, 'elrepexp','supm_elmope')
   end subroutine  
   
   
   subroutine FinOssSplitConstitutive_VE
      implicit none

      call a%Memor%dealloc(auxtens,gpconvsigma, 'gpconvsigma','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elrepconvsigma, 'elrepconvsigma','supm_elmope')
      call a%Memor%dealloc(auxtens,gpdeform, 'gpdeform','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elrepdeform, 'elrepdeform','supm_elmope')          
      
   end subroutine   
         
   
    
   subroutine GatherOssSplit
      implicit none
      call e%gather(e%ndime,elrepdivs,a%reproSDiv)
      call e%gather(e%ndime,elrepconv,a%reproUGradU)
      call e%gather(e%ndime,elrepgrap,a%reproGradp)
      
     
   end subroutine
   
   subroutine GatherOssSplitConstitutive
      implicit none
      
      !S-OSS Term by term
      call e%gather(1,elrepdivu,a%reproDivU)
      call e%gather(auxtens,elrepgrau,a%reproGradU)

   end subroutine   
   
   
     subroutine GatherOssSplitConstitutive_LCR
      implicit none

      call e%gather(auxtens,elrepexp,a%reproExpS)
   end subroutine   
   
   
    subroutine GatherOssSplitConstitutive_VE
      implicit none
      call e%gather(auxtens,elrepconvsigma,a%reproUGradS)
      call e%gather(auxtens,elrepdeform,a%reproSGradU)
    end subroutine  
   

   subroutine InterpolateOssSplit
      implicit none
      !Interpolate
      call e%interpg(e%ndime,elrepdivs,gpdivs)
      call e%interpg(e%ndime,elrepconv,gpconv)
      call e%interpg(e%ndime,elrepgrap,gpgrap)
         
   end subroutine  
   
   subroutine InterpolateOssSplitConstitutive
      implicit none
      
      call e%interpg(1,elrepdivu,gpdivu)
      call e%interpg(auxtens,elrepgrau,gpgrau)
      

   end subroutine  
   
    subroutine InterpolateOssSplitConstitutive_LCR
      implicit none
     
      call e%interpg(auxtens,elrepexp,gpexp)

   end subroutine   
   
   subroutine InterpolateOssSplitConstitutive_VE
      implicit none

      call e%interpg(auxtens,elrepconvsigma,gpconvsigma)
      call e%interpg(auxtens,elrepdeform,gpdeform)
   end subroutine   
   
    
   
   subroutine AllocOssSplitNL
      implicit none
      call a%Memor%alloc(e%ndime,gplapl, 'gplapl','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,elreplapl, 'elreplapl','supm_elmope')    
   end subroutine

   subroutine DeallocOssSplitNL
      implicit none
      call a%Memor%dealloc(e%ndime,gplapl, 'gplapl','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elreplapl, 'elreplapl','supm_elmope')  
   end subroutine
   
   subroutine GatherOssSplitNL
      implicit none
      call e%gather(e%ndime,elreplapl,a%reproLapla)
   end subroutine

   subroutine InterpolateOssSplitNL
      implicit none
      call e%interpg(e%ndime,elreplapl,gplapl)
   end subroutine   

end module


!**********************************************************************************
! SUBGRID SPACE RESIDUAL 
!**********************************************************************************

module Mod_supm_SubgridSpaceSplitTerms
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_ComputeGpSplitTerms
   use Mod_supm_InterpolateSplitTermsProjection
   implicit none
   private
   public SetPointersComputeSubgridSpaceSplitTerms, gpSubscaleSpaceSplitV1
   public gpSubscaleSpaceSplitV2,  gpSubscaleSpaceSplitV3
   public gpSubscaleSpaceSplitS1, gpSubscaleSpaceSplitS2, gpSubscaleSpaceSplitS3 
   real(rp), allocatable :: gpSubscaleSpaceSplitV1(:), gpSubscaleSpaceSplitV2(:)
   real(rp), allocatable :: gpSubscaleSpaceSplitV3(:)
   real(rp), allocatable :: gpSubscaleSpaceSplitS1(:), gpSubscaleSpaceSplitS2(:)
   real(rp), allocatable :: gpSubscaleSpaceSplitS3(:)
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeSubgridSpaceSplitTerms(itask)
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
    
            if (a%kfl_repro >= 2) then
               !Allocations
               call ConcatenateProcedures(ProcHook_Initializations,AllocSGRes)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocSGres)
               
               !Compute the Residual at the Gauss point 
               call SetPointersComputeGpSplitTerms(1)
               call SetPointersInterpolateSplitTermsProjection(1)
               !Now we can obtain the projection of the residual onto the subscale space (gpres-gprep)
               call ConcatenateProcedures(ProcHook_InGaussElmats,SubgridSpaceSplitOSS)
               
               if (a%kfl_repro == 4 .and. a%MatProp(imat)%lawvi<0) then
                  call ConcatenateProcedures(ProcHook_Initializations,AllocSGResVE)
                  call ConcatenateProcedures(ProcHook_Finalizations,DeallocSGresVE)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,SubgridSpaceSplitOSS_VE)
               end if
               
            endif 
         endif 
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select   
   end subroutine   

   subroutine AllocSGRes
      implicit none
      
      call a%Memor%alloc(e%ndime,gpSubscaleSpaceSplitV1,'gpSubscaleSpaceSplitV1','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpSubscaleSpaceSplitV2,'gpSubscaleSpaceSplitV2','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpSubscaleSpaceSplitV3,'gpSubscaleSpaceSplitV3','supm_EnditeElmope')
   end subroutine
   
   subroutine DeallocSGRes
      implicit none
      
      call a%Memor%dealloc(e%ndime,gpSubscaleSpaceSplitV1,'gpSubscaleSpaceSplitV1','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpSubscaleSpaceSplitV2,'gpSubscaleSpaceSplitV2','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpSubscaleSpaceSplitV3,'gpSubscaleSpaceSplitV3','supm_EnditeElmope')
   end subroutine

   subroutine SubgridSpaceSplitOSS
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceSplitV1 = gpresconv- gpconv
      gpSubscaleSpaceSplitV2 = gpresgrap- gpgrap
      gpSubscaleSpaceSplitV3 = gpresdivs- gpdivs
      
   end subroutine
   
   
   subroutine AllocSGResVE
      implicit none
      
      call a%Memor%alloc(e%ndime,gpSubscaleSpaceSplitS1,'gpSubscaleSpaceSplitS1','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpSubscaleSpaceSplitS2,'gpSubscaleSpaceSplitS2','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpSubscaleSpaceSplitS3,'gpSubscaleSpaceSplitS3','supm_EnditeElmope')
   end subroutine
   
   subroutine DeallocSGResVE
      implicit none
      
      call a%Memor%dealloc(e%ndime,gpSubscaleSpaceSplitS1,'gpSubscaleSpaceSplitS1','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpSubscaleSpaceSplitS2,'gpSubscaleSpaceSplitS2','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpSubscaleSpaceSplitS3,'gpSubscaleSpaceSplitS3','supm_EnditeElmope')
   end subroutine

   subroutine SubgridSpaceSplitOSS_VE
      implicit none
      
      !Substract the residual projection to gpres
      gpSubscaleSpaceSplitS1=gpresgrau       - gpgrau
      gpSubscaleSpaceSplitS2=gpresconvsigma  - gpconvsigma
      gpSubscaleSpaceSplitS3=gpresdeform     - gpdeform
      
   end subroutine
   

end module



   
  
