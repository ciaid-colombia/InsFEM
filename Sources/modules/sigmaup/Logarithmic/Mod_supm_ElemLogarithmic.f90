 module Mod_supm_ElemLogarithmic
   use Mod_supm_BaseElmope
   use Mod_Debugging
   implicit none
   private
   public SetPointersElemLCR
   integer(ip), allocatable :: kfl_IsSet
   type(AuxTimer) :: TimerLCR
contains
   
   subroutine SetPointersElemLCR(itask,task)
      implicit none
      integer(ip) :: itask
      character(6)::task
      !!procedure() :: NULL()
      select case (itask)
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            if(a%LogFormulation==1) then
               if (task .eq. 'elmope') then
                  ProcPointer%ViscoelasticGalerkin => elmatGalerk_LCR
                  ProcPointer%ViscoelasticEstab1 => elmatEstab1_LCR_Split
                  ProcPointer%ViscoelasticEstab2 => elmatEstab2_LCR
                  
                  if (a%kfl_splitOSSMomentum/=0) then
                     !Add Cross terms
                     call ConcatenateProcedures(Procpointer%ViscoelasticEstab1,elmatEstab1_LCR_Cross)
                  end if
                  
               else if (task .eq. 'Endite') then
                  ProcPointer%FEResPro   => FEResPro2dLCR
               end if
               
            end if  
         end if   
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------

!SUBROUTINES      
    subroutine elmatGalerk_LCR
      implicit none
       
      !Momentum
      call supm_LCR_elmbsvGal(e,dvol,auxtens,auxL,DivExpPsi,ExpGpPsi_Matrix,elmsv)
      call supm_elmbuvGal(e,dvol,acvis,beta,acden,AGradV,grvel,LHSdtinv,a%kfl_linearConvectiveTerm,elmuv)        
      call supm_LCR_elmbpvGal(e,dvol,auxL,elmpv)
      call supm_LCR_elmrhuGal(e,auxtens,auxL,acden,dvol,elext,grvel,gpvel,gpsig,ExpGpPsi_Matrix,RHS1momGal,RHS2momGal,elrhu)
      
      !Constitutive 
       call supm_LCR_elmbstGal(e,lambda,lambda0,auxG,LHSdtinv,dvol,auxtens,AGradV_VE,gpvel,grvel,&
                                 ExpGpPsi_Matrix,ExpGiesekus_Matrix,GrPsiMatrix,GrExpMatrix,ConvExpMatrix,GrvelExp,elmst)
       call supm_LCR_elmbutGal(e,dvol,lambda,lambda0,auxtens,gpsig,ExpGpPsi_Matrix,GrExpMatrix,a%kfl_linearConstitutiveTerms,elmut)
       call supm_LCR_elmrhcGal(e,auxG,dvol,auxtens,elextS,ExpGiesekus_Matrix,ExpGpPsi_Matrix,gpsig,RHSconst,RHSGiesekus,elrhc)
                                 
      !Continuity
      call supm_elmbuqGal(e,dvol,elmuq)       
      call supm_elmrhpGal(e,dvol,elextC,elrhp)
   end subroutine
   
   
   subroutine elmatEstab1_LCR_Cross
   implicit none
      !Momentum
      call supm_LCR_elmbsvEst1(e,timom,dvol,acden,auxL,AGradV,DivExpPsi,ExpGpPsi_Matrix,auxtens,elmsv)
      call supm_elmbpvEst1(e,acden,timom,AGradV,dvol,elmpv) 
      
      !Constitutive
      call supm_LCR_elmbutEst1(e,acden,timom,LHSdtinv,dvol,auxtens,a%kfl_linearConvectiveTerm,AGradV,grvel,elmut)   
      call supm_LCR_elmbptEst1(e,timom,dvol,auxtens,elmpt)
      
      !Continuity
      call supm_LCR_elmbsqEst1(e,auxL,timom,dvol,auxtens,DivExpPsi,ExpGradV,elmsq)
      call supm_elmbuqEst1(e,timom,dvol,acden,LHSdtinv,grvel,AGradV,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elmuq) 
   end subroutine 
  
  
  subroutine elmatEstab1_LCR_Split
      implicit none
      !Momentum
      call supm_elmbuvEst1(e,dvol,beta,acden,timom,tidiv,AGradV,grvel,LHSdtinv,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elmuv) 
      call supm_LCR_elmrhuEst1(e,timom,tidiv,acden,auxL,auxtens,dvol,elextC,elext,gpvel,grvel,AGradV,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,RHSmom,elrhu)

      !Constitutive
      call supm_LCR_elmbstEst1(e,timom,dvol,auxL,gpsig,auxtens,DivExpPsi,ExpGpPsi_Matrix,elmst)  
      call supm_LCR_elmrhcEst1(e,acden,timom,dvol,elext,auxtens,gpvel,grvel,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,RHSmom,elrhc)
                              
      !Continuity 
      call supm_elmbpqEst1(e,timom,dvol,elmpq)
      call supm_LCR_elmrhpEst1(e,auxtens,acden,timom,dvol,elext,gpvel,grvel,gpsig,a%kfl_linearConvectiveTerm,a%kfl_splitOSSMomentum,RHSmom,elrhp)
                               
      !only in splitt OSS with gravity          
      elext2=0.0_rp      
      !Compute vector of external forces
      call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext2)   
      call supm_elmrhpForce(e,a%kfl_splitOSSMomentum,timom,dvol,elext2,elrhp)
   end subroutine
   
   subroutine elmatEstab2_LCR   
      implicit none

      !Constitutive 
      call supm_LCR_elmbstEst2(e,acvis,beta,lambda,lambda0,auxG,tisig,LHSdtinv*a%kfl_reproTemporalTermZero,AGradV_VE,gpvel,grvel,gpsig,&
         ExpGpPsi_Matrix,ExpGiesekus_Matrix,GrExpMatrix,ConvExpMatrix,GrvelExp,dvol,auxtens,a%kfl_splitOSSConstitutive,elmst) 
      call supm_LCR_elmbutEst2(e,lambda,lambda0,auxG,acvis,beta,tisig,gpsig,grvel,dvol,auxtens,ExpGpPsi_Matrix,GrExpMatrix,AGradV_VE,&
         a%kfl_splitOSSConstitutive,a%kfl_linearConstitutiveTerms,elmut) 
      
!       if (a%kfl_repro/=4) then
         call supm_LCR_elmrhcEst2(e,lambda0,lambda,auxG,beta,acvis,tisig,gpsig,grsig,grvel,AGradV_VE,ExpGpPsi_Matrix,ExpGiesekus_Matrix,&
            dvol,auxtens,elextS*a%kfl_reproTemporalTermZero,RHSconstEstab,RHSGiesekus,elrhc)
!       else
!          call supm_LCR_elmrhcEst2_split(e,lambda0,lambda,beta,acvis,tisig,grvel,AGradV_VE, dvol,auxtens,RHSconst_conv,RHSconst_deform,elrhc)         
!       end if

      !Momentum 
      call supm_LCR_elmbsvEst2(e,auxtens,lambda,lambda0,auxG,LHSdtinv*a%kfl_reproTemporalTermZero,tisig,gpvel,grvel,AGradV_VE,dvol,&
         ExpGpPsi_Matrix,ExpGiesekus_Matrix,GrExpMatrix,ConvExpMatrix,a%kfl_splitOSSConstitutive,elmsv) !
      call supm_LCR_elmbuvEst2(e,lambda,lambda0,tisig,gpsig,dvol,auxtens,ExpGpPsi_Matrix,GrExpMatrix,a%kfl_splitOSSConstitutive,&
         a%kfl_linearConstitutiveTerms,elmuv)
      call supm_LCR_elmrhuEst2(e,auxG,tisig,gpsig,ExpGiesekus_Matrix,ExpGpPsi_Matrix,dvol,auxtens,elextS*a%kfl_reproTemporalTermZero,&
         RHSconstEstab*a%kfl_splitOSSConstitutive,RHSGiesekus,elrhu) 
   end subroutine
   
   
   subroutine FEResPro2dLCR
      implicit none 
      call supm_elmrfeLCR_oto(e,LHSdtinv,LHSdtinv*a%kfl_reproTemporalTermZero,acden,acvis,lambda,lambda0,auxG,auxL,gpadvec,&
           gpvel,gpsig,grsig,grpre,grvel,ExpGpPsi_Matrix,GrExpMatrix,auxtens,elext,elextC,a%kfl_reproTemporalTermZero*elextS,gpres)
   end subroutine
   
end module
