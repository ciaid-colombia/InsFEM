 module Mod_supm_ElemVElastic
   use Mod_supm_BaseElmope
   implicit none
   private
   public setPointersElemVElastic
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine setPointersElemVElastic(itask,task)
      implicit none
      integer(ip) :: itask
      character(6)::task
      !procedure() :: NULL()
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            if(a%MatProp(imat)%lawvi < 0)then
               if (a%LogFormulation==0) then
                  if (task .eq. 'Elmope') then
                     ProcPointer%ViscoelasticGalerkin => elmatGalerk
                     ProcPointer%ViscoelasticEstab1 => elmatEstab1 
                     if (a%kfl_repro/=4) then
                        ProcPointer%ViscoelasticEstab2 => elmatEstab2 
                     else
                        ProcPointer%ViscoelasticEstab2 => elmatEstab2_split
                     end if
                     ProcPointer%DiscontinuityCapturing => discontinuity2d
                  else if (task .eq. 'Endite') then
                     ProcPointer%FEResPro   => FEResPro2dVE 
                  end if
               end if
            end if
           
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------
!SUBROUTINES      

   subroutine elmatGalerk
      !Constitutive
      call supm_elmbstGal(e,acvis,auxVE,auxG,a%PTT_model,LHSdtinv,dvol,auxtens,AgradV_VE,grvel,gpsig,elmst)   
      call supm_elmbutGal(e,auxVE,dvol,auxtens,beta,grvel,gpsig,grsig,elmut)
      call supm_elmrhcGal(e,auxVE,auxG,a%PTT_model,acvis,dvol,auxtens,elextS,gpvel,grvel,gpsig,grsig,elrhc)
      !Momentum
      call supm_elmbsvGal(e,dvol,auxtens,elmsv)
      call supm_elmbuvGal(e,dvol,acvis,beta,acden,AgradV,grvel,LHSdtinv,a%kfl_linearConvectiveTerm,elmuv)    
      call supm_elmbpvGal(e,dvol,a%LogFormulation,elmpv)     
      call supm_elmrhuGal(e,acden,dvol,elext,gpvel,grvel,a%kfl_linearConvectiveTerm,elrhu)
      !Continuity
      call supm_elmbuqGal(e,dvol,elmuq)       
      call supm_elmrhpGal(e,dvol,elextC,elrhp)
   end subroutine
  
  
   subroutine elmatEstab1
      implicit none
      
      !Constitutive
      call supm_elmbstEst1(e,beta,timom,dvol,auxtens,elmst)
      call supm_elmbutEst1(e,beta,acden,timom,LHSdtinv,dvol,auxtens,AgradV,grvel,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elmut)
      call supm_elmbptEst1(e,timom,dvol,auxtens,beta,a%kfl_splitOSSMomentum,elmpt)     
      call supm_elmrhcEst1(e,beta,acden,timom,dvol,elext,auxtens,gpvel,grvel,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elrhc) 
      !Momentum
      call supm_elmbsvEst1(e,timom,dvol,acden,AgradV,auxtens,a%kfl_splitOSSMomentum,elmsv)
      call supm_elmbuvEst1(e,dvol,beta,acden,timom,tidiv,AgradV,grvel,LHSdtinv,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elmuv)
      if (a%kfl_splitOSSMomentum/=0) call supm_elmbpvEst1(e,acden,timom,AgradV,dvol,elmpv)     
      call supm_elmrhuEst1(e,timom,tidiv,acden,auxtens,dvol,elextC,elext,gpvel,grvel,AgradV,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elrhu)     
      !Continuity
      call supm_elmbsqEst1(e,timom,dvol,auxtens,a%kfl_splitOSSMomentum,elmsq)
      call supm_elmbuqEst1(e,timom,dvol,acden,LHSdtinv,grvel,AgradV,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elmuq)
      call supm_elmbpqEst1(e,timom,dvol,elmpq)
      call supm_elmrhpEst1(e,acden,timom,dvol,elext,gpvel,grvel,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,elrhp)
      !only in splitt OSS with gravity         
      elext2=0.0_rp      
      !Compute vector of external forces
      call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext2)  
      call supm_elmrhpForce(e,a%kfl_splitOSSMomentum,timom,dvol,elext2,elrhp)
   end subroutine
   
   
    subroutine elmatEstab2 
      implicit none
      !Constitutive
      call supm_elmbstEst2(e,auxVE,auxG,a%PTT_model,acvis,tisig,a%kfl_reproTemporalTermZero*LHSdtinv,AgradV_VE,gpsig,dvol,auxtens,grvel,elmst)
      call supm_elmbutEst2(e,beta,auxVE,auxG,a%PTT_model,acvis,tisig,gpsig,dvol,auxtens,grvel,grsig,AgradV_VE,elmut)  
      call supm_elmrhcEst2(e,auxVE,auxG,a%PTT_model,tisig,acvis,gpsig,grsig,gpvel,AgradV_VE,grvel,dvol,auxtens,a%kfl_reproTemporalTermZero*elextS,elrhc) 
      !Momentum
      call supm_elmbsvEst2(e,beta,auxVE,auxG,a%PTT_model,acvis,tisig,a%kfl_reproTemporalTermZero*LHSdtinv,AgradV_VE,gpsig,dvol,auxtens,grvel,elmsv)
      call supm_elmbuvEst2(e,beta,auxVE,tisig,gpsig,dvol,auxtens,grvel,grsig,elmuv)
      call supm_elmrhuEst2(e,beta,auxVE,auxG,a%PTT_model,tisig,acvis,gpsig,grsig,gpvel,grvel,dvol,auxtens,a%kfl_reproTemporalTermZero*elextS,elrhu)  
   end subroutine
   
   
   subroutine elmatEstab2_split
      implicit none
      !This subroutine is used for kfl_repro==4 only (Split-OSS + Dynamic)
      !PTT and Giesekus models not implemented yet
      !Constitutive
      call supm_elmbstEst2_split(e,auxVE,tisig,AgradV_VE,gpsig,dvol,auxtens,grvel,elmst) 
      call supm_elmbutEst2_split(e,auxVE,tisig,gpsig,dvol,auxtens,grvel,grsig,AgradV_VE,elmut) 
      call supm_elmrhcEst2_split(e,auxVE,tisig,gpsig,grsig,gpvel,AgradV_VE,grvel,dvol,auxtens,elrhc) 
      !Momentum
      call supm_elmbuvEst2_split(e,beta,tisig,gpsig,dvol,auxtens,grvel,grsig,elmuv) 
   end subroutine
 
   
   subroutine discontinuity2d
      implicit none   
      call supm_elmbstDC(e,auxtens,kdisc,dvol,elmst)   
   end subroutine
   
   
   subroutine FEResPro2dVE
      implicit none   
      call supm_elmrfeVE_oto2(e,LHSdtinv,a%kfl_reproTemporalTermZero*LHSdtinv,acden,acvis,beta,lambda,auxVE,auxG,a%PTT_model,gpadvec, &
         gpvel,gpsig,grpre,grvel,grsig,auxtens,elext,elextC,a%kfl_reproTemporalTermZero*elextS,gpres)
   end subroutine
   
   
end module
