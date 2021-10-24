module Mod_supm_ElemNoVElastic
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateResidualProjection
   implicit none
   private
   public SetPointersNoVElastic
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersNoVElastic(itask, task)
      implicit none
      integer(ip) :: itask
      character(6) :: task
      !!procedure() :: NULL() 
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1

            if (task .eq. 'Elmope') then

               !if (a%kfl_tacsg==1 .and. (a%kfl_repro==4 .or. a%kfl_repro==2)) then !DynamicSS - SOSS
                if (a%kfl_repro==4) then
                  !Pointers inherited from the NSTINC module
                  ProcPointer%supm_elmbuv => nsm_elmbuv_split
                  ProcPointer%supm_elmrhu => nsm_elmrhu_split
                  ProcPointer%supm_elmrhp => nsm_elmrhp_split
                  ProcPointer%supm_elmbuq => nsm_elmbuq_split
                  ProcPointer%supm_elmbpv => nsm_elmbpv_split
               end if
               
               ProcPointer%TwoFieldTerms          => elmatTwoFieldUP 
               
               ProcPointer%ConstitutiveComponents => elementalMatrixNonElastic
               
            else if (task .eq. 'Endite') then  !Compute Residual Three Field Problem
      
               ProcPointer%FEResPro   => FEResPro2d
            end if   
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------
!SUBROUTINES      

  !Equations implemented in nstinc module and inherited by sigmaup
  subroutine elmatTwoFieldUP
      real(rp) :: eltemp(e%ndime)
      eltemp = 0.0_rp
       
      !Compute contributions to RHS : Block U
      call  ProcPointer%supm_elmrhu(e,dvol,testf,elext,eltemp,elrhu)               
      !Compute contributions to elemental matrix : Block U,V
      call  ProcPointer%supm_elmbuv(e,dvol,acden,LHSdtinv,AGradV,testf,wrmat1)         
      !Compute contributions to elemental matrix : Block V,P
      call  ProcPointer%supm_elmbpv(e,dvol,testf,elmpv)            
      !Compute contributions to RHS : Block P
      call  ProcPointer%supm_elmrhp(e,timom,dvol,elext,eltemp,elrhp)            
      !Compute contributions to elemental matrix : Block U,Q
      call  ProcPointer%supm_elmbuq(e,timom,dvol,acden,LHSdtinv,AGradV,elmuq)
  end subroutine
   
  !Equations that belong to the non elastic case
  subroutine elementalMatrixNonElastic
      implicit none      
      !Compute contribution to elemental matrix : Block S,T
      call supm_elmbst(e,timom,tisig,dvol,acvis,auxtens,a%kfl_splitOSSConstitutive,elmst)
      
      !Compute contribution to elemental matrix : Block U,T
      call supm_elmbut(e,timom*a%kfl_splitOSSMomentum,tisig,LHSdtinv,AGradV,dvol,acden,acvis,auxtens,a%kfl_splitOSSConstitutive,elmut)
            
      !Compute contribution to elemental matrix : Block P,T
      call supm_elmbpt(e,timom*a%kfl_splitOSSMomentum,dvol,auxtens,elmpt) 
      !Compute contribution to elemental matrix : RHConstitutive
      call supm_elmrhc(e,timom*a%kfl_splitOSSMomentum,tisig,acvis,dvol,elext,auxtens,elextS,elrhc)  
      
      !Compute contribution to elemental matrix : Block S,V
      call supm_elmbsv(e,timom*a%kfl_splitOSSMomentum,tisig,dvol,acvis,acden,AGradV,auxtens,a%kfl_splitOSSConstitutive,elmsv)
      !Compute contribution to elemental matrix : Block U,V (Grad_Sym,Grad_Sym) + beta*acvis(Grad_v,Grad_u) 
      call supm_elmbuv2(e,tisig,dvol,wrmat2)
      

      !Compute contribution to elemental matrix : RHmomentum
      call supm_elmrhu(e,tidiv,tisig,auxtens,dvol,elextC,elextS,elrhu)             
      !Compute contribution to elemental matrix : Block S,Q
      call supm_elmbsq(e,timom*a%kfl_splitOSSMomentum,dvol,auxtens,elmsq)       
      !Compute contribution to elemental matrix : RHcontinuity
      call supm_elmrhp(e,dvol,elextC,elrhp)  
   
   end subroutine
   
   
   subroutine FEResPro2d
      implicit none      
      call supm_elmrfe_oto(e,LHSdtinv,acden,acvis,gpadvec, &
        gpvel,gpsig,grpre,grvel,grsig,auxtens,elext,gpres)
   end subroutine   
      
end module
