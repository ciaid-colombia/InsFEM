 module Mod_supm_NonLinearElem
   use Mod_supm_BaseElmope
   use Mod_supm_InGaussNonLinear
   implicit none
   private
   public SetPointersNonLinearElem
   integer(ip) :: kfl_nonlinear
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersNonLinearElem(itask,task)
      implicit none
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
            if (task .eq. 'Elmope') then

               if(a%MatProp(imat)%lawvi>=0)then
                  call ConcatenateProcedures(ProcHook_InGaussElmats,ProcPointer%PostGaussElmats_sup)  
               end if   
               
               ProcPointer%PostGaussElmats_sup => NULLSUB

               if(a%MatProp(imat)%lawvi<0 )then !.and. a%kfl_splitOSSMomentum==1
                  if (a%LogFormulation==0) then
                    ProcPointer%ViscoelasticEstab1Lapla => elmatEstab1Lapla
                  else if (a%LogFormulation==1) then
                    ProcPointer%ViscoelasticEstab1Lapla => elmatEstab1Lapla_LCR
                  end if  
                  call ConcatenateProcedures(ProcHook_HessianComponents,InGaussElmatsNonLinear) 
               end if

            else if (task .eq. 'Endite') then
               if (a%MatProp(imat)%lawvi<0 ) then
                  call ConcatenateProcedures(ProcHook_InGauss,InGaussNonLinear_Endite)
               endif
            endif  
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------
!SUBROUTINES
   subroutine InGaussNonLinear_Endite
      implicit none
      call e%elmder
      call e%elmhes
   end subroutine
   
   subroutine InGaussElmatsNonLinear
      implicit none 
      call  ProcPointer%ViscoelasticEstab1Lapla   
   end subroutine
   
   
   subroutine elmatEstab1Lapla     
      !Constitutive
      call supm_elmbutEst1Lapla(e,beta,acvis,timom,dvol,auxtens,a%kfl_splitOSSMomentum,elmut) 
      
      !Momentum
      call supm_elmbsvEst1Lapla(e,beta,acvis,timom,dvol,auxtens,a%kfl_splitOSSMomentum,elmsv)
      call supm_elmbuvEst1Lapla(e,beta,acvis,acden,timom,dvol,auxtens,Agradv,grvel,LHSdtinv,a%kfl_splitOSSMomentum,elmuv)
      call supm_elmbpvEst1Lapla(e,beta,acvis,timom,dvol,a%kfl_splitOSSMomentum,elmpv)
      call supm_elmrhuEst1Lapla(e,beta,timom,acden,acvis,dvol,elext,grvel,gpvel,a%kfl_splitOSSMomentum,elrhu)
      
      !Continuity
      call supm_elmbuqEst1Lapla(e,beta,acvis,timom,dvol,a%kfl_splitOSSMomentum,elmuq)
      if(a%kfl_repro==1)then
        call supm_elmrhu_ossNL(e,auxtens,acvis,beta,timom,dvol,gprep,a%kfl_splitOSSMomentum,elrhu)
      end if      
   end subroutine


   
   subroutine elmatEstab1Lapla_LCR    
      implicit none
      
      !Constitutive
      call supm_elmbutEst1Lapla(e,beta,acvis,timom,dvol,auxtens,a%kfl_splitOSSMomentum,elmut) 
      
      !Momentum
      call supm_LCR_elmbsvEst1Lapla(e,beta,acvis,auxL,timom,dvol,auxtens,a%kfl_splitOSSMomentum,DivExpPsi,ExpGpPsi_Matrix,elmsv)
      call supm_elmbuvEst1Lapla(e,beta,acvis,acden,timom,dvol,auxtens,Agradv,grvel,LHSdtinv,a%kfl_splitOSSMomentum,elmuv)
      call supm_elmbpvEst1Lapla(e,beta,acvis,timom,dvol,a%kfl_splitOSSMomentum,elmpv)
      call supm_LCR_elmrhuEst1Lapla(e,beta,auxtens,timom,acden,acvis,dvol,elext,grvel,gpvel,a%kfl_splitOSSMomentum,a%kfl_linearConvectiveTerm,RHSmom,elrhu) 
      
      !Continuity
      call supm_elmbuqEst1Lapla(e,beta,acvis,timom,dvol,a%kfl_splitOSSMomentum,elmuq)
      
      if(a%kfl_repro==1)then
       call supm_elmrhu_ossNL(e,auxtens,acvis,beta,timom,dvol,gprep,a%kfl_splitOSSMomentum,elrhu)
      end if      
   end subroutine
      
   
end module
