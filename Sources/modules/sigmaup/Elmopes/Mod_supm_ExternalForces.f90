 module Mod_supm_ExternalForces
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersExternalForces
   integer(ip), allocatable :: kfl_IsSet 
   real(rp) :: acceleration_module
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersExternalForces(itask)
      implicit none
      integer(ip) :: itask
      real(rp), external      :: funcre
               
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
             ProcPointer%ExternalForces_sup  => nsiForces      
            if (a%kfl_exacs/=0) then  
                ProcPointer%ExternalForces_sup => ExactSolutionForces
            endif 
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine  nsiForces
      implicit none    
      call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext)  
   end subroutine   
        
   subroutine  ExactSolutionForces    
      use Mod_Mesh    
      implicit none      
      real(rp)    :: gpcod(e%ndime), expre, exprg(e%ndime)
      call e%interpg(e%ndime,e%elcod,gpcod) 
      call exacso%sup_ComputeSolution(e%ndime,gpcod,a%ctime,a%LogFormulation,a)
      call exacso%sup_getPressure(e%ndime,expre,exprg)
      call exacso%sup_GetForce(e%ndime,a%LogFormulation,a%kfl_LogFormulation,elext,elextC,elextS,elextEstab,elextEstab2, elextSEstab,elextEstab3,&
                              elextSEstab2,elextEstab4,elextSEstab3,elextSEstab4,elextEstab5,elextSMat,elextSEstabMat,a)       
   end subroutine  
   
   
end module

