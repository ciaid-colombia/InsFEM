module Mod_nsc_pr_ComputeVariables
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersVariables
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersVariables(itask)
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
         
            !Non-linear subscales
            if (a%kfl_nolsg == 1) then
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeVariables,ComputeNonlinearVariables)
            endif

            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeVariables,ComputeAuxiliaryVariables)

            !State law
            if (a%lawde /= 0) then
               if (a%lawde == 1) then
                  call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeVariables,ComputeIdealAuxVariables)
               else if (a%lawde /= 1) then
                  call runend('Nsc_elmope: Non-ideal state law not ready')
               endif
            endif
      
         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   !Actual Computations

   subroutine ComputeNonlinearVariables
      implicit none
      
      gpadp = gpadp       + a%cosgs(ielem)%a(1,e%igaus)
      gpadv(:) = gpadv(:) + a%mosgs(ielem)%a(:,1,e%igaus)
      gpadt = gpadt       + a%ensgs(ielem)%a(1,e%igaus)

   end subroutine

   subroutine ComputeAuxiliaryVariables
      implicit none
      
      !Advection velocity norm 
      call vecnor(gpadv,e%ndime,gpvno,2)

   end subroutine

   ! Ideal Gas Law auxiliary variables
   subroutine ComputeIdealAuxVariables
      implicit none
    
      gpden = gpadp/((accph-accvh)*gpadt)
      aux = gpden*(accvh*gpadt+gpvno*gpvno/2.0_rp)+gpadp
      acalpha = 1.0_rp/gpadt
      acbeta = 1.0_rp/gpadp

   end subroutine
end module
