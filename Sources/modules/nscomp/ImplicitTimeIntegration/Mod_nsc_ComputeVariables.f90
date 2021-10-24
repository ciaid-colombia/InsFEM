module Mod_nsc_ComputeVariables
   use typre
   use Mod_nsc_BaseElmope
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
               call ConcatenateProcedures(ProcPointer_nsc_ComputeVariables,ComputeNonlinearVariables)
            endif

            call ConcatenateProcedures(ProcPointer_nsc_ComputeVariables,ComputeAuxiliaryVariables)

            !State law
            if (a%lawde /= 0) then
               if (a%lawde == 1) then
                  call ConcatenateProcedures(ProcPointer_nsc_ComputeVariables,ComputeIdealAuxVariables)
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
      
      gpadd = gpadd + a%cosgs(ielem)%a(1,e%igaus)
      gpadm(:) = gpadm(:) + a%mosgs(ielem)%a(:,1,e%igaus)
      gpade = gpade + a%ensgs(ielem)%a(1,e%igaus)

   end subroutine

   subroutine ComputeAuxiliaryVariables
      implicit none
      
      !Advection momentum norm 
      call vecnor(gpadm,e%ndime,gpmno,2)

      invcvh = 1.0_rp / accvh 
      invgpd = 1.0_rp / gpadd
      gpadv = gpadm / gpadd
      sqinvgpd = invgpd * invgpd
      sqgpmn = gpmno * gpmno
      sqgpvn = sqgpmn * sqinvgpd

   end subroutine

   ! Primitive auxiliary variables
   subroutine ComputeIdealAuxVariables
      implicit none
     
      acgamma = accph * invcvh
      actcn = invgpd * invcvh
      aux = acgamma - 1.0_rp
      aux_t = gpade - (sqgpmn*invgpd/2.0_rp) 
      aux_d = aux * sqgpvn / 2.0_rp
      gptem = invgpd * invcvh * aux_t 
      gppre = (acgamma - 1.0_rp) * aux_t
      
   end subroutine
end module
