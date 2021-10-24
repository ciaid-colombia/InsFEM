 module Mod_nsc_pr_ComputeTaus
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersComputeTaus, timom_static
   
   !External Procedures
   procedure() :: NULLSUB

   !Dynamic subscales
   real(rp) :: timom_static(3)
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersComputeTaus(itask)
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
            
            !Taus are computed as usual
            ProcHook_nsc_pr_ComputeTaus => ComputeTaus
            
            if ( a%kfl_tacsg == 1) then
               call ConcatenateProcedures(ProcHook_nsc_pr_ComputeTaus,TransientTaus)
            endif

         endif  

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   !Compute Tau values
   subroutine ComputeTaus
     implicit none
     real(rp)                   :: auxtm, auxte
      
     call a%ComputeNscompNstincVelocity(gpvno,gpspd,gpvst)
     call ComputeTau(e,gpden,acvis,gpvst,a%staco,chale,auxtm)
     call ComputeTauCDR(e,gpden*accvh,actco,gpden*a%grnor,gpvst,a%staco,chale,auxte) 

     timom(2) = auxtm
     timom(3) = auxte
     timom(1) = a%staco(4)*chale(2)*chale(2)/(gpden*auxtm*e%npol*e%npol)

   end subroutine
   
   !Computes the transient stabilization parameter
   subroutine TransientTaus
     implicit none
     timom_static = timom
     !timom is the transient one, the steady one goes to timom_static
     call ComputeTransientTau(a%dtinv,(gpden*chale(2))/((gpspd+gpvno)*acvis),timom_static(1),timom(1)) 
     call ComputeTransientTau(a%dtinv,gpden,timom_static(2),timom(2))
     call ComputeTransientTau(a%dtinv,gpden*accph,timom_static(3),timom(3))
   end subroutine

end module

 
