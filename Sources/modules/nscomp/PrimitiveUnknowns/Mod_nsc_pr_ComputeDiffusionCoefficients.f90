module Mod_nsc_pr_ComputeDiffusionCoefficients
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersDiffusionCoefficients
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersDiffusionCoefficients(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateDiffusionCoefficients)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,DiffusionCoefficientsToZero)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ComputeDiffusionCoefficients)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateDiffusionCoefficients)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateDiffusionCoefficients

      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Kmm,'Kmm','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Kem,'Kem','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Kee,'Kee','nsc_pr_elmope')

   end subroutine

   !Diffusion matrix coefficients
   subroutine DiffusionCoefficientsToZero

     !Stabilization Contribution                  
      Kmm = 0.0_rp !Kmm(i,d,p)
      Kem = 0.0_rp !Kem(d,p)
      Kee = 0.0_rp !Kee(p)

   end subroutine   

   subroutine ComputeDiffusionCoefficients
      implicit none    

      Kem = Kem + acvis*matmul(grvel,e%cartd)!mu*d_j(vel_d)*d_j
      Kem = Kem + acvis*matmul(transpose(grvel),e%cartd)!mu*d_d(vel_j)*d_j
      Kem = Kem - 2*acvis*divvel*e%cartd/3!-(2mu/3)*d_k(vel_k)*d_dj*d_j

   end subroutine   

   subroutine DeallocateDiffusionCoefficients

      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Kmm,'Kmm','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Kem,'Kem','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Kee,'Kee','nsc_pr_elmope')

   end subroutine

end module
