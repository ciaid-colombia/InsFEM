module Mod_nsc_pr_ComputeNonLinearDiffusionCoefficients
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersNonLinearDiffusionCoefficients
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersNonLinearDiffusionCoefficients(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateHighOrderDerivatives)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ComputeHessian)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ComputeNonLinearDiffusionCoefficients)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateHighOrderDerivatives)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateHighOrderDerivatives
      
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,hess,'hess','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,lapl,'lapl','nsc_pr_elmope')

   end subroutine

   !----------------------------------------------------------
   subroutine DeallocateHighOrderDerivatives
      
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,hess,'hess','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,lapl,'lapl','nsc_pr_elmope')

   end subroutine

   subroutine ComputeHessian

      !Hessian and Laplacian
      !hess(i,j,p)
      !lapl(p)
      call nsc_HessLap(e,hess,lapl)

   end subroutine

   !Diffusion matrix coefficients
   subroutine ComputeNonLinearDiffusionCoefficients
      implicit none    
      integer(ip)    :: idime,jdime


     !Momentum equation
      do idime=1,e%ndime

         Kmm(idime,idime,:) = Kmm(idime,idime,:) + acvis*lapl(:)!mu*d_id*d_kj*d^2_kj
         do jdime=1,e%ndime
            Kmm(idime,jdime,:) = Kmm(idime,jdime,:) + acvis*hess(jdime,idime,:)!mu*d_ij*d_kd*d^2_kj
            Kmm(idime,jdime,:) = Kmm(idime,jdime,:) - 2*acvis*hess(idime,jdime,:)/3!-(2mu/3)*d_ik*d_dj*d^2_kj

            !Energy equation

            Kem(idime,:) = Kem(idime,:) + acvis*hess(idime,jdime,:)*gpadv(jdime)!mu*vel_j*d_kd*d^2_kj
            Kem(idime,:) = Kem(idime,:) - 2*acvis*gpadv(jdime)*hess(jdime,idime,:)/3!-(2mu/3)*vel_k*d_dj*d^2_kj
         end do
         
         Kem(idime,:) = Kem(idime,:) + acvis*gpadv(idime)*lapl(:)!mu*vel_d*d_kj*d^2_kj

      end do

      Kee = Kee + actco*lapl!tdiff*d_kj*d^2_kj

   end subroutine   

   !-------------------------------------------------------------------


end module
