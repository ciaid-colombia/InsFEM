module Mod_nsc_ComputeNonLinearDiffusionCoefficients
   use typre
   use Mod_nsc_BaseElmope
!   use Mod_NSCompressibleImplicitElement
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
         
            call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateHighOrderDerivatives)
            call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeHessian)
            call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeNonLinearDiffusionCoefficients)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateHighOrderDerivatives)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateHighOrderDerivatives
      
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,hess,'hess','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,lapl,'lapl','nsc_elmope_im')

   end subroutine

   !----------------------------------------------------------
   subroutine DeallocateHighOrderDerivatives
      
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,hess,'hess','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,lapl,'lapl','nsc_elmope_im')

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
      real(rp)       :: kvis,tdiff

      kvis  = acvis*invgpd
      tdiff = actco*invgpd*invcvh

     !Momentum equation
      do idime=1,e%ndime
         Kmd(idime,:) = Kmd(idime,:) - kvis*gpadv(idime)*lapl(:)!-(mu/rho)*vel_i*d_kj*d^2_kj
!
         Kmm(idime,idime,:) = Kmm(idime,idime,:) + kvis*lapl(:)!(mu/rho)*d_id*d_kj*d^2_kj
         do jdime=1,e%ndime
            Kmd(idime,:) = Kmd(idime,:) - kvis*gpadv(jdime)*hess(jdime,idime,:)!-(mu/rho)*vel_k*d_ij*d^2_kj
            Kmd(idime,:) = Kmd(idime,:) + 2*kvis*hess(idime,jdime,:)*gpadv(jdime)/3! +(2mu/3rho)*vel_j*d_ik*d^2_ik
            Kmm(idime,jdime,:) = Kmm(idime,jdime,:) + kvis*hess(jdime,idime,:)!(mu/rho)*d_ij*d_kd*d^2_kj
            Kmm(idime,jdime,:) = Kmm(idime,jdime,:) - 2*kvis*hess(idime,jdime,:)/3!-(2mu/3rho)*d_ik*d_dj*d^2_kj

            !Energy equation
            Ked(:) = Ked(:) - kvis*gpadv(idime)*gpadv(jdime)*hess(idime,jdime,:)/3 !-(mu/3rho)*vel_k*vel_j*d^2_kj

            Kem(idime,:) = Kem(idime,:) + kvis*hess(idime,jdime,:)*gpadv(jdime)!(mu/rho)*vel_j*d_kd*d^2_kj
            Kem(idime,:) = Kem(idime,:) - 2*kvis*gpadv(jdime)*hess(jdime,idime,:)/3!-(2mu/3rho)*vel_k*d_dj*d^2_kj
         end do
         
         Kem(idime,:) = Kem(idime,:) + kvis*gpadv(idime)*lapl(:)!(mu/rho)*vel_d*d_kj*d^2_kj
         Kem(idime,:) = Kem(idime,:) - tdiff*gpadv(idime)*lapl(:)!-tdiff*vel_d*d_kj*d^2_kj

      end do

      Ked = Ked - kvis*sqgpvn*lapl !-(mu/rho)*vel_i^2*d_kj*d^2_kj
      Ked = Ked - tdiff*gpade*invgpd*lapl !-tdiff*(ene/rho)*d_kj*d^2_kj
      Ked = Ked + tdiff*sqgpvn*lapl !tdiff*vel_i^2*d_kj*d^2_kj

      Kee = Kee + tdiff*lapl!tdiff*d_kj*d^2_kj

   end subroutine   

   !-------------------------------------------------------------------


end module
