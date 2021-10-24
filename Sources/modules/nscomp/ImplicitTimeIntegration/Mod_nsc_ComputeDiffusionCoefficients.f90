module Mod_nsc_ComputeDiffusionCoefficients
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersDiffusionCoefficients
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   real(rp)       :: kvis,tdiff
   
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
         
            call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateDiffusionCoefficients)
            call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,DiffusionCoefficientsToZero)
            call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeDiffusionCoefficients)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateDiffusionCoefficients)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateDiffusionCoefficients

      call a%Memor%alloc(e%ndime,e%mnode,Kmd,'Kmd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Kmm,'Kmm','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,Ked,'Ked','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Kem,'Kem','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,Kee,'Kee','nsc_elmope_im')

   end subroutine

   !Diffusion matrix coefficients
   subroutine DiffusionCoefficientsToZero

     !Stabilization Contribution                  
      Kmd = 0.0_rp !Kmd(i,p)
      Kmm = 0.0_rp !Kmm(i,d,p)
      Ked = 0.0_rp !Ked(p)
      Kem = 0.0_rp !Kem(d,p)
      Kee = 0.0_rp !Kee(p)

      kvis  = acvis*invgpd
      tdiff = actco*invgpd*invcvh

   end subroutine   

   subroutine ComputeDiffusionCoefficients
      implicit none    
      integer(ip)         :: idime,jdime

      !Momentum equation
      Kmd = Kmd + 2*kvis*invgpd*vgden*e%cartd!2*(nu/rho)*d_ij*vel_k*d_k(rho)*d_j
      Kmd = Kmd - kvis*invgpd*divmom*e%cartd!-(nu/rho)*d_ij*d_k(mom_k)*d_j
      Kmd = Kmd - kvis*invgpd*matmul(grmom,e%cartd)!-(nu/rho)*d_j(mom_i)*d_j
      Kmd = Kmd + 2*kvis*invgpd*matmul(transpose(grmom),e%cartd)/3!(2nu/3rho)*d_i(mom_j)*d_j

      do idime=1,e%ndime
         Kmd(idime,:) = Kmd(idime,:) + 2*kvis*invgpd*gpadv(idime)*matmul(grden,e%cartd)!2*(nu/rho)*vel_i*d_j(rho)*d_j
         Kmd(idime,:) = Kmd(idime,:) - 4*kvis*invgpd*grden(idime)*AGradV(:)/3!-(4nu/3rho)*vel_j*d_i(rho)*d_j

         Kmm(idime,idime,:) = Kmm(idime,idime,:) - kvis*invgpd*matmul(grden,e%cartd)!-(nu/rho)*d_id*d_j(rho)*d_j
         do jdime=1,e%ndime
            Kmm(idime,jdime,:) = Kmm(idime,jdime,:) - kvis*invgpd*e%cartd(idime,:)*grden(jdime)!-(nu/rho)*d_ij*d_d(rho)*d_j
            Kmm(idime,jdime,:) = Kmm(idime,jdime,:) + 2*kvis*invgpd*e%cartd(jdime,:)*grden(idime)/3!(2nu/3rho)*d_dj*d_i(rho)*d_j

         end do
         !Energy equation
         Kem(idime,:) = Kem(idime,:) + 2*(tdiff-kvis)*invgpd*gpadv(idime)*matmul(grden,e%cartd) !(2(tdiff-nu)/rho)*vel_d*d_j(rho)*d_j
         Kem(idime,:) = Kem(idime,:) - 2*kvis*invgpd*grden(idime)*AGradV(:)!-2(nu/rho)*vel_j*d_d(rho)*d_j

      end do
      Ked = Ked + 3*(kvis-tdiff)*invgpd*sqgpvn*matmul(grden,e%cartd)!(3(nu-tdiff)/rho)*|vel|^2*d_j(rho)*d_j
      Ked = Ked + 3*tdiff*sqinvgpd*gpade*matmul(grden,e%cartd)!(3(tdiff)/rho)*(ene/rho)*d_j(rho)*d_j
      Ked = Ked + 2*(tdiff-kvis)*invgpd*matmul(gpadv,matmul(grmom,e%cartd))!(2(tdiff-nu)/rho)*vel_k*d_j(mom_k)*d_j
      Ked = Ked + kvis*invgpd*AGradV*vgden!(mu/rho)*vel_j*vel_k*d_k(rho)*d_j
      Ked = Ked - kvis*invgpd*divmom*AGradV/3 !-(mu/3rho)*vel_j*d_k(mom_k)*d_j
      Ked = Ked - kvis*invgpd*matmul(vgmom,e%cartd)/3!(mu/3rho)*vel_k*d_k(mom_j)*d_j
      Ked = Ked - tdiff*invgpd*matmul(grene,e%cartd)!-(tdiff/rho)*d_j(ene)*d_j

      Kem = Kem + (kvis-tdiff)*invgpd*matmul(grmom,e%cartd)!((nu-tdiff)/rho)*d_j(mom_d)*d_j
      Kem = Kem + kvis*invgpd*matmul(transpose(grmom),e%cartd)!(nu/rho)*d_d(mom_j)*d_j
      Kem = Kem + 4*kvis*invgpd*vgden*e%cartd/3!(4nu/3rho)*vel_k*d_k(rho)*d_dj*d_j
      Kem = Kem - 2*kvis*invgpd*divmom*e%cartd/3!-(2nu/3rho)*d_k(mom_k)*d_dj*d_j

      Kee = Kee - tdiff*invgpd*matmul(grden,e%cartd)!(tdiff/rho)*d_j(rho)*d_j

   end subroutine   


   subroutine DeallocateDiffusionCoefficients

      call a%Memor%dealloc(e%ndime,e%mnode,Kmd,'Kmd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Kmm,'Kmm','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,Ked,'Ked','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Kem,'Kem','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,Kee,'Kee','nsc_elmope_im')

   end subroutine

end module
