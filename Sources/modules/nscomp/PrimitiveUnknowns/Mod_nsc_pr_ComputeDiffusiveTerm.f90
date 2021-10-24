module Mod_nsc_pr_ComputeShockCapturing
   use typre
   use Mod_NSCompressibleElement
!   use Mod_NSCompressiblePrimitiveElement
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_InterpolateGradientProjection
   implicit none
   private
   public SetPointersShockCapturing,DCvis,DCtdf,STvis,STtdf
   
   real(rp)       :: DCvis,DCtdf,STvis,STtdf
   real(rp), allocatable :: ortht(:), orthv(:,:)
   
   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersShockCapturing(itask)
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
         
            !Shock Capturing
            if (a%kfl_shock /= 0 ) then
                 if (a%kfl_shock == 1) then
                    call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ResidualDiffusionShockCapturing)
                 else if (a%kfl_shock == 2) then
                    call SetPointersInterpolateGradientProjection(1)
                    call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateGradProjDiffusion)
                    call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateGradProjDiffusion)
                    call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,GradientProjectionDiffusionShockCapturing)
                 else if (a%kfl_shock > 2) then
                    call runend('Nsc_EndElmope: Other shock capturing methods not ready')
                 endif

            end if

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   !Residual Diffusion Shock Capturing

   subroutine ResidualDiffusionShockCapturing

      call nsc_ArtificialVis(e,a%shock,chale,a%moResGP(ielem)%a(:,e%igaus),grvel,DCvis)
      call nsc_ArtificialCnd(e,a%shock,chale,a%enResGP(ielem)%a(e%igaus),grtem,DCtdf)

   end subroutine

   !-----------------------------------------------------------------------
   !Gradient Orthogonal Projection Shock Capturing
   subroutine AllocateGradProjDiffusion
      
      call a%Memor%alloc(e%ndime,e%ndime,orthv,'orthv','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,ortht,'ortht','nsc_pr_elmope')

   end subroutine

   subroutine DeallocateGradProjDiffusion
      
      call a%Memor%dealloc(e%ndime,e%ndime,orthv,'orthv','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,ortht,'ortht','nsc_pr_elmope')

   end subroutine

   subroutine GradientProjectionDiffusionShockCapturing

      orthv = grvel - gprjv
      ortht = grtem - gprjt
      
      call nsc_GradOrthVis(e,a%shock,chale,gpvno,orthv,grvel,DCvis)
      call nsc_GradOrthCnd(e,a%shock,chale,gpvno,ortht,grtem,DCtdf)

   end subroutine

end module

module Mod_nsc_pr_ComputeDiffusiveTerm
   use typre
   use Mod_nsc_pr_BaseElmope
!   use Mod_NSCompressiblePrimitiveElement
   use Mod_nsc_pr_ComputeShockCapturing
   implicit none
   private
   public SetPointersDiffusiveTerm
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   real(rp), allocatable :: untmv(:)
   real(rp), allocatable :: Anid(:,:), secstr(:,:), secort(:,:)
   real(rp), allocatable :: Aniv(:,:), forstr(:,:), forort(:,:)
   real(rp)       :: Kvis,Ktdf
   real(rp)       :: DCstv,DCstd,dummr
   real(rp)       :: dampcd
   
contains

   subroutine SetPointersDiffusiveTerm(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateDiffusiveTerm)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,DiffusiveTermToZero)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateDiffusiveTerm)

            !Shock Capturing
            if (a%kfl_shock /= 0 ) then
                call SetPointersShockCapturing(1)
                if (a%kfl_sctyp == 0) then
                   call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ComputeIsotropicAddedDiffusion)
                elseif (a%kfl_sctyp > 0) then
                    call runend('Nsc_pr_Elmope: Anisotropic shock capturing method not ready')
                !elseif (a%kfl_sctyp == 1) then
!                   call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateAnisotropic)
!                   call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,AnisotropicToZero)
!                   call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateAnisotropic)
!                   call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeAnisotropicAddedDiffusion)
!                   call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeAnisotropicTensors)
!                   if(e%ndime.eq.2) then
!                      call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeAnisotropicDiffusiveTerm_Bidim)
!                   elseif(e%ndime.eq.3) then
!                      call ConcatenateProcedures(ProcPointer_nsc_ComputeDiffusionCoefficients,ComputeAnisotropicDiffusiveTerm_Tridim)
!                   endif
                endif
            endif

            !Radial Damping Exists (with added diffusion)
            if (a%rdamp > 0) then
!               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ComputeRadialDampingDiffusion)
            endif

            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeDiffusionCoefficients,ComputeIsotropicDiffusiveTerm)

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine

   !-----------------------------------------------------------------------
   subroutine AllocateDiffusiveTerm

      call a%Memor%alloc(e%ndime,e%mnode,KGdd,'KGdd','nsc_pr_elmope') !Radial Damping Diffusion
      call a%Memor%alloc(e%ndime,e%ndime,e%ndime,e%mnode,KGmm,'KGmm','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,KGem,'KGem','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,KGee,'KGee','nsc_pr_elmope')

   end subroutine

   !Diffusive Term coefficients
   subroutine DiffusiveTermToZero

      Kvis = acvis
      Ktdf = actco
      
     !Galerkin Contribution                  
      KGdd = 0.0_rp !KGdd(k,p) !Radial Damping Diffusion
      KGmm = 0.0_rp !KGmm(k,i,d,p)
      KGem = 0.0_rp !KGem(k,d,p)
      KGee = 0.0_rp !KGee(k,p)

   end subroutine   

   subroutine DeAllocateDiffusiveTerm

      call a%Memor%dealloc(e%ndime,e%mnode,KGdd,'KGdd','nsc_pr_elmope') !Radial Damping Diffusion
      call a%Memor%dealloc(e%ndime,e%ndime,e%ndime,e%mnode,KGmm,'KGmm','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,KGem,'KGem','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,KGee,'KGee','nsc_pr_elmope')

   end subroutine

   subroutine ComputeIsotropicAddedDiffusion

      Kvis  = Kvis + gpden*DCvis 
      Ktdf  = Ktdf + gpden*accvh*DCtdf

   end subroutine

   subroutine ComputeIsotropicDiffusiveTerm
      implicit none    
      integer(ip)         :: idime,jdime

      !Momentum equation
      do idime=1,e%ndime
         KGmm(idime,idime,:,:) = KGmm(idime,idime,:,:) - 2*Kvis*e%cartd(:,:)/3!-(2mu/3)*d_ik*d_dj*d_j
         
         do jdime=1,e%ndime
            KGmm(idime,jdime,idime,:) = KGmm(idime,jdime,idime,:) + Kvis*e%cartd(jdime,:)!mu*d_ij*d_kd*d_j
            KGmm(idime,jdime,jdime,:) = KGmm(idime,jdime,jdime,:) + Kvis*e%cartd(idime,:)!mu*d_id*d_kj*d_j

            !Energy equation
            KGem(idime,jdime,:) = KGem(idime,jdime,:) + Kvis*gpadv(jdime)*e%cartd(idime,:)!mu*vel_d*d_kj*d_j
         end do
         KGem(idime,idime,:) = KGem(idime,idime,:) + Kvis*AGradV(:)!nu*vel_j*d_kd*d_j
         KGem(idime,:,:) = KGem(idime,:,:) - 2*Kvis*gpadv(idime)*e%cartd(:,:)/3!-(2nu/3)*vel_k*d_dj*d_j
      
      end do

      KGee = KGee + Ktdf*e%cartd !tdiff*d_kj*d_j


   end subroutine   

   subroutine ComputeRadialDampingDiffusion
     use Mod_Mesh    
     implicit none      
     real(rp)    :: gpcod(e%ndime)

     !Interpolate
     call e%interpg(e%ndime,e%elcod,gpcod)

     !Compute coefficient of damping
     call damping%nsc_ComputeRadialDamping(e%ndime,gpcod,a)
     call damping%nsc_GetDampingCoefficient(dampcd)
     
     !Mass equation
      KGdd  = KGdd + gpden*dampcd*acvis*e%cartd/gpadp !tdiff*d_kj*d_j At the order of acvis
      Kvis  = Kvis + dampcd*acvis/gpadp  !At the order of pressure
      Ktdf  = Ktdf + dampcd*actco/gpadp  !At the order of pressure
 
   end subroutine   

!   subroutine AllocateAnisotropic
!      implicit none    
!      integer(ip)         :: ntens
!     
!      ntens = (e%ndime-1)*(e%ndime-1)+2  
!      call a%Memor%alloc(e%ndime,e%ndime,secstr,'secstr','nsc_pr_elmope')
!      call a%Memor%alloc(e%ndime,e%ndime,secort,'secort','nsc_pr_elmope')
!      call a%Memor%alloc(e%ndime,e%ndime,Anid,'Anid','nsc_pr_elmope')
!      call a%Memor%alloc(ntens,ntens,forstr,'forstr','nsc_pr_elmope')
!      call a%Memor%alloc(ntens,ntens,forort,'forort','nsc_pr_elmope')
!      call a%Memor%alloc(ntens,ntens,Aniv,'Aniv','nsc_pr_elmope')
!      call a%Memor%alloc(e%ndime,untmv,'untmv','nsc_pr_elmope')
!
!   end subroutine
!
!   subroutine AnisotropicToZero
!
!      untmv = gpadm
!      call vecuni(e%ndime,untmv,dummr)
!      forstr = 0.0_rp
!      forort = 0.0_rp
!      Aniv = 0.0_rp
!      secstr = 0.0_rp
!      secort = 0.0_rp
!      Anid = 0.0_rp
!
!   end subroutine
!
!   subroutine DeallocateAnisotropic
!      implicit none    
!      integer(ip)         :: ntens
!     
!      ntens = (e%ndime-1)*(e%ndime-1)+2  
!      call a%Memor%dealloc(e%ndime,e%ndime,secstr,'secstr','nsc_pr_elmope')
!      call a%Memor%dealloc(e%ndime,e%ndime,secort,'secort','nsc_pr_elmope')
!      call a%Memor%dealloc(e%ndime,e%ndime,Anid,'Anid','nsc_pr_elmope')
!      call a%Memor%dealloc(ntens,ntens,forstr,'forstr','nsc_pr_elmope')
!      call a%Memor%dealloc(ntens,ntens,forort,'forort','nsc_pr_elmope')
!      call a%Memor%dealloc(ntens,ntens,Aniv,'Aniv','nsc_pr_elmope')
!      call a%Memor%dealloc(e%ndime,untmv,'untmv','nsc_pr_elmope')
!
!   end subroutine
!
!   subroutine ComputeAnisotropicAddedDiffusion
!
!      DCstv = max(0.0_rp,DCvis-sqgpvn*timom(2))
!      DCstd = max(0.0_rp,DCtdf-sqgpvn*timom(3))
!
!   end subroutine
!
!   subroutine ComputeAnisotropicTensors
!
!      call nsc_FourthStreamOrthogonalTensor(e%ndime,untmv,forstr,forort)
!      Aniv = DCstv*forstr+DCvis*forort
!
!      call nsc_StreamlineTensor(e%ndime,untmv,secstr)
!      call nsc_OrthogonalTensor(e%ndime,secstr,secort)
!      Anid = DCstd*secstr+DCtdf*secort
!
!   end subroutine
!
!   subroutine ComputeAnisotropicDiffusiveTerm_Bidim
!      implicit none    
!      integer(ip)         :: idime
!
!      !Momentum equation
!      !- nu*vel_k*d_ij*d_j - nu*vel_i*d_kj*d_j + (2nu/3)*vel_j*d_ij*d_j
!      KGmd(1,1,:) = KGmd(1,1,:) + Aniv(1,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3)&
!                                + Aniv(1,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) 
!      KGmd(2,2,:) = KGmd(2,2,:) + Aniv(2,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3)&
!                                + Aniv(2,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) 
!      KGmd(1,2,:) = KGmd(1,2,:) + Aniv(3,3)*(-gpadv(1)*e%cartd(2,:)-gpadv(2)*e%cartd(1,:)) 
!      KGmd(2,1,:) = KGmd(2,1,:) + Aniv(3,3)*(-gpadv(1)*e%cartd(2,:)-gpadv(2)*e%cartd(1,:)) 
!
!      ! nu*d_ij*d_kd*d_j + nu*d_id*d_kj*d_j - nu*(2/3)*d_ik*d_dj*d_j
!      KGmm(1,1,1,:) = KGmm(1,1,1,:) + Aniv(1,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + Aniv(1,2)*(-2*e%cartd(1,:)/3)+Aniv(1,3)*e%cartd(2,:) 
!      KGmm(1,1,2,:) = KGmm(1,1,2,:) + Aniv(1,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + Aniv(1,1)*(-2*e%cartd(2,:)/3)+Aniv(1,3)*e%cartd(1,:) 
!      KGmm(2,2,1,:) = KGmm(2,2,1,:) + Aniv(2,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + Aniv(2,2)*(-2*e%cartd(1,:)/3)+Aniv(2,3)*e%cartd(2,:) 
!      KGmm(2,2,2,:) = KGmm(2,2,2,:) + Aniv(2,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + Aniv(2,1)*(-2*e%cartd(2,:)/3)+Aniv(2,3)*e%cartd(1,:) 
!      KGmm(1,2,1,:) = KGmm(1,2,1,:) + Aniv(3,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + Aniv(3,2)*(-2*e%cartd(1,:)/3)+Aniv(3,3)*e%cartd(2,:) 
!      KGmm(1,2,2,:) = KGmm(1,2,2,:) + Aniv(3,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + Aniv(3,1)*(-2*e%cartd(2,:)/3)+Aniv(3,3)*e%cartd(1,:) 
!      KGmm(2,1,1,:) = KGmm(2,1,1,:) + Aniv(3,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + Aniv(3,2)*(-2*e%cartd(1,:)/3)+Aniv(3,3)*e%cartd(2,:) 
!      KGmm(2,1,2,:) = KGmm(2,1,2,:) + Aniv(3,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + Aniv(3,1)*(-2*e%cartd(2,:)/3)+Aniv(3,3)*e%cartd(1,:) 
!
!      !Energy equation
!      !- vel_i*nu*vel_k*d_ij*d_j - vel_i*nu*vel_i*d_kj*d_j + vel_i*(2nu/3)*vel_j*d_ij*d_j
!      KGed(1,:) = KGed(1,:) + gpadv(1)*Aniv(1,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3)&
!                            + gpadv(1)*Aniv(1,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3)&
!                            + gpadv(2)*Aniv(3,3)*(-gpadv(1)*e%cartd(2,:)-gpadv(2)*e%cartd(1,:)) 
!      KGed(2,:) = KGed(2,:) + gpadv(1)*Aniv(3,3)*(-gpadv(1)*e%cartd(2,:)-gpadv(2)*e%cartd(1,:))& 
!                            + gpadv(2)*Aniv(2,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3)&
!                            + gpadv(2)*Aniv(2,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) 
!      KGed = KGed + sqgpvn*matmul(Anid,e%cartd)!tdiff*|vel|^2*d_kj*d_j
!      KGed = KGed - gpade*invgpd*matmul(Anid,e%cartd) !-tdiff*(ene/rho)*d_kj*d_j
!
!      ! vel_i*nu*d_ij*d_kd*d_j + vel_i*nu*d_id*d_kj*d_j - vel_i*nu*(2/3)*d_ik*d_dj*d_j
!      KGem(1,1,:) = KGem(1,1,:) + gpadv(1)*Aniv(1,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(1)*Aniv(1,2)*(-2*e%cartd(1,:)/3)+Aniv(1,3)*e%cartd(2,:)& 
!                                    + gpadv(2)*Aniv(3,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(2)*Aniv(3,2)*(-2*e%cartd(1,:)/3)+Aniv(3,3)*e%cartd(2,:) 
!      KGem(1,2,:) = KGem(1,2,:) + gpadv(1)*Aniv(1,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(1)*Aniv(1,1)*(-2*e%cartd(2,:)/3)+Aniv(1,3)*e%cartd(1,:)& 
!                                    + gpadv(2)*Aniv(3,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(2)*Aniv(3,1)*(-2*e%cartd(2,:)/3)+Aniv(3,3)*e%cartd(1,:) 
!      KGem(2,1,:) = KGem(2,1,:) + gpadv(1)*Aniv(3,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(1)*Aniv(3,2)*(-2*e%cartd(1,:)/3)+Aniv(3,3)*e%cartd(2,:)& 
!                                    + gpadv(2)*Aniv(2,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(2)*Aniv(2,2)*(-2*e%cartd(1,:)/3)+Aniv(2,3)*e%cartd(2,:) 
!      KGem(2,2,:) = KGem(2,2,:) + gpadv(1)*Aniv(3,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(1)*Aniv(3,1)*(-2*e%cartd(2,:)/3)+Aniv(3,3)*e%cartd(1,:)& 
!                                    + gpadv(2)*Aniv(2,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(2)*Aniv(2,1)*(-2*e%cartd(2,:)/3)+Aniv(2,3)*e%cartd(1,:) 
!      do idime=1,e%ndime
!      KGem(:,idime,:) = KGem(:,idime,:) - gpadv(idime)*matmul(Anid,e%cartd)!-tdiff*vel_d*d_kj*d_j
!      end do
!      
!      KGee = KGee + matmul(Anid,e%cartd) !tdiff*d_kj*d_j
!
!   end subroutine   
!
!   subroutine ComputeAnisotropicDiffusiveTerm_Tridim
!      implicit none    
!      integer(ip)         :: idime
!
!      !Momentum equation
!      !- nu*vel_k*d_ij*d_j - nu*vel_i*d_kj*d_j + (2nu/3)*vel_j*d_ij*d_j
!      KGmd(1,1,:) = KGmd(1,1,:) + Aniv(1,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3) &
!                                + Aniv(1,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) &
!                                + Aniv(1,3)*(-2*gpadv(3)*e%cartd(3,:)+2*AGradV(:)/3) 
!      KGmd(2,2,:) = KGmd(2,2,:) + Aniv(2,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3) &
!                                + Aniv(2,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) &
!                                + Aniv(2,3)*(-2*gpadv(3)*e%cartd(3,:)+2*AGradV(:)/3) 
!      KGmd(3,3,:) = KGmd(3,3,:) + Aniv(3,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3) &
!                                + Aniv(3,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) &
!                                + Aniv(3,3)*(-2*gpadv(3)*e%cartd(3,:)+2*AGradV(:)/3) 
!      KGmd(1,2,:) = KGmd(1,2,:) - Aniv(6,6)*(gpadv(1)*e%cartd(2,:)+gpadv(2)*e%cartd(1,:)) 
!      KGmd(2,1,:) = KGmd(2,1,:) - Aniv(6,6)*(gpadv(1)*e%cartd(2,:)+gpadv(2)*e%cartd(1,:)) 
!      KGmd(1,3,:) = KGmd(1,3,:) - Aniv(5,5)*(gpadv(1)*e%cartd(3,:)+gpadv(3)*e%cartd(1,:)) 
!      KGmd(3,1,:) = KGmd(3,1,:) - Aniv(5,5)*(gpadv(1)*e%cartd(3,:)+gpadv(3)*e%cartd(1,:)) 
!      KGmd(2,3,:) = KGmd(2,3,:) - Aniv(4,4)*(gpadv(2)*e%cartd(3,:)+gpadv(3)*e%cartd(2,:)) 
!      KGmd(3,2,:) = KGmd(3,2,:) - Aniv(4,4)*(gpadv(2)*e%cartd(3,:)+gpadv(3)*e%cartd(2,:)) 
!
!      ! nu*d_ij*d_kd*d_j + nu*d_id*d_kj*d_j - nu*(2/3)*d_ik*d_dj*d_j
!      KGmm(1,1,1,:) = KGmm(1,1,1,:) + Aniv(1,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + (Aniv(1,2) + Aniv(1,3))*(-2*e%cartd(1,:)/3)
!      KGmm(1,1,2,:) = KGmm(1,1,2,:) + Aniv(1,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + (Aniv(1,1) + Aniv(1,3))*(-2*e%cartd(2,:)/3)
!      KGmm(1,1,3,:) = KGmm(1,1,3,:) + Aniv(1,3)*(2*e%cartd(3,:) - 2*e%cartd(3,:)/3)&
!                                    + (Aniv(1,1) + Aniv(1,2))*(-2*e%cartd(3,:)/3)
!      KGmm(2,2,1,:) = KGmm(2,2,1,:) + Aniv(2,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + (Aniv(2,2) + Aniv(2,3))*(-2*e%cartd(1,:)/3)
!      KGmm(2,2,2,:) = KGmm(2,2,2,:) + Aniv(2,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + (Aniv(2,1) + Aniv(2,3))*(-2*e%cartd(2,:)/3)
!      KGmm(2,2,3,:) = KGmm(2,2,3,:) + Aniv(2,3)*(2*e%cartd(3,:) - 2*e%cartd(3,:)/3)&
!                                    + (Aniv(2,1) + Aniv(2,2))*(-2*e%cartd(3,:)/3)
!      KGmm(3,3,1,:) = KGmm(3,3,1,:) + Aniv(3,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + (Aniv(3,2) + Aniv(3,3))*(-2*e%cartd(1,:)/3)
!      KGmm(3,3,2,:) = KGmm(3,3,2,:) + Aniv(3,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + (Aniv(3,1) + Aniv(3,3))*(-2*e%cartd(2,:)/3)
!      KGmm(3,3,3,:) = KGmm(3,3,3,:) + Aniv(3,3)*(2*e%cartd(3,:) - 2*e%cartd(3,:)/3)&
!                                    + (Aniv(3,1) + Aniv(3,2))*(-2*e%cartd(3,:)/3)
!      KGmm(1,2,1,:) = KGmm(1,2,1,:) + Aniv(6,6)*e%cartd(2,:) 
!      KGmm(1,2,2,:) = KGmm(1,2,2,:) + Aniv(6,6)*e%cartd(1,:) 
!      KGmm(2,1,1,:) = KGmm(2,1,1,:) + Aniv(6,6)*e%cartd(2,:) 
!      KGmm(2,1,2,:) = KGmm(2,1,2,:) + Aniv(6,6)*e%cartd(1,:) 
!      KGmm(1,3,1,:) = KGmm(1,3,1,:) + Aniv(5,5)*e%cartd(3,:) 
!      KGmm(1,3,3,:) = KGmm(1,3,3,:) + Aniv(5,5)*e%cartd(1,:) 
!      KGmm(3,1,1,:) = KGmm(3,1,1,:) + Aniv(5,5)*e%cartd(3,:) 
!      KGmm(3,1,3,:) = KGmm(3,1,3,:) + Aniv(5,5)*e%cartd(1,:) 
!      KGmm(2,3,2,:) = KGmm(2,3,2,:) + Aniv(4,4)*e%cartd(3,:) 
!      KGmm(2,3,3,:) = KGmm(2,3,3,:) + Aniv(4,4)*e%cartd(2,:) 
!      KGmm(3,2,2,:) = KGmm(3,2,2,:) + Aniv(4,4)*e%cartd(3,:) 
!      KGmm(3,2,3,:) = KGmm(3,2,3,:) + Aniv(4,4)*e%cartd(2,:) 
!
!      !Energy equation
!      !- vel_i*nu*vel_k*d_ij*d_j - vel_i*nu*vel_i*d_kj*d_j + vel_i*(2nu/3)*vel_j*d_ij*d_j
!      KGed(1,:) = KGed(1,:) + gpadv(1)*Aniv(1,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3) &
!                                + gpadv(1)*Aniv(1,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) &
!                                + gpadv(1)*Aniv(1,3)*(-2*gpadv(3)*e%cartd(3,:)+2*AGradV(:)/3) &
!                                - gpadv(2)*Aniv(6,6)*(gpadv(1)*e%cartd(2,:)+gpadv(2)*e%cartd(1,:)) &
!                                - gpadv(3)*Aniv(5,5)*(gpadv(1)*e%cartd(3,:)+gpadv(3)*e%cartd(1,:)) 
!      KGed(2,:) = KGed(2,:) - gpadv(1)*Aniv(6,6)*(gpadv(1)*e%cartd(2,:)+gpadv(2)*e%cartd(1,:)) &
!                                + gpadv(2)*Aniv(2,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3) &
!                                + gpadv(2)*Aniv(2,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) &
!                                + gpadv(2)*Aniv(2,3)*(-2*gpadv(3)*e%cartd(3,:)+2*AGradV(:)/3) & 
!                                - gpadv(3)*Aniv(4,4)*(gpadv(2)*e%cartd(3,:)+gpadv(3)*e%cartd(2,:)) 
!      KGed(3,:) = KGed(3,:) - gpadv(1)*Aniv(5,5)*(gpadv(1)*e%cartd(3,:)+gpadv(3)*e%cartd(1,:)) &
!                                - gpadv(2)*Aniv(4,4)*(gpadv(2)*e%cartd(3,:)+gpadv(3)*e%cartd(2,:)) &
!                                + gpadv(3)*Aniv(3,1)*(-2*gpadv(1)*e%cartd(1,:)+2*AGradV(:)/3) &
!                                + gpadv(3)*Aniv(3,2)*(-2*gpadv(2)*e%cartd(2,:)+2*AGradV(:)/3) &
!                                + gpadv(3)*Aniv(3,3)*(-2*gpadv(3)*e%cartd(3,:)+2*AGradV(:)/3) 
!     KGed = KGed + Ktdf*sqgpvn*matmul(Anid,e%cartd)!tdiff*|vel|^2*d_kj*d_j
!     KGed = KGed - gpade*invgpd*matmul(Anid,e%cartd) !-tdiff*(ene/rho)*d_kj*d_j
!
!     ! vel_i*nu*d_ij*d_kd*d_j + vel_i*nu*d_id*d_kj*d_j - vel_i*nu*(2/3)*d_ik*d_dj*d_j
!      KGem(1,1,:) = KGem(1,1,:) + gpadv(1)*Aniv(1,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(1)*(Aniv(1,2) + Aniv(1,3))*(-2*e%cartd(1,:)/3)&
!                                    + gpadv(2)*Aniv(6,6)*e%cartd(2,:) &
!                                    + gpadv(3)*Aniv(5,5)*e%cartd(3,:) 
!      KGem(1,2,:) = KGem(1,2,:) + gpadv(1)*Aniv(1,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(1)*(Aniv(1,1) + Aniv(1,3))*(-2*e%cartd(2,:)/3)&
!                                    + gpadv(2)*Aniv(6,6)*e%cartd(1,:) 
!      KGem(1,3,:) = KGem(1,3,:) + gpadv(1)*Aniv(1,3)*(2*e%cartd(3,:) - 2*e%cartd(3,:)/3)&
!                                    + gpadv(1)*(Aniv(1,1) + Aniv(1,2))*(-2*e%cartd(3,:)/3)&
!                                    + gpadv(3)*Aniv(5,5)*e%cartd(1,:) 
!      KGem(2,1,:) = KGem(2,1,:) + gpadv(1)*Aniv(6,6)*e%cartd(2,:) &
!                                    + gpadv(2)*Aniv(2,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(2)*(Aniv(2,2) + Aniv(2,3))*(-2*e%cartd(1,:)/3)
!      KGem(2,2,:) = KGem(2,2,:) + gpadv(1)*Aniv(6,6)*e%cartd(1,:) &
!                                    + gpadv(2)*Aniv(2,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(2)*(Aniv(2,1) + Aniv(2,3))*(-2*e%cartd(2,:)/3)&
!                                    + gpadv(3)*Aniv(4,4)*e%cartd(3,:) 
!      KGem(2,3,:) = KGem(2,3,:) + gpadv(2)*Aniv(2,3)*(2*e%cartd(3,:) - 2*e%cartd(3,:)/3)&
!                                    + gpadv(2)*(Aniv(2,1) + Aniv(2,2))*(-2*e%cartd(3,:)/3)&
!                                    + gpadv(3)*Aniv(4,4)*e%cartd(2,:) 
!      KGem(3,1,:) = KGem(3,1,:) + gpadv(1)*Aniv(5,5)*e%cartd(3,:) &
!                                    + gpadv(3)*Aniv(3,1)*(2*e%cartd(1,:) - 2*e%cartd(1,:)/3)&  
!                                    + gpadv(3)*(Aniv(3,2) + Aniv(3,3))*(-2*e%cartd(1,:)/3)
!      KGem(3,2,:) = KGem(3,2,:) + gpadv(2)*Aniv(4,4)*e%cartd(3,:)& 
!                                    + gpadv(3)*Aniv(3,2)*(2*e%cartd(2,:) - 2*e%cartd(2,:)/3)&
!                                    + gpadv(3)*(Aniv(3,1) + Aniv(3,3))*(-2*e%cartd(2,:)/3)
!      KGem(3,3,:) = KGem(3,3,:) + gpadv(1)*Aniv(5,5)*e%cartd(1,:) &
!                                    + gpadv(2)*Aniv(4,4)*e%cartd(2,:)& 
!                                    + gpadv(3)*Aniv(3,3)*(2*e%cartd(3,:) - 2*e%cartd(3,:)/3)&
!                                    + gpadv(3)*(Aniv(3,1) + Aniv(3,2))*(-2*e%cartd(3,:)/3)
!     do idime=1,e%ndime
!     KGem(:,idime,:) = KGem(:,idime,:) - gpadv(idime)*matmul(Anid,e%cartd)!-tdiff*vel_d*d_kj*d_j
!     end do
!
!     KGee = KGee + matmul(Anid,e%cartd) !tdiff*d_kj*d_j
!
!   end subroutine   


end module

