  module Mod_nsc_pr_ComputeGpResidual
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersComputeGpResidual,gpred,gprem,gpree
   
   real(rp), allocatable :: gpred(:), gprem(:), gpree(:)
   real(rp), allocatable :: Ldd(:),   Ldm(:,:),   Lde(:)
   real(rp), allocatable :: Lmd(:,:), Lmm(:,:,:), Lme(:,:)
   real(rp), allocatable :: Led(:),   Lem(:,:),   Lee(:)
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SetPointersComputeGpResidual(itask) 
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
         
         call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocGpRes)
         call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocGpRes)
         call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,GpResToZero)

         !Convection
         if (a%kfl_advec == 1) then
            call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,ConvectionGpRes)
         endif
             
         !Diffusion Exists
         if (a%kfl_visco == 1) then
             call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,DiffusionGpRes)
         endif

         !Reaction Exists
         if (a%kfl_react == 1) then
            call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,ReactionGpRes)
         endif

         !Residual for OSS skipping FE components
         call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,ComputeGpRes)

         !Store complete FEM residual for shock capturing and output
         if ((a%kfl_shock == 1) .or. (a%npp_stepi(13) /= 0)) then
            call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,StoreGpRes)
         endif  

      endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select

   end subroutine   
   
   
   
  !-------------------------------------------------------------------
   !Residual computation
   subroutine AllocGpRes
      implicit none
      
      call a%Memor%alloc(1,gpred,'gpred','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,gprem,'gprem','nsc_pr_EndElmope')
      call a%Memor%alloc(1,gpree,'gpree','nsc_pr_EndElmope')
      call a%Memor%alloc(e%mnode,Ldd,'Ldd','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,Ldm,'Ldm','nsc_pr_EndElmope')
      call a%Memor%alloc(e%mnode,Lde,'Lde','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,Lmd,'Lmd','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Lmm,'Lmm','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,Lme,'Lme','nsc_pr_EndElmope')
      call a%Memor%alloc(e%mnode,Led,'Led','nsc_pr_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,Lem,'Lem','nsc_pr_EndElmope')
      call a%Memor%alloc(e%mnode,Lee,'Lee','nsc_pr_EndElmope')

   end subroutine
   
   
   subroutine DeAllocGpRes
      implicit none
      
      call a%Memor%dealloc(1,gpred,'gpred','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,gprem,'gprem','nsc_pr_EndElmope')
      call a%Memor%dealloc(1,gpree,'gpree','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%mnode,Ldd,'Ldd','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Ldm,'Ldm','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%mnode,Lde,'Lde','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Lmd,'Lmd','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Lmm,'Lmm','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Lme,'Lme','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%mnode,Led,'Led','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Lem,'Lem','nsc_pr_EndElmope')
      call a%Memor%dealloc(e%mnode,Lee,'Lee','nsc_pr_EndElmope')

   end subroutine
   
   !Residual matrices to zero
   subroutine GpResToZero

      gpred=0.0_rp
      gprem=0.0_rp
      gpree=0.0_rp
      Ldd = 0.0_rp !Ldd(p)
      Ldm = 0.0_rp !Ldm(d,p)
      Lde = 0.0_rp !Lde(p)
      Lmd = 0.0_rp !Lmd(i,p)
      Lmm = 0.0_rp !Lmm(i,d,p)
      Lme = 0.0_rp !Lme(i,p)
      Led = 0.0_rp !Led(p)
      Lem = 0.0_rp !Lem(d,p)
      Lee = 0.0_rp !Lee(p)
                        
   end subroutine   

   subroutine ConvectionGpRes
      implicit none
   
      !A_j*d_j U

      !Mass equation
      Ldd = Ldd + Add 
      Ldm = Ldm + Adm 
      Lde = Ldd + Ade 
      !Momentum equation
      Lmd = Lmd + Amd
      Lmm = Lmm + Amm
      Lme = Lme + Ame
      !Energy equation
      Led = Led + Aed  
      Lem = Lem + Aem  
      Lee = Lee + Aee  

   end subroutine
      
   subroutine DiffusionGpRes
      implicit none
   
      !-d_U(K_kj)·d_k(Uh)*d_j U
      !if NonLinear elements:  -K_kj*d^2_kj U 

      !Momentum equation
      Lmm = Lmm - Kmm
      !Energy equation
      Lem = Lem - Kem  
      Lee = Lee - Kee  

   end subroutine
   
   subroutine ReactionGpRes
      implicit none
   
      !-S U
      !Mass equation
      Ldd = Ldd - Sdd
      !Momentum equation
      Lmd = Lmd - Smd
      Lmm = Lmm - Smm
      !Energy equation
      Led = Led - Sed  
      Lem = Lem - Sem  
      Lee = Lee - See  

   end subroutine

   subroutine ComputeGpRes
      implicit none
      integer(ip) :: idime,jdime
      
     gpred(1) = gpred(1) + elexp(1) & 
      + gpden*acbeta*(eltemp(1)-gppre(1)*LHSdtinv)&
      - gpden*acalpha*(eltemt(1)-gptem(1)*LHSdtinv)
     gpred(1) = gpred(1) - dot_product(Ldd(:),elpre(:,1))
     gpred(1) = gpred(1) - dot_product(Lde(:),eltem(:,1))
     do idime = 1,e%ndime
        gpred(1) = gpred(1) - dot_product(Ldm(idime,:),elvel(idime,:,1))
        gprem(idime) = gprem(idime) + elexv(idime) &
         + gpden*acbeta*gpadv(idime)*(eltemp(1)-gppre(1)*LHSdtinv) &
         + gpden*(eltemv(idime)-gpvel(idime,1)*LHSdtinv) &
         - gpden*acalpha*gpadv(idime)*(eltemt(1)-gptem(1)*LHSdtinv)
        gprem(idime) = gprem(idime) - &
                       dot_product(Lmd(idime,:),elpre(:,1)) -&
                       dot_product(Lme(idime,:),eltem(:,1))
        do jdime = 1,e%ndime
           gprem(idime) = gprem(idime) - &
                          dot_product(Lmm(idime,jdime,:),elvel(jdime,:,1)) 
        end do
        gpree(1) = gpree(1) - dot_product(Lem(idime,:),elvel(idime,:,1))
     end do
     gpree(1) = gpree(1) - dot_product(Led(:),elpre(:,1))
     gpree(1) = gpree(1) - dot_product(Lee(:),eltem(:,1))
     gpree(1) = gpree(1) + elext(1) &
      + (aux*acbeta-acalpha*gpadt)*(eltemp(1)-gppre(1)*LHSdtinv) &
      + gpden*dot_product(gpadv,eltemv-gpvel(:,1)*LHSdtinv) &
      + (gpden*accph-aux*acalpha)*(eltemt(1)-gptem(1)*LHSdtinv)
                
      
   end subroutine

   subroutine StoreGpRes
      implicit none
      
      a%coResGP(ielem)%a(e%igaus) = gpred(1) 

      a%moResGP(ielem)%a(:,e%igaus) = gprem(:)

      a%enResGP(ielem)%a(e%igaus) = gpree(1) 

   end subroutine
   
end module


