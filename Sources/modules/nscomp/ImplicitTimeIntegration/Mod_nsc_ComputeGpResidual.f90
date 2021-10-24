  module Mod_nsc_ComputeGpResidual
   use typre
   use Mod_nsc_BaseElmope
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
         
         call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocGpRes)
         call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocGpRes)
         call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,GpResToZero)

         !Convection
         if (a%kfl_advec == 1) then
            call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,ConvectionGpRes)
         endif
             
         !Diffusion Exists
         if (a%kfl_visco == 1) then
             call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,DiffusionGpRes)
         endif

         !Reaction Exists
         if (a%kfl_react == 1) then
            call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,ReactionGpRes)
         endif

         !Residual for OSS skipping FE components
         call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,ComputeGpRes)

         !Store complete FEM residual for shock capturing and output
         if ((a%kfl_shock == 1) .or. (a%npp_stepi(13) /= 0)) then
            call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,StoreGpRes)
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
      
      call a%Memor%alloc(1,gpred,'gpred','nsc_EndElmope_im')
      call a%Memor%alloc(e%ndime,gprem,'gprem','nsc_EndElmope_im')
      call a%Memor%alloc(1,gpree,'gpree','nsc_EndElmope_im')
      call a%Memor%alloc(e%mnode,Ldd,'Ldd','nsc_EndElmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Ldm,'Ldm','nsc_EndElmope_im')
      call a%Memor%alloc(e%mnode,Lde,'Lde','nsc_EndElmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Lmd,'Lmd','nsc_EndElmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Lmm,'Lmm','nsc_EndElmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Lme,'Lme','nsc_EndElmope_im')
      call a%Memor%alloc(e%mnode,Led,'Led','nsc_EndElmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Lem,'Lem','nsc_EndElmope_im')
      call a%Memor%alloc(e%mnode,Lee,'Lee','nsc_EndElmope_im')

   end subroutine
   
   
   subroutine DeAllocGpRes
      implicit none
      
      call a%Memor%dealloc(1,gpred,'gpred','nsc_EndElmope_im')
      call a%Memor%dealloc(e%ndime,gprem,'gprem','nsc_EndElmope_im')
      call a%Memor%dealloc(1,gpree,'gpree','nsc_EndElmope_im')
      call a%Memor%dealloc(e%mnode,Ldd,'Ldd','nsc_EndElmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Ldm,'Ldm','nsc_EndElmope_im')
      call a%Memor%dealloc(e%mnode,Lde,'Lde','nsc_EndElmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Lmd,'Lmd','nsc_EndElmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Lmm,'Lmm','nsc_EndElmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Lme,'Lme','nsc_EndElmope_im')
      call a%Memor%dealloc(e%mnode,Led,'Led','nsc_EndElmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Lem,'Lem','nsc_EndElmope_im')
      call a%Memor%dealloc(e%mnode,Lee,'Lee','nsc_EndElmope_im')

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
      Lmd = Lmd - Kmd
      Lmm = Lmm - Kmm
      !Energy equation
      Led = Led - Ked  
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
      
     gpred(1) = gpred(1) + dot_product(Ldd(:),elden(:,1))
     do idime = 1,e%ndime
        gpred(1) = gpred(1) + dot_product(Ldm(idime,:),elmom(idime,:,1))
        gprem(idime) = gprem(idime) + &
                       dot_product(Lmd(idime,:),elden(:,1)) +&
                       dot_product(Lme(idime,:),elene(:,1))
        do jdime = 1,e%ndime
           gprem(idime) = gprem(idime) + &
                          dot_product(Lmm(idime,jdime,:),elmom(jdime,:,1)) 
        end do
        gpree(1) = gpree(1) + dot_product(Lem(idime,:),elmom(idime,:,1))
     end do
     gpree(1) = gpree(1) + dot_product(Led(:),elden(:,1))
     gpree(1) = gpree(1) + dot_product(Lee(:),elene(:,1))
                
      
   end subroutine

   subroutine StoreGpRes
      implicit none
      
      a%coResGP(ielem)%a(e%igaus) = - gpred(1) - LHSdtinv*gpden(1) + elexd(1) + eltemd(1) 
      a%moResGP(ielem)%a(:,e%igaus) = - gprem(:) - LHSdtinv*gpmom(:,1) + elexm(:) + eltemm(:)
      a%enResGP(ielem)%a(e%igaus) = - gpree(1) - LHSdtinv*gpene(1) + elexe(1) + elteme(1)

   end subroutine
   
end module


