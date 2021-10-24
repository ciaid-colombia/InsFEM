module Mod_nsc_ComputeTestf
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersComputeTestf

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeTestf(itask)
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
            
            !the elemental stabilization matrix takes the form
            !************************************************
            !q· |L*dd.Td  L*dm.Tm  L*de.Te| |Ldd  Ldm  Lde| |rho|
            !n· |L*md.Td  L*mm.Tm  L*me.Te|.|Lmd  Lmm  Lme|.|mom| 
            !g· |L*ed.Td  L*em.Tm  L*ee.Te| |Led  Lem  Lee| |ene| 
            !************************************************

            !Allocations
            call ConcatenateProcedures(ProcHook_nsc_Initializations, AllocTestf)
            call ConcatenateProcedures(ProcHook_nsc_ComputeTestf,TestfToZero)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations, DeallocTestf)
         
            call ConcatenateProcedures(ProcHook_nsc_ComputeTestf,ConvectionTestf)
               
            !Convection Jacobian Gradient Matrix in testf
            if (a%kfl_jacgr == 1) then
               call ConcatenateProcedures(ProcHook_nsc_ComputeTestf,ConvectionJacobianGradientTestf)
            endif
      
           !Diffusion Exists
            if (a%kfl_visco == 1) then
                if (abs(a%kfl_stabm) > 0) then
                   call ConcatenateProcedures(ProcHook_nsc_ComputeTestf,DiffusionTestf)
                endif
            endif
            !Reaction Exists
            if (a%kfl_react == 1) then
                if (abs(a%kfl_stabm) > 0) then
                   call ConcatenateProcedures(ProcHook_nsc_ComputeTestf,ReactionTestf)
                endif
            endif

         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
   !-------------------------------------------------------------------
   !Compute Testf values
   subroutine AllocTestf
      implicit none
      
      call a%Memor%alloc(e%mnode,LTdd,'LTdd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,LTdm,'LTdm','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,LTde,'LTde','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,LTmd,'LTmd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,LTmm,'LTmm','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,LTme,'LTme','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,LTed,'LTed','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,LTem,'LTem','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,LTee,'LTee','nsc_elmope_im')

   end subroutine
   
   subroutine DeallocTestf
      implicit none
      
      call a%Memor%dealloc(e%mnode,LTdd,'LTdd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,LTdm,'LTdm','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,LTde,'LTde','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,LTmd,'LTmd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,LTmm,'LTmm','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,LTme,'LTme','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,LTed,'LTed','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,LTem,'LTem','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,LTee,'LTee','nsc_elmope_im')

   end subroutine
   
   subroutine TestfToZero
      implicit none    
      LTdd = 0.0_rp !LTdd(p)
      LTdm = 0.0_rp !LTdm(d,p)
      LTde = 0.0_rp !LTde(p)
      LTmd = 0.0_rp !LTmd(i,p)
      LTmm = 0.0_rp !LTmm(i,d,p)
      LTme = 0.0_rp !LTme(i,p)
      LTed = 0.0_rp !LTed(p)
      LTem = 0.0_rp !LTem(d,p)
      LTee = 0.0_rp !LTee(p)
                        
   end subroutine   

   subroutine ConvectionTestf
      implicit none
   
      !Adjoint Test Function
      !Stabilization terms : -tau A_j*d_j V

      !Mass equation
      LTdd = LTdd + timom(1)*Add 
      LTdm = LTdm + timom(2)*Adm 
      !Momentum equation
      LTmd = LTmd + timom(1)*Amd
      LTmm = LTmm + timom(2)*Amm
      LTme = LTme + timom(3)*Ame
      !Energy equation
      LTed = LTed + timom(1)*Aed  
      LTem = LTem + timom(2)*Aem  
      LTee = LTee + timom(3)*Aee  

   end subroutine
      
   subroutine ConvectionJacobianGradientTestf
      implicit none
   
      !Adjoint Test Function
      !Stabilization terms : -tau d_U(A_j)·d_j(Uh)* V

      !Momentum equation
      LTmd = LTmd + timom(1)*JAmd
      LTmm = LTmm + timom(2)*JAmm
      !Energy equation
      LTed = LTed + timom(1)*JAed  
      LTem = LTem + timom(2)*JAem  
      LTee = LTee + timom(3)*JAee  

   end subroutine
   
   subroutine DiffusionTestf
      implicit none
   
      !Adjoint Test Function
      !Stabilization terms :   +tau d_U(K_kj)·d_k(Uh)*d_j V
      !if NonLinear elements:  +tau K_kj*d^2_kj V 

      !Momentum equation
      LTmd = LTmd + real(a%kfl_stabm)*timom(1)*Kmd
      LTmm = LTmm + real(a%kfl_stabm)*timom(2)*Kmm
      !Energy equation
      LTed = LTed + real(a%kfl_stabm)*timom(1)*Ked  
      LTem = LTem + real(a%kfl_stabm)*timom(2)*Kem  
      LTee = LTee + real(a%kfl_stabm)*timom(3)*Kee  

   end subroutine

   subroutine ReactionTestf
      implicit none
   
      !Adjoint Test Function
      !Stabilization terms :   +tau S V

      !Mass equation
      LTdd = LTdd + real(a%kfl_stabm)*timom(1)*Sdd 
      !Momentum equation
      LTmd = LTmd + real(a%kfl_stabm)*timom(1)*Smd
      LTmm = LTmm + real(a%kfl_stabm)*timom(2)*Smm
      !Energy equation
      LTed = LTed + real(a%kfl_stabm)*timom(1)*Sed  
      LTem = LTem + real(a%kfl_stabm)*timom(2)*Sem  
      LTee = LTee + real(a%kfl_stabm)*timom(3)*See  

   end subroutine

end module 
