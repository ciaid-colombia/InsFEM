module Mod_lmn_ComputeTestf
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_ComputeTaus
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

      integer(ip) :: kfl_nonlinear
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1) 
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            !Allocations
            call ConcatenateProcedures(ProcHook%Initializations,AllocTestf)
            call ConcatenateProcedures(ProcHook%Finalizations,DeallocTestf)
         
            ProcHook%ComputeTestf => ComputeTestf
            call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
            if (kfl_nonlinear == 1) then
               call ConcatenateProcedures(ProcHook%ComputeTestf,ComputeTestfNonlinear)
            endif
            if (a%kfl_stabm == 1) then
               if (a%kfl_tacsg == 1) then
                  if (a%kfl_repro == 0) then
                     call ConcatenateProcedures(ProcHook%Initializations,AllocDSS)
                     call ConcatenateProcedures(ProcHook%Finalizations,DeallocDSS)
                     call ConcatenateProcedures(ProcHook%ComputeTestf,ComputeDSSTestf)
                  endif
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
      
      call a%Memor%alloc(e%ndime,e%mnode,testf_mom,'testf_mom','lmn_elmope')
      call a%Memor%alloc(        e%mnode,testf_ene,'testf_ene','lmn_elmope')
   end subroutine
   
   subroutine DeallocTestf
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,testf_mom,'testf_mom','lmn_elmope')
      call a%Memor%dealloc(        e%mnode,testf_ene,'testf_ene','lmn_elmope')
   end subroutine

   subroutine AllocDSS
      implicit none

      call a%Memor%alloc(e%ndime,e%mnode,testf_mom_DSS,'testf_mom_DSS','lmn_elmope')
      call a%Memor%alloc(        e%mnode,testf_ene_DSS,'testf_ene_DSS','lmn_elmope')
   end subroutine

   subroutine DeallocDSS
      implicit none

      call a%Memor%dealloc(e%ndime,e%mnode,testf_mom_DSS,'testf_mom_DSS','lmn_elmope')
      call a%Memor%dealloc(        e%mnode,testf_ene_DSS,'testf_ene_DSS','lmn_elmope')
   end subroutine
   
   subroutine ComputeTestf
      implicit none
   
      !Adjoint Test Function
      !Stabilization terms : -tau L*v
      call lmn_ComputeTestf_mom(e,acden,timom,AGradN,testf_mom)
      call lmn_ComputeTestf_ene(e,acden,tiene,AGradN,testf_ene)
   end subroutine
      
   subroutine ComputeTestfNonLinear
      implicit none 

      call lmn_ComputeTestfNonLinear_mom(e,acvis,timom,a%kfl_stabm,testf_mom)
      call lmn_ComputeTestfNonLinear_ene(e,actco,tiene,a%kfl_stabm,testf_ene)
   end subroutine

   subroutine ComputeDSSTestf
      implicit none      
      integer(ip)    :: inode,idime

      forall (idime=1:e%ndime, inode=1:e%pnode)
         testf_mom_DSS(idime,inode) = testf_mom(idime,inode) + &
         e%shape(inode,e%igaus)*(timom/timom_static)
      end forall
      testf_ene_DSS(1:e%pnode) = testf_ene(1:e%pnode) + &
      e%shape(1:e%pnode,e%igaus)*(tiene/tiene_static)

      call lmn_ComputeTestfDSS_mom(e,acden,timom,a%dtinv,testf_mom)
      call lmn_ComputeTestfDSS_ene(e,acden,tiene,a%dtinv,testf_ene)
   end subroutine
end module 
