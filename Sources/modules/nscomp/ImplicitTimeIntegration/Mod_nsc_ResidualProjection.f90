  module Mod_nsc_ComputeResidualProjection
   use typre
   use Mod_nsc_BaseElmope
   use Mod_nsc_ComputeGpResidual
   implicit none
   private
   public SetPointersComputeResidualProjection

   real(rp), allocatable :: elres(:,:)
   real(rp), allocatable :: wrepro(:,:)
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeResidualProjection(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateRep)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateRep)
         
            call ConcatenateProcedures(ProcHook_nsc_ElmatsToZero,RepElmatsToZero)
            
            !We need to compute the residual at the gauss point
            call SetPointersComputeGpResidual(1)
            
            !Now we assembly everything for computing the projection
            call ConcatenateProcedures(ProcHook_nsc_InGaussElmatsAssembly,GpresToElres)
            
            call ConcatenateProcedures(ProcHook_nsc_AssemblyEndite,AssemblyResidual)
            call ConcatenateProcedures(ProcHook_nsc_PostLoop,ProjectResidual)
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !-------------------------------------------------------------------
   !Residual Projection Computation
   
   !Residual Projection 
   
   subroutine AllocateRep
      implicit none
      
      integer(ip) :: ndime,npoin
      
      !We allocate wrepro for the computation of the residual
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      
      call a%Memor%alloc(a%ndofn,npoin,wrepro,'wrepro','nsc_EndElmope_im')
      call a%Memor%alloc(a%ndofn,e%mnode,elres,'elres','nsc_EndElmope_im')

      wrepro(:,:) = 0.0_rp

   end subroutine
   
   subroutine DeallocateRep
      implicit none
      call a%Memor%dealloc(a%ndofn,e%mnode,elres,'elres','nsc_EndElmope_im')
   end subroutine
   
   subroutine RepElmatsToZero
      implicit none
      
      elres = 0.0_rp
   end subroutine
   
   subroutine GpresToElres
      implicit none
      
      integer(ip) :: idime
   
      elres(1,1:e%pnode) = elres(1,1:e%pnode) + e%shape(1:e%pnode,e%igaus)*gpred(1)*dvol
      do idime = 1,e%ndime
         elres(idime+1,1:e%pnode) = elres(idime+1,1:e%pnode) + e%shape(1:e%pnode,e%igaus)*gprem(idime)*dvol
      enddo
      elres(e%ndime+2,1:e%pnode) = elres(e%ndime+2,1:e%pnode) + e%shape(1:e%pnode,e%igaus)*gpree(1)*dvol

   end subroutine
   
   subroutine AssemblyResidual
      implicit none
      
      call a%Mesh%AssemblyToArray(e,size(wrepro,1),elres,wrepro) 
   end subroutine

   subroutine ProjectResidual
      implicit none
      
      integer(ip) :: ndime,npoin
      
      a%repro = wrepro
      call move_alloc(wrepro,a%repro)
      
      !Dealloc wrepro
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%deallocObj(0,'wrepro','nsc_EndElmope_im',a%ndofn*npoin*rp)
      
   end subroutine
   
end module


 
