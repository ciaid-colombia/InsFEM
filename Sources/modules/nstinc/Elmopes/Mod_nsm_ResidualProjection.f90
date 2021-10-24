  module Mod_nsm_ComputeResidualProjection
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_ComputeGpResidual
   implicit none
   private
   public SetPointersComputeResidualProjection

   type, extends(PointerSetter) :: SPComputeResidualProjection
contains
      procedure :: SpecificSet => SpecificSetComputeResidualProjection
   end type
   type(SPComputeResidualProjection) :: SetPointersComputeResidualProjection
   
   real(rp), allocatable :: elres(:,:)
   real(rp), allocatable :: wrepro(:,:)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SpecificSetComputeResidualProjection(d)
      implicit none
      class(SPComputeResidualProjection) :: d
         
      call ConcatenateProcedures(ProcHook_Initializations,AllocRep)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocRep)
      
      call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroRep)
      
      !We need to compute the residual at the gauss point
      call SetPointersComputeGpResidual%Set
      
      !Now we assembly everything for computing the projection
      call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,GpresToElres)
      
      call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyResidual)
      call ConcatenateProcedures(ProcHook_PostLoop,ProjectResidual)
   end subroutine   
   
   !-------------------------------------------------------------------
   !Residual Projection Computation
   
   subroutine AllocRep
      implicit none
      integer(ip) :: npoin
      
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(a%ResidualSize,npoin,wrepro,'wrepro','nsm_EndElmope')
      call a%Memor%alloc(a%ResidualSize,e%mnode,elres,'elres','nsm_EnditeElmope')
   end subroutine
   
   subroutine DeAllocRep
      implicit none
      call a%Memor%dealloc(a%ResidualSize,e%mnode,elres,'elres','nsm_EnditeElmope')
   end subroutine
   
   subroutine ElmatsToZeroRep
      implicit none
      elres = 0.0_rp
   end subroutine
   
   subroutine GpresToElres
      implicit none
      call nsm_elmrep(e,a%ResidualSize,dvol,gpres,elres)
   end subroutine
   
   subroutine AssemblyResidual
      implicit none
      call a%Mesh%AssemblyToArray(e,a%ResidualSize,elres,wrepro) 
   end subroutine

   subroutine ProjectResidual
      implicit none
      integer(ip) :: npoin
      
      if (a%kfl_SwitchOff==1) then
         call a%Project_disc(a%ResidualSize,wrepro,a%kfl_fixno(1,:)) 
      else
         call a%Project(a%ResidualSize,wrepro) 
      endif
      a%repro = wrepro
      call move_alloc(wrepro,a%repro)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%deallocObj(0,'wrepro','nsm_EndElmope',a%ResidualSize*npoin*rp)
   end subroutine
   
end module
