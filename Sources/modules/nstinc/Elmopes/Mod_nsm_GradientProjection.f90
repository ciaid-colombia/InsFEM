  module Mod_nsm_ComputeGradientProjection
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersComputeGradientProjection

   type, extends(PointerSetter) :: SPComputeGradientProjection
contains
      procedure :: SpecificSet => SpecificSetComputeGradientProjection
   end type
   type(SPComputeGradientProjection) :: SetPointersComputeGradientProjection
 
   real(rp), allocatable :: wgrpro(:,:)
   real(rp), allocatable :: elrhgr(:,:)
   real(rp), allocatable :: elgrv(:,:,:)
   
contains
   
   subroutine SpecificSetComputeGradientProjection(d)
      implicit none
      class(SPComputeGradientProjection) :: d
         
      call ConcatenateProcedures(ProcHook_Initializations,AllocateGradProj)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocateGradProj)
      
      call ConcatenateProcedures(ProcHook_ElmatsToZero,GradProjElmatsToZero)
      
      !Now we assembly everything for computing the projection
      call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,GradProjToElGradProj)
      
      call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyGradProj)
      call ConcatenateProcedures(ProcHook_PostLoop,ProjectGradient)

   end subroutine   
   
   !-------------------------------------------------------------------
   !Gradient Orthogonal Projection calculation 

   subroutine AllocateGradProj
      implicit none
      integer(ip) :: npoin,ngrdf

      call a%Mesh%GetNpoin(npoin)
      ngrdf = (e%ndime)*e%ndime
      call a%Memor%alloc(ngrdf,npoin,wgrpro,'wgrpro','nsm_EndElmope')
      call a%Memor%alloc(ngrdf,e%mnode,elrhgr,'elrhgr','nsm_EndElmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgrv,'elgrv','nsm_EndElmope')
      wgrpro = 0.0_rp

   end subroutine

   subroutine DeallocateGradProj
      implicit none
      integer(ip) :: ngrdf
      
      ngrdf = (e%ndime)*e%ndime
      call a%Memor%dealloc(ngrdf,e%mnode,elrhgr,'elrhgr','nsm_EndElmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgrv,'elgrv','nsm_EndElmope')

   end subroutine

   subroutine GradProjElmatsToZero
      !ElmatsToZero
      elrhgr=0.0_rp
      elgrv=0.0_rp
   end subroutine

   subroutine GradProjToElGradProj
      call grad_proj(e,e%ndime,dvol,grvel,elgrv)
   end subroutine

   subroutine AssemblyGradProj
      implicit none
      integer(ip) :: idime,jdime
      
      do jdime=1,e%ndime      
         do idime=1,e%ndime      
            elrhgr(e%ndime*(idime-1)+jdime,1:e%pnode) = elgrv(idime,jdime,1:e%pnode) &
                  + elrhgr(e%ndime*(idime-1)+jdime,1:e%pnode)
         end do
      end do
      call a%Mesh%AssemblyToArray(e,size(wgrpro,1),elrhgr,wgrpro) 
      
   end subroutine

   subroutine ProjectGradient
      implicit none
      integer(ip) :: idime,jdime
      integer(ip) :: ndime,npoin
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Project(size(wgrpro,1),wgrpro) 
      do idime=1,ndime
         do jdime=1,ndime      
            a%grprj(idime,jdime,:) = wgrpro(ndime*(idime-1)+jdime,:)
         end do
      end do
      call a%Memor%dealloc(ndime*ndime,npoin,wgrpro,'wgrpro','nsm_EndElmope')

   end subroutine
  
end module
