  module Mod_nsc_ComputeGradientProjection
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersComputeGradientProjection

   real(rp), allocatable :: wgrpro(:,:)
   real(rp), allocatable :: elrhgr(:,:)
   real(rp), allocatable :: elgrd(:,:), elgrm(:,:,:), elgre(:,:)
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeGradientProjection(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateGradProj)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateGradProj)
         
            call ConcatenateProcedures(ProcHook_nsc_ElmatsToZero,GradProjElmatsToZero)
            
            !Now we assembly everything for computing the projection
            call ConcatenateProcedures(ProcHook_nsc_InGaussElmatsAssembly,GradProjToElGradProj)
            
            call ConcatenateProcedures(ProcHook_nsc_AssemblyEndite,AssemblyGradProj)
            call ConcatenateProcedures(ProcHook_nsc_PostLoop,ProjectGradient)
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !-------------------------------------------------------------------
   !Gradient Orthogonal Projection calculation 

   subroutine AllocateGradProj
      implicit none
      
      integer(ip) :: npoin,ngrdf
      call a%Mesh%GetNpoin(npoin)
     
      ngrdf = (e%ndime+2)*e%ndime
      call a%Memor%alloc(ngrdf,npoin,wgrpro,'wgrpro','nsc_EndElmope')
      call a%Memor%alloc(ngrdf,e%mnode,elrhgr,'elrhgr','nsc_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elgrd,'elgrd','nsc_EndElmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgrm,'elgrm','nsc_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elgre,'elgre','nsc_EndElmope')

      wgrpro = 0.0_rp

   end subroutine

   subroutine DeallocateGradProj
      implicit none
      integer(ip) :: ngrdf
      
      ngrdf = (e%ndime+2)*e%ndime
      call a%Memor%dealloc(ngrdf,e%mnode,elrhgr,'elrhgr','nsc_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elgrd,'elgrd','nsc_EndElmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgrm,'elgrm','nsc_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elgre,'elgre','nsc_EndElmope')

   end subroutine

   !----------------------------------------------------------
   subroutine GradProjElmatsToZero
      
      !ElmatsToZero
      elrhgr=0.0_rp
      elgrd=0.0_rp
      elgrm=0.0_rp
      elgre=0.0_rp

   end subroutine

   !----------------------------------------------------------

   subroutine GradProjToElGradProj
      
       call grad_proj(e,1_ip,dvol,grden,elgrd)
       call grad_proj(e,e%ndime,dvol,grmom,elgrm)
       call grad_proj(e,1_ip,dvol,grene,elgre)
      
   end subroutine

   !----------------------------------------------------------
   subroutine AssemblyGradProj
      implicit none
      integer(ip) :: idime,jdime
      
      do jdime=1,e%ndime      
         do idime=1,e%ndime      
            elrhgr(e%ndime*(idime-1)+jdime,1:e%pnode) = elgrm(idime,jdime,1:e%pnode) &
                  + elrhgr(e%ndime*(idime-1)+jdime,1:e%pnode)
         end do
         elrhgr(e%ndime*e%ndime+jdime,1:e%pnode) = elgre(jdime,1:e%pnode) &
              + elrhgr(e%ndime*e%ndime+jdime,1:e%pnode)
         elrhgr(e%ndime*(e%ndime+1)+jdime,1:e%pnode) = elgrd(jdime,1:e%pnode) &
              + elrhgr(e%ndime*(e%ndime+1)+jdime,1:e%pnode)
      end do
 
 
      !GlobalAssembly

      call a%Mesh%AssemblyToArray(e,size(wgrpro,1),elrhgr,wgrpro) 
      
   end subroutine

   !----------------------------------------------------------
   subroutine ProjectGradient
      implicit none
      integer(ip) :: idime,jdime
      integer(ip) :: ndime,npoin
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)

      call a%Project(size(wgrpro,1),wgrpro) 
      
      do idime=1,ndime+1      
         do jdime=1,ndime      
            a%grprj(idime,jdime,:) = wgrpro(ndime*(idime-1)+jdime,:)
         end do
      end do
      do jdime=1,ndime      
         a%dgrprj(jdime,:) = wgrpro(ndime*(ndime+1)+jdime,:)
      end do
      
      !Dealloc wrepro
      call a%Memor%dealloc((ndime+2)*ndime,npoin,wgrpro,'wgrpro','nsc_EndElmope')
      
   end subroutine
  
end module


 
