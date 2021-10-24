  module Mod_tem_ComputeGradientProjection
   use typre
   use Mod_tem_BaseElmope
   implicit none
   private
   public SetPointersComputeGradientProjection

   real(rp), allocatable :: wgrpro(:,:)
   real(rp), allocatable :: elrhgr(:,:)
   real(rp), allocatable :: elgrt(:,:),grtem(:,:)
   
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
         
            call ConcatenateProcedures(ProcHook%Initializations,AllocateGradProj)
            call ConcatenateProcedures(ProcHook%Finalizations,DeallocateGradProj)
         
            call ConcatenateProcedures(ProcHook%ElmatsToZero,GradProjElmatsToZero)
            
            !Now we assembly everything for computing the projection
            call ConcatenateProcedures(ProcHook%InGaussElmats,GradProjToElGradProj)
            
            call ConcatenateProcedures(ProcHook%Assembly,AssemblyGradProj)
            call ConcatenateProcedures(ProcHook%PostLoop,ProjectGradient)

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
     
      ngrdf = e%ndime
      call a%Memor%alloc(ngrdf,npoin,wgrpro,'wgrpro','tem_EndElmope')
      call a%Memor%alloc(ngrdf,e%mnode,elrhgr,'elrhgr','tem_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elgrt,'elgrt','tem_EndElmope')
      call a%Memor%alloc(1_ip   ,e%ndime,grtem,'grtem','tem_EndElmope')

      wgrpro = 0.0_rp

   end subroutine

   subroutine DeallocateGradProj
      implicit none
      integer(ip) :: ngrdf
      
      ngrdf = e%ndime
      call a%Memor%dealloc(ngrdf,e%mnode,elrhgr,'elrhgr','tem_EndElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elgrt,'elgrt','tem_EndElmope')
      call a%Memor%dealloc(1_ip   ,e%ndime,grtem,'grtem','tem_EndElmope')

   end subroutine

   !----------------------------------------------------------
   subroutine GradProjElmatsToZero
      
      !ElmatsToZero
      elrhgr=0.0_rp
      elgrt=0.0_rp

   end subroutine

   !----------------------------------------------------------

   subroutine GradProjToElGradProj
      
       call grad_proj(e,1_ip,dvol,grtem,elgrt)
      
   end subroutine

   !----------------------------------------------------------
   subroutine AssemblyGradProj
      implicit none
      integer(ip) :: idime,jdime
      
      do jdime=1,e%ndime      
         elrhgr(jdime,1:e%pnode) = elgrt(jdime,1:e%pnode) &
              + elrhgr(jdime,1:e%pnode)
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
      
      do jdime=1,ndime      
         a%grprj(jdime,:) = wgrpro(jdime,:)
      end do
      
      !Dealloc wrepro
      call a%Memor%dealloc(ndime,npoin,wgrpro,'wgrpro','tem_EndElmope')
      
   end subroutine
  
end module


 
