module Mod_tem_TempeGradientProjection
   use typre
   use Mod_tem_BaseElmope
   use Mod_tem_TempeGradient
   implicit none
   private
   public SetPointersComputeTempeGradientProjection

   
   integer(ip), allocatable :: kfl_IsSet
   
   real(rp), allocatable :: elgrtem(:,:)
   integer(ip) :: ndime
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeTempeGradientProjection(itask)
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
   
            call ConcatenateProcedures(ProcHook%Initializations,GradientProjectionToZero)
            call SetPointersComputeTempeGradient(1)
            call ConcatenateProcedures(ProcHook%PreGauss,GradientProjectionElgrtemToZero)
            call ConcatenateProcedures(ProcHook%InGaussElmats,GradientProjectionAssemblyToElgrtem)
            call ConcatenateProcedures(ProcHook%Assembly,GradientProjectionAssembly)
            call ConcatenateProcedures(ProcHook%Finalizations,GradientProjectionFinalizations)
         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
         
   end subroutine   
   
  !--------------------------------------------------------------------
   !Temperature Gradient interpolation
   subroutine GradientProjectionToZero
      implicit none
      
      call a%Mesh%GetNdime(ndime)
      
      a%gradient = 0.0_rp
      call a%Memor%alloc(ndime,e%mnode,elgrtem,'elgrtem','Tem_GradientProjection')
   end subroutine
   
   subroutine GradientProjectionElgrtemToZero
      implicit none
      
      elgrtem = 0.0_rp
   end subroutine
   
   subroutine GradientProjectionAssemblyToElGrtem
      implicit none
      
      do inode = 1,e%pnode
         elgrtem(:,inode) = elgrtem(:,inode) + e%shape(inode,e%igaus)*grtem*dvol
      enddo
      
   end subroutine
   
   subroutine GradientProjectionAssembly
      implicit none
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%AssemblyToArray(e,ndime,elgrtem,a%gradient)
   end subroutine
   
   
   
   subroutine GradientProjectionFinalizations
      implicit none
      
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%Project(ndime,a%gradient)
      
      call a%Memor%dealloc(ndime,e%mnode,elgrtem,'elgrtem','Tem_GradientProjection')
   end subroutine
   
end module 
