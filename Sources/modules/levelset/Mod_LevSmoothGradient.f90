module Mod_LevSmoothGradient
   use typre
   use Mod_Lev_BaseElmope
   private
   public SetPointersSmoothGradient
   
   integer(ip), allocatable :: kfl_IsSet
   
   real(rp), allocatable :: elrhsSG(:,:)

   

contains
   subroutine SetPointersSmoothGradient(itask)
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
            
            call ConcatenateProcedures(ProcHook%Initializations,SmoothGradientInitializations)
            call ConcatenateProcedures(ProcHook%PreGauss,SmoothGradientPreGauss)
            call ConcatenateProcedures(ProcHook%InGaussElmats,SmoothGradientInGaussElmats)
            call ConcatenateProcedures(ProcHook%PostGaussElmats,SmoothGradientPostGaussElmats)
            call ConcatenateProcedures(ProcHook%Finalizations,SmoothGradientFinalizations)

         
            
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !-------------------------------------------------------------------
   subroutine SmoothGradientInitializations
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elrhsSG,'elrhsSG','SmoothGradientInitializations')
      
      a%SmoothGradient = 0.0_rp
   end subroutine
   
   subroutine SmoothGradientPreGauss
      implicit none
      
      elrhsSG = 0.0_rp
   end subroutine
   
   subroutine SmoothGradientInGaussElmats
      implicit none
      
      integer(ip) :: inode,idime
      real(rp) :: grlev(3)
      
      call e%gradient(1_ip,ellev(:,1),grlev)
      
      do inode = 1,e%pnode
         do idime = 1,e%ndime
            elrhsSG(idime,inode) = elrhsSG(idime,inode) + e%shape(inode,e%igaus)*dvol*grlev(idime)
         enddo
      enddo
   end subroutine   
   
   subroutine SmoothGradientPostGaussElmats
      implicit none
      
      call a%Mesh%AssemblyToArray(e,e%ndime,elrhsSG,a%SmoothGradient)
   end subroutine   
      
   subroutine SmoothGradientFinalizations
      
      call a%Memor%dealloc(e%ndime,e%mnode,elrhsSG,'elrhsSG','SmoothGradientInitializations')
      call a%Mesh%Smooth(e%ndime,a%SmoothGradient)
   
   end subroutine
end module