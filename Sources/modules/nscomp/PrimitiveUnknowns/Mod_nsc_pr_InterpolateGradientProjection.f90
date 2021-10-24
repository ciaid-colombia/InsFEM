module Mod_nsc_pr_InterpolateGradientProjection
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersInterpolateGradientProjection, gprjv, gprjt
   
   real(rp), allocatable :: elgpv(:,:,:), elgpt(:,:)
   real(rp), allocatable :: gprjv(:,:), gprjt(:)

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersInterpolateGradientProjection(itask)
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
         
            if (a%kfl_shock == 2) then
               call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateGradientProjection)
               call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateGradientProjection)
               call ConcatenateProcedures(ProcHook_nsc_pr_Gathers,GatherGradientProjection)
               call ConcatenateProcedures(ProcHook_nsc_pr_Interpolates,InterpolateGradientProjection)
            endif

         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
   
   
   !---------------------------------------------------------------------------
   !Interpolation Subroutines
   subroutine AllocateGradientProjection
      implicit none
      
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,elgpt,'elgpt','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,gprjv,'gprjv','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,gprjt,'gprjt','nsc_pr_elmope')

   end subroutine
   
   subroutine DeallocateGradientProjection
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elgpt,'elgpt','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,gprjv,'gprjv','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,gprjt,'gprjt','nsc_pr_elmope')

   end subroutine
   
   subroutine GatherGradientProjection
      implicit none
      integer(ip)                :: idime

      do idime=1,e%ndime
         elgpv(idime,:,1:e%pnode) = a%grprj(idime,:,e%lnods(1:e%pnode))
      enddo
      elgpt(:,1:e%pnode) = a%grprj(e%ndime+1,:,e%lnods(1:e%pnode))

   end subroutine

   subroutine InterpolateGradientProjection
      implicit none
      integer(ip)                :: idime
      
      do idime=1,e%ndime
         gprjv(idime,:) = matmul(elgpv(idime,:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
      enddo
      gprjt(:) = matmul(elgpt(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
      
   end subroutine
   
end module
