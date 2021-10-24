module Mod_nsc_InterpolateGradientProjection
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersInterpolateGradientProjection, gprjm, gprje
   
   real(rp), allocatable :: elgpm(:,:,:), elgpe(:,:)
   real(rp), allocatable :: gprjm(:,:), gprje(:)

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
               call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateGradientProjection)
               call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateGradientProjection)
               call ConcatenateProcedures(ProcHook_nsc_Gathers,GatherGradientProjection)
               call ConcatenateProcedures(ProcHook_nsc_Interpolates,InterpolateGradientProjection)
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
      
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgpm,'elgpm','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,elgpe,'elgpe','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,gprjm,'gprjm','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,gprje,'gprje','nsc_elmope_im')

   end subroutine
   
   subroutine DeallocateGradientProjection
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgpm,'elgpm','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,elgpe,'elgpe','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,gprjm,'gprjm','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,gprje,'gprje','nsc_elmope_im')

   end subroutine
   
   subroutine GatherGradientProjection
      implicit none
      integer(ip)                :: idime

      do idime=1,e%ndime
         elgpm(idime,:,1:e%pnode) = a%grprj(idime,:,e%lnods(1:e%pnode))
      enddo
      elgpe(:,1:e%pnode) = a%grprj(e%ndime+1,:,e%lnods(1:e%pnode))

   end subroutine

   subroutine InterpolateGradientProjection
      implicit none
      integer(ip)                :: idime
      
      do idime=1,e%ndime
         gprjm(idime,:) = matmul(elgpm(idime,:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
      enddo
      gprje(:) = matmul(elgpe(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
      
   end subroutine
   
end module
