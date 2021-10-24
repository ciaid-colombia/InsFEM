  module Mod_supm_ComputeVorticity
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateGradients
   implicit none
   private
   public SetPointersComputeVorticity

   !Vorticity
   real(rp), allocatable :: elvort(:,:), gpvort(:)
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeVorticity(itask)
      implicit none
      integer(ip) :: itask
      
      logical     :: aux_logic
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            !Logics for deciding if we need to compute vorticity
            aux_logic = .false.
            if (a%npp_stepi(10) /= 0) then   
               if (mod(a%istep,a%npp_stepi(10))==0) then
                  aux_logic = .true.
               endif
            endif
         
            !If vorticity needs to be computed
            if (aux_logic .eqv. .true.) then
               call SetPointersInterpolateGradients(1)
               call ConcatenateProcedures(ProcHook_Initializations,AllocVorti)
               call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeVorticity)
               call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyVorticity)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocVorti)
               call ConcatenateProcedures(ProcHook_PostLoop,SmoothVorticity)
            endif  
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !Vorticity 
   subroutine AllocVorti
      implicit none
      
      call a%Memor%alloc(3,e%mnode,elvort,'elvort','nsm_EndElmope')
      call a%Memor%alloc(3,gpvort,'gpvort','nsm_EndElmope')
      a%vorti(:,:) = 0_rp
   end subroutine
   
   subroutine ComputeVorticity
      implicit none
      integer(ip)  :: inode

      gpvort(3) = grvel(2,1)-grvel(1,2)
      
      if (e%ndime == 3) then
         gpvort(1) = grvel(3,2)-grvel(2,3)
         gpvort(2) = grvel(1,3)-grvel(3,1)
      endif      
      do inode=1,e%pnode
         elvort(:,inode)=elvort(:,inode)+e%shape(inode,e%igaus)*gpvort*dvol
      end do
   end subroutine

   subroutine AssemblyVorticity
      implicit none
      integer(ip)  :: inode
      
      call a%Mesh%AssemblyToArray(e,size(elvort,1),elvort,a%vorti) 
      elvort(:,:)=0_rp
   end subroutine

   subroutine SmoothVorticity
      implicit none
      
      call a%Mesh%Smooth(3,a%vorti)
      
   end subroutine

   subroutine DeallocVorti
      implicit none

      call a%Memor%dealloc(3,e%mnode,elvort,'elvort','nsm_EndElmope')
      call a%Memor%dealloc(3,gpvort,'gpvort','nsm_EndElmope')
   end subroutine
end module


 
