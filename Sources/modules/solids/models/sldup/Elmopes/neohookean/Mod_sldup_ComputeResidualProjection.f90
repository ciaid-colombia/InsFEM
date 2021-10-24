  module Mod_sldup_ComputeResidualProjection
   use typre
   use Mod_sld_BaseElmope
   use Mod_sldup_ComputeGpResidual
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
         
            call ConcatenateProcedures(ProcHook%Initializations,AllocRep)
            call ConcatenateProcedures(ProcHook%Finalizations,DeallocRep)
         
            call ConcatenateProcedures(ProcHook%ElmatsToZero,ElmatsToZeroRep)
            
            !We need to compute the residual at the gauss point
            call SetPointersComputeGpResidual(1)
            
            call ConcatenateProcedures(ProcHook%Initializations,resProAllocate)

            !Now we assemble everything for computing the projection
            call ConcatenateProcedures(ProcHook%EndGauss,GpresToElres)
            
            call ConcatenateProcedures(ProcHook%ToLinearSystem,AssemblyResidual)

            call ConcatenateProcedures(ProcHook%Finalizations,ProjectResidual)
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !-------------------------------------------------------------------
   !Residual Projection Computation
   subroutine resProAllocate
      implicit none
      integer(ip) :: npoin
    
      !We allocate wrepro for the computation of the residual
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(a%ndofn,npoin,wrepro,'wrepro','sldup_EndElmope')
   end subroutine
   
   subroutine AllocRep
      implicit none
      call a%Memor%alloc(a%ndofn,e%mnode,elres,'elres','sldup_EnditeElmope')
   end subroutine
   
   subroutine DeallocRep
      implicit none
      call a%Memor%dealloc(a%ndofn,e%mnode,elres,'elres','sldup_EnditeElmope')
   end subroutine
   
   subroutine ElmatsToZeroRep
      implicit none
      
      elres = 0.0_rp
   end subroutine
   
   subroutine GpresToElres
      implicit none
      integer(ip) :: inode
   
      do inode = 1,e%pnode
         elres(:,inode) = elres(:,inode) + e%shape(inode,e%igaus)*gpres*dvol
      enddo
   end subroutine
   
   subroutine AssemblyResidual
      implicit none
      
      call a%Mesh%AssemblyToArray(e,size(wrepro,1),elres,wrepro)

   end subroutine

   subroutine ProjectResidual
      implicit none
      integer(ip) :: npoin
      
      call up%Project(a%ndofn,wrepro) 
      
      up%repro = wrepro
      call move_alloc(wrepro,up%repro)
      
      !Dealloc wrepro
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%deallocObj(0,'wrepro','sldup_EndElmope',(a%ndofn)*npoin*rp)
      
   end subroutine
   
end module

