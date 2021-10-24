module Mod_lev_ComputeHeightGauges
   use typre
   use Mod_lev_BaseElmope
   use Mod_LevelSet
   implicit none
   private
   public SetPointersComputeHeightGauges
   
   integer(ip), allocatable :: kfl_IsSet
   
   
   real(rp), allocatable :: xglob(:,:)
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersComputeHeightGauges(itask)
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
         
            !--------------------------------------------------------------------
            !Advection (nonsense without convective term)
            if(a%nHeightGauges /=0) then
               call ConcatenateProcedures(ProcHook%Initializations,Initializations)
               call ConcatenateProcedures(ProcHook%PreGauss,PreGauss)
               call ConcatenateProcedures(ProcHook%Finalizations,Finalizations)
            end if
            
            
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   

   subroutine Initializations
      integer(ip) :: ndime, igauge
      
      call a%Mesh%GetNdime(ndime)
      
      
      
      
      
      do igauge = 1,a%nHeightGauges
         call a%HeightGauges(igauge)%ResetGauge
      enddo
      
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,4,xglob,'xglob','EndsteElmope')
   
   end subroutine
   
   subroutine PreGauss
      integer(ip) :: ninters,elemstatus, igauge
   
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      
      if(elemStatus==0)then
         call a%CutMesh%GetNinters(ielem,ninters)
         call a%CutMesh%GetRealIntersectionPoints(ielem,e,xglob)
      
         do igauge = 1,a%nHeightGauges
            call a%HeightGauges(igauge)%CheckHeight(ninters,xglob)
         enddo
         
      endif 
   end subroutine
   
   subroutine Finalizations
      integer(ip) :: igauge
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      
      do igauge = 1,a%nHeightGauges
         call a%HeightGauges(igauge)%CollectHeight
      enddo
      
      call a%Memor%dealloc(ndime,4,xglob,'xglob','EndsteElmope')
   end subroutine





end module