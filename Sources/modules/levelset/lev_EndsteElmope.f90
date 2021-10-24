module Mod_lev_EndsteElmope
   use typre
   use Mod_LevelSet
   use Mod_lev_BaseElmope

contains   
   
   subroutine SetPointers
      use Mod_lev_ComputeHeightGauges
      !External Procedures
      procedure() :: NULLSUB  


      call ResetProcedureComposition
      
      !--------------------------------------------------------------------
      !Defaults
      
      !Pointers
      ProcPointer%ComputeTau => NULLSUB
      ProcPointer%ExternalForces => NULLSUB
      
      !Hooks
      ProcHook%Initializations => NULLSUB
      ProcHook%Gathers         => NULLSUB
      ProcHook%OnIeltyChange   => NULLSUB
      ProcHook%PreGauss        => NULLSUB
      ProcHook%Interpolates    => NULLSUB
      ProcHook%Elext           => NULLSUB
      ProcHook%InGauss         => NULLSUB
      ProcHook%InGaussElmats   => NULLSUB
      ProcHook%PostGaussElmats => NULLSUB
      ProcHook%Testf           => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      ProcHook%PhysicalProp    => NULLSUB
   
      !SetPointers
      call SetPointersComputeHeightGauges(0)
      
      call SetPointersComputeHeightGauges(1)
      
      call SetPointersComputeHeightGauges(100)
   
   
   end subroutine


end module



subroutine lev_EndsteElmope(b)
   use typre
   use Mod_LevelSet
   use Mod_lev_BaseElmope
   use Mod_lev_EndsteElmope
   implicit none
   class(LevelSetProblem), target :: b
   integer(ip) :: igauge,ndime
   
   a=>b 
   
!    call a%Mesh%GetNdime(ndime)
!    a%nHeightGauges = 1
!    allocate(a%HeightGauges(1))
!    a%HeightGauges(1)%Origin(1) = 0.05
!    a%HeightGauges(1)%Origin(2) = 0
!    a%HeightGauges(1)%Origin(3) = 0.1
!    
!    a%HeightGauges(1)%DirectionVector(1) = 0.0
!    a%HeightGauges(1)%DirectionVector(2) = 1.0
!    a%HeightGauges(1)%DirectionVector(3) = 0.0
!    
!    do igauge = 1,a%nHeightGauges
!       call a%HeightGauges(igauge)%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%Mpirank)
!       call a%HeightGauges(igauge)%SetNdime(ndime)
!    enddo
      
   !SetPointers
   call SetPointers
   
   !Initializations
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_elmope')   
   !Hook
   call ProcHook%Initializations
   
   
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)   
      
      !Hook
      call ProcHook%PreGauss
   
   enddo elements
   
   !Hook
   call ProcHook%Finalizations
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_elmope')
   
end subroutine