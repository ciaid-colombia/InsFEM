module Mod_lev_EnditeElmope
   use typre
   use Mod_lev_BaseElmope
   use Mod_LevSmoothGradient
   use Mod_Lev_MassCorrection
   implicit none
   
contains

   !SetPointers
   subroutine SetPointers
      use typre
      implicit none
      
      integer(ip) :: kfl_nonlinear,nelty    
      
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
      
      call SetPointersSmoothGradient(0)
      call SetPointersMassCorrection%Initialize
      
      
      if (a%kfl_SmoothGradient == 1) then
         call SetPointersSmoothGradient(1)
      endif
      call SetPointersMassCorrection%Set
      
      
      
      !Finalize setting pointers
      call SetPointersSmoothGradient(100)
      call SetPointersMassCorrection%Finalize
      
      
   end subroutine

end module


subroutine lev_EnditeElmope(b)
   use typre
   use Mod_LevelSet
   use Mod_lev_BaseElmope
   use Mod_lev_EnditeElmope
   implicit none
   class(LevelSetProblem), target :: b
   
   logical :: DoSomething 
   
   a=>b
   
   DoSomething = .false.
   if (a%kfl_SmoothGradient == 1) DoSomething = .true.
   if (DoSomething .eqv. .false.) return

   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Allocate elements
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_elmope')
   
   !Allocation of arrays
   call a%Memor%alloc(e%mnode,a%ncomp-1,ellev,'ellev','lev_elmope')
   call a%Memor%alloc(a%ncomp-1,gplev,'gplev','lev_elmope')   
   
   
   !Initializations
   call ProcHook%Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem) 
   
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e) 
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
      
      call e%gather(1,ellev(:,1),a%level(:,1))
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss
         
         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call e%interpg(1,ellev(:,1),gplev(1))
         
         !Hook
         call ProcHook%InGaussElmats
      enddo gauss_points
      
      call ProcHook%PostGaussElmats
      
      
   enddo elements
   
   !Finalizations
   !Hook
   call ProcHook%Finalizations
   
   !DeAllocation of arrays
   call a%Memor%dealloc(a%ncomp-1,gplev,'gplev','lev_elmope')   
   call a%Memor%dealloc(e%mnode,a%ncomp-1,ellev,'ellev','lev_elmope')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_elmope')
         
         

end subroutine
