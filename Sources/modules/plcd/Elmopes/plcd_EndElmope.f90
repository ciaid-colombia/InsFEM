subroutine plcd_EndsteElmope(b)
   use Mod_plcd_BaseElmope
   use Mod_PLCD
   use Mod_plcd_BMatrix
   use Mod_plcd_BMatrixFactory 
   use Mod_plcd_ExacsoError
   implicit none
   class(PLCDProblem), target :: b
   
   real(rp) :: ielem2
   
   a=>b


   
   
   !ToNullsub
   call SetPointersAndHooksToNULLSUB
   !ResetProcedureComposition
   call ResetProcedureComposition
   
   
   !SetPointers
   !Initialize
   call SetPointersExacsoError%Initialize
   
   !Set
   call SetPointersExacsoError%Set
   
   !Finalize
   call SetPointersExacsoError%Finalize
   
   
   
   
  
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_EnditeElmope')
   call a%Mesh%GetNelem(nelem)
   
   call ProcHook%Initializations
   elements : do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      
      ielem2 = ielem
   
      !Move history data to converged
      ElementMatData => a%ElementMaterialsData(ielem)%p
      call ElementMatData%MoveHistoryVariablesToConverged
      
      !Derivatives at the center of gravity, detjm etc
      call e%elmdcg
      
      call ProcHook%PreGauss
      
      do igaus = 1,e%pgaus
         e%igaus = igaus
         
         
         
         call ProcHook%InGaussElmats
         
      enddo
      
      
   enddo elements

   call ProcHook%Finalizations
   
   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')
   
end subroutine
