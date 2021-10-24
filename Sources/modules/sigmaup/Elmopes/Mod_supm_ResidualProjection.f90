module Mod_supm_ComputeResidualProjection
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_ComputeGpResidual
   use Mod_supm_InterpolateGradients
   implicit none
   private
   public SetPointersComputeResidualProjection
   real(rp), allocatable :: wrepro(:,:)
   integer(ip), allocatable :: kfl_IsSet
   integer(ip) :: kfl_nonlinear
 
contains
   
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

            call ConcatenateProcedures(ProcHook_Initializations,InitializationsRep)
            call ConcatenateProcedures(ProcHook_Finalizations,FinalizationsRep)
            
            if (a%kfl_repro >= 1 .and. a%kfl_repro/=4) then
               call ConcatenateProcedures(ProcHook_PreLoop,PreLoopRep)
               call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroRep)
            
               !We need to compute the residual at the gauss point
               call SetPointersComputeGpResidual(1)

               !Now we assembly everything for computing the projection
               call ConcatenateProcedures(ProcHook_InGaussElmats,GaussPointAssemblyRep)         
               call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyResidual)     
               
               if(a%MatProp(imat)%lawvi<0)then
                  call ConcatenateProcedures(ProcHook_PostLoop,ResidualBoundary)
               endif
               
               call ConcatenateProcedures(ProcHook_PostLoop,ProjectResidual)
            endif
         endif  

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select
      
   end subroutine  
   
!------------------------------------------------
!SUBROUTINES

   !------------------------------------------
   !Specific subroutines for supm_EnditeElmope
   !------------------------------------------
   subroutine PreLoopRep
      implicit none
       a%repro = 0.0_rp
   end subroutine
   
   subroutine InitializationsRep
      implicit none
      !Matrices alloc
      call a%Memor%alloc(a%ndofn,e%mnode,elres,'elres','supm_EnditeElmope')
   end subroutine

   subroutine ElmatsToZeroRep
      implicit none
        elres = 0.0_rp
   end subroutine
    
   subroutine FinalizationsRep
      implicit none
      !Matrices dealloc
      call a%Memor%dealloc(a%ndofn,e%mnode,elres,'elres','supm_EnditeElmope')
   end subroutine
 
   subroutine GaussPointAssemblyRep
      implicit none
      integer(ip) :: ntens,elemStatus     
      !this is crucial when we have two fluids and enriched pressure
      if(a%kfl_colev==1)then
         call a%CutMesh%GetElementType(ielem,elemStatus)
         if(elemStatus ==0)then 
            gpres=0.0_rp
         end if         
      end if
      !Residual gauss point to elemental residual
      call supm_elmrep(e,dvol,a%ndofn,gpres,elres)
   end subroutine
   
   subroutine AssemblyResidual
      implicit none
      a%repro(:,e%lnods(1:e%pnode)) = a%repro(:,e%lnods(1:e%pnode)) + elres(:,1:e%pnode)    
   end subroutine
   
   subroutine ResidualBoundary
      call a%Mesh%Getnpoin(npoin)      
      do ipoin=1,npoin      
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if(ibopo>0 .and. a%kfl_reproBoundzero==1) then
            do idime=1,auxtens+e%ndime+1
               a%repro(idime,ipoin) = 0.0_rp
            end do           
         end if
      end do    
   end subroutine

   subroutine ProjectResidual
      implicit none
      call a%Project(a%ndofn,a%repro) 
   end subroutine
   
end module
