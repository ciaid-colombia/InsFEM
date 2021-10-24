module Mod_nsm_OrthogonalSubscales
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_nsm_ComputeAdvectionVelocity
   use Mod_nsm_InterpolateResidualProjection
   implicit none
   public

   !For elmopes
   type, extends(PointerSetter) :: SPOrthogonalSubscales
contains
      procedure :: SpecificSet => SpecificSetOrthogonalSubscales
   end type
   type(SPOrthogonalSubscales) :: SetPointersOrthogonalSubscales
   
contains

   subroutine SpecificSetOrthogonalSubscales(d)
      implicit none
      class(SPOrthogonalSubscales) :: d
   
      !ResidualProjection
      if (a%kfl_repro /= 0) then
         !Interpolate Residual Projection
         call SetPointersInterpolateResidualProjection%Set
         
         if (a%kfl_repro == 1) then
            !Matrices
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRep)
            !Skip the finite element part of the residual both in the left and right hand side
            if (a%kfl_repro_SkipFE == 1) then
               ProcPointer_nsm_elmbuv => nsm_elmbuv_trm
               ProcPointer_nsm_elmrhu => nsm_elmrhu_trm
               ProcPointer_nsm_elmrhp => nsm_elmrhp_trm
               ProcPointer_nsm_elmbuq => nsm_elmbuq_trm
               ProcPointer_nsm_elmbpv => nsm_elmbpv
            endif
         elseif (a%kfl_repro == 2) then
            !Matrices
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepSplit)
            !Split OSS   
            ProcPointer_nsm_elmbuv => nsm_elmbuv_split
            ProcPointer_nsm_elmrhu => nsm_elmrhu_split
            ProcPointer_nsm_elmrhp => nsm_elmrhp_split
            
            ProcPointer_nsm_elmbuq => nsm_elmbuq_split
            ProcPointer_nsm_elmbpv => nsm_elmbpv_split
         elseif (a%kfl_repro == 3) then !Split Oss with forces
             !Matrices
             call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepSplit)
             !Split OSS   
             ProcPointer_nsm_elmbuv => nsm_elmbuv_split
             ProcPointer_nsm_elmrhu => nsm_elmrhu_split
             ProcPointer_nsm_elmrhp => nsm_elmrhp_split2
             
             ProcPointer_nsm_elmbuq => nsm_elmbuq_split
             ProcPointer_nsm_elmbpv => nsm_elmbpv_split  
            
         endif
 
         if (a%kfl_FORAxesRotation == 1) then
            !Matrices
            call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepCoriolis)
         endif

      endif
   end subroutine

  !For ResidualProjection
   subroutine InGaussElmatsRep
      use Mod_nsm_InterpolateResidualProjection
      implicit none
      
      !Compute contributions to RHS : Block U
      call nsm_elmrhu_oss(e,tidiv,dvol,testf,gprep,elrhu)
      !Compute contributions to RHS : Block P
      call nsm_elmrhp_oss(e,timom,dvol,gprep,elrhp)     
   end subroutine
   
   subroutine InGaussElmatsRepSplit
      use Mod_nsm_InterpolateResidualProjection
      implicit none    

      !Compute contributions to RHS : Block U
      call nsm_elmrhu_oss(e,tidiv,dvol,testf,gprep,elrhu)
      !Compute contributions to RHS : Block P
      call nsm_elmrhp_oss(e,timom,dvol,gprep(e%ndime+2),elrhp)     
   end subroutine
   
   subroutine InGaussElmatsRepCoriolis
      use Mod_nsm_InterpolateResidualProjection
      implicit none
      
      !Compute contributions to RHS : Block U
      call nsm_elmrhu_oss_coriolis(e,dvol,acden,timom,gprep,a%FORAxesAngularVeloc,elrhu)
   end subroutine
   
end module
