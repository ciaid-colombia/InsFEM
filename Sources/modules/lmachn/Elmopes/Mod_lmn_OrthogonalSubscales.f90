module Mod_lmn_OrthogonalSubscales
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_InterpolateGradients
   use Mod_lmn_InterpolateResidualProjection
   implicit none
   private
   public SetPointersOrthogonalSubscales
   
   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersOrthogonalSubscales(itask)
      implicit none
      procedure() :: NULLSUB
      integer(ip) :: itask
   
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         !ResidualProjection
         if (a%kfl_repro /= 0) then
            !Interpolate Residual Projection
            call SetPointersInterpolateResidualProjection(1)
            
            !Matrices
            call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsRep)
            
            !Skip the finite element part of the residual both in the left and right hand side
            ProcPointer%lmn_elmbuv => lmn_elmbuv_trm
            ProcPointer%lmn_elmbuq => lmn_elmbuq_trm
            ProcPointer%lmn_elmbtw => lmn_elmbtw_trm
            ProcPointer%lmn_elmrhu => lmn_elmrhu_trm
            ProcPointer%lmn_elmrht => lmn_elmrht_trm
            ProcPointer%TimeIntegrationToElext => NULLSUB
         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine

  !For ResidualProjection
   subroutine InGaussElmatsRep
      implicit none
      
      !Compute contributions to RHS : Block U
      call lmn_elmrhu_oss(e,ticon,dvol,testf_mom,gprep(1:e%ndime+1),elrhu)
      !Compute contributions to RHS : Block T
      call lmn_elmrht_oss(e,dvol,testf_ene,gprep(e%ndime+2),elrht)
      !Compute contributions to RHS : Block P
      call lmn_elmrhp_oss(e,timom,dvol,acden,gprep(1:e%ndime),elrhp) 
   end subroutine

end module
