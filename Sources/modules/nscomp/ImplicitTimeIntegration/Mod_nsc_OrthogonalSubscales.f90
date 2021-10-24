module Mod_nsc_OrthogonalSubscales
   use typre
   use Mod_nsc_BaseElmope
   use Mod_nsc_InterpolateResidualProjection
   implicit none
   private
   public SetPointersOrthogonalSubscales
   
   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersOrthogonalSubscales(itask)
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
         
            !ResidualProjection
            if (a%kfl_repro /= 0) then
               !Interpolate Residual Projection
               call SetPointersInterpolateResidualProjection(1)
               
               if (a%kfl_repro == 1) then
                  !Matrices
                  call ConcatenateProcedures(ProcHook_nsc_InGaussElmats,InGaussElmatsRep)
               endif
            endif

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine

  !For ResidualProjection
   subroutine InGaussElmatsRep
      implicit none

      !Compute contributions to RHS : Block pressure
      call nsc_elmrhd_oss(e,dvol,gprep,LTdd,LTdm,LTde,elrhd)     
      
      !Compute contributions to RHS : Block velocity
      call nsc_elmrhm_oss(e,dvol,gprep,LTmd,LTmm,LTme,elrhm)
      
      !Compute contributions to RHS : Block temperature
      call nsc_elmrhe_oss(e,dvol,gprep,LTed,LTem,LTee,elrhe)     

   end subroutine
   
   subroutine nsc_elmrhd_oss(e,dvolu,gprep,LTdd,LTdm,LTde,elrhd)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for mass equation
      !    - Tau_d(q·L*dd, proj(Resd)) 
      !    - Tau_m(q.L*dm, proj(Resm))
      !    - Tau_e(q.L*de, proj(Rese))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gprep(e%ndime+2)
      real(rp),    intent(in)    :: LTdd(e%mnode)
      real(rp),    intent(in)    :: LTdm(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTde(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhd(e%mnode)

      real(rp)                   :: tmp(e%ndime)
      real(rp)                   :: aux(e%mnode)

      tmp(1:e%ndime) = gprep(2:e%ndime+1)
      aux = matmul(tmp,LTdm)

      elrhd(1:e%pnode) = elrhd(1:e%pnode) - &
       (LTdd(1:e%pnode)*gprep(1) + &
        aux(1:e%pnode) + LTde(1:e%pnode)*gprep(e%ndime+2))*dvolu 

   end subroutine nsc_elmrhd_oss

   subroutine nsc_elmrhm_oss(e,dvolu,gprep,LTmd,LTmm,LTme,elrhm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for momentum equation
      !    - Tau_d(n·L*md, proj(Resd)) 
      !    - Tau_m(n.L*mm, proj(Resm))
      !    - Tau_e(n.L*me, proj(Rese))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gprep(e%ndime+2)
      real(rp),    intent(in)    :: LTmd(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTmm(e%ndime,e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTme(e%ndime,e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhm(e%ndime,e%mnode)

      integer(ip)                :: idime
      real(rp)                   :: tmp(e%ndime)
      real(rp)                   :: aux(e%mnode)

      tmp(1:e%ndime) = gprep(2:e%ndime+1)

      do idime=1,e%ndime
         aux(:) =  matmul(transpose(LTmm(idime,:,:)),tmp) 
         elrhm(idime,1:e%pnode) = elrhm(idime,1:e%pnode) - &
          (LTmd(idime,1:e%pnode)*gprep(1) +&
           LTme(idime,1:e%pnode)*gprep(e%ndime+2) + &
           aux(1:e%pnode))*dvolu 
      end do

   end subroutine nsc_elmrhm_oss

   subroutine nsc_elmrhe_oss(e,dvolu,gprep,LTed,LTem,LTee,elrhe)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for energy
      !    - Tau_d(g·L*ed, proj(Resd)) 
      !    - Tau_m(g.L*em, proj(Resm))
      !    - Tau_e(g.L*ee, proj(Rese))
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gprep(e%ndime+2)
      real(rp),    intent(in)    :: LTed(e%mnode)
      real(rp),    intent(in)    :: LTem(e%ndime,e%mnode)
      real(rp),    intent(in)    :: LTee(e%mnode)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhe(e%mnode)

      real(rp)                   :: tmp(e%ndime)
      real(rp)                   :: aux(e%mnode)

      tmp(1:e%ndime) = gprep(2:e%ndime+1)
      aux =  matmul(tmp,LTem) 

      elrhe(1:e%pnode) = elrhe(1:e%pnode) - &
       (LTee(1:e%pnode)*gprep(e%ndime+2) + &
       LTed(1:e%pnode)*gprep(1) + aux(1:e%pnode))*dvolu 

   end subroutine nsc_elmrhe_oss

end module
