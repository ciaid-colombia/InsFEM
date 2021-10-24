 module Mod_nsc_pr_ConfinedFLow
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersConfinedFlow
   
   !External Procedures
   procedure() :: NULLSUB

   !Dynamic subscales
   real(rp) :: timom_static(3)
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersConfinedFlow(itask)
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
            
            !Confined flow
            if (a%kfl_confi == 1) call ConcatenateProcedures(ProcHook_nsc_pr_InGaussElmats,InGaussConfined)

         endif  

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   !----------------------------------------------------------
   !Confined   
   subroutine InGaussConfined
      implicit none
      call nsc_pr_elmbdq_cnf(e,a%epspe,gpden,acvis,dvol,elmdq)
      call nsc_pr_elmrhd_cnf(e,a%epspe,gpden,acvis,dvol,gppre(1),elrhd)
   end subroutine

   subroutine nsc_pr_elmbdq_cnf(e,epspe,acden,acvis,dvol,elmdq)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: acden,acvis,dvol,epspe
      real(rp),    intent(inout) :: elmdq(e%mnode,e%mnode)

      integer(ip)                :: inode,jnode

      forall(inode=1:e%pnode,jnode=1:e%pnode)
            elmdq(inode,jnode) = epspe*acden/acvis*dvol*e%shape(jnode,e%igaus)*e%shape(inode,e%igaus) + elmdq(inode,jnode)
      end forall

   end subroutine nsc_pr_elmbdq_cnf

   subroutine nsc_pr_elmrhd_cnf(e,epspe,acden,acvis,dvolu,gppre,elrhd)
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gppre,acden,acvis,epspe,dvolu
      real(rp),    intent(inout) :: elrhd(e%mnode)
      integer(ip)                :: inode
      
      do inode=1,e%pnode
         elrhd(inode) = elrhd(inode) + epspe*acden/acvis*e%shape(inode,e%igaus)*dvolu
      end do

   end subroutine nsc_pr_elmrhd_cnf
   
end module

 
