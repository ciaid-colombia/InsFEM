module Mod_sup_buoup0
    use typre
!     use Mod_supm_BaseElmope
    use Mod_supm_BoundaryConditionsHooksAndSubroutines
   implicit none

contains
   subroutine SetBCPointers
      implicit none
      call ResetProcedureComposition

      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetHooksToNULL
      call SetPointersBoundaryConditions(0)
      call SetPointersBoundaryConditions(1)
      call SetPointersBoundaryConditions(100)
   end subroutine
end module

subroutine sup_bouop0(b)
!-----------------------------------------------------------------------
!****f* Nstinc/nsm_bouope
! NAME 
!    nsm_bouop0
! DESCRIPTION
!    Navier-Stokes boundary elemental operations:
!    1. Compute elemental matrix and RHS 
!    2. Impose Dirichlet boundary conditions
!    3. Assemble them
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_NavierStokes
   use Mod_ThreeField
   use Mod_nsm_elmdir
   use Mod_php_elmdir
   use Mod_supm_BaseElmope
   use Mod_sup_buoup0
   implicit none
   class(ThreeFieldNSProblem), target :: b 
  integer(ip) :: nboun,ndime
   
   interface    
      !subroutine nsm_bouope(a) !(*)
      subroutine sup_bouope_Velocity(a)
         use typre
         use Mod_ThreeField
         implicit none
         class(ThreeFieldNSProblem)    :: a
      end subroutine

      subroutine sup_bouope(a)   
         !This subroutine computes the boundary contribution, incompressible Navier-Stokes equations in Three Field Case for Stress
         !-----------------------------------------------------------------------
         use typre       
         use Mod_ThreeField                
         implicit none
         class(ThreeFieldNSProblem) :: a
      end subroutine  
   end interface  
   
      
   a=> b
   
   !Pointers
   call SetBCPointers
 
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
      
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','supm_bouop0')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','supm_bouop0')
   call a%Memor%alloc(e%ndime,e%pnodb,bovel,'bovel','nsm_bouop0')

   !Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !call ProcHookE%DoLoopInitialize
      !Initialize
      elmat=0.0_rp
      elrhs=0.0_rp
      bcstart=(e%ndime-1)*(e%ndime-1)+2
   
      !Compute boundary terms
      !call nsm_bouope(a)  !Velocity conditions (*)
      
      call sup_bouope_Velocity(a) !Velocity conditions
      call sup_bouope(a) ! Stress Conditions

      !Boundary conditions
      call nsm_elmdir(a,e,elmat,elrhs)
      
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   end do boundaries

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','sub_bouelem')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','sub_bouelem')
   call a%Memor%dealloc(e%ndime,e%pnodb,bovel,'bovel','nsm_bouop0')
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   
   !(*)Note: in order to use nsm_boupe we need to use nsm_BaseElmope
   !therefore supm_BaseElmope must use nsm_BaseElmope.
  
   
end subroutine sup_bouop0
