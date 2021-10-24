module Mod_SUPSolids_lin
   use typre
   use Mod_SUPSolids
   implicit none
   private
   public SUPSolidsProblem_lin,SUPSolidsProblem_lin_Const

   type, extends(SUPSolidsProblem) ::  SUPSolidsProblem_lin

   contains

      !procedure :: SpecificBegite           => sldsup_NULLSUB
      procedure :: Elmope                   => sldsup_Elmope_lin
      procedure :: EndElmope                => sldsup_EndElmope_lin
      procedure :: SUPGetPhysicalParameters => sldsup_GetPhysicalParameters_lin

   end type

   interface

      subroutine sldsup_GetPhysicalParameters_lin(a,ndime,sz,densi,K,G,C_dev,D,D_dev,D_comp,P,ielem)
         import 
         implicit none
         class(SUPSolidsProblem_lin), target :: a
         integer(ip),intent(in) :: ndime,sz
         real(rp),   intent(out):: densi,K,G,D_comp
         real(rp),   intent(out):: C_dev(sz,sz),D_dev(sz,sz),D(sz,sz),P(sz,sz)
         integer(ip),optional   :: ielem
      end subroutine
     
      subroutine sldsup_Elmope_lin(a)
         import 
         implicit none
         class(SUPSolidsProblem_lin) :: a      
      end subroutine

      subroutine sldsup_EndElmope_lin(a,task)
         import 
         implicit none
         class(SUPSolidsProblem_lin) :: a
         character(6) :: task
      end subroutine
      
  end interface

  interface SUPSolidsProblem_lin_Const
      procedure constructor
  end interface SUPSolidsProblem_lin_Const

  contains

      function constructor()
          class(SUPSolidsProblem_lin), pointer :: constructor

          allocate(constructor)

      end function constructor

   subroutine sldsup_NULLSUB(a)
      use typre
      implicit none
      class(SUPSolidsProblem_lin) :: a
   
   end subroutine


end module Mod_SUPSolids_lin

