module Mod_SUPSolids_NH
   use typre
   use Mod_SUPSolids
   implicit none
   private
   public SUPSolidsProblem_NH,SUPSolidsProblem_NH_Const

   type, extends(SUPSolidsProblem) ::  SUPSolidsProblem_NH

   contains

      procedure :: SpecificBegite           => sldsup_begite
      procedure :: Elmope                   => sldsup_Elmope_NH
      procedure :: EndElmope                => sldsup_EndElmope_NH
      procedure :: SUPModelTurnon           => sldsup_turnonModel_NH
      procedure :: SUPGetPhysicalParameters => sldsup_GetPhysicalParameters_NH

   end type

   interface

       subroutine sldsup_begite(a)
           import 
           implicit none
       class(SUPSolidsProblem_NH) :: a      
       end subroutine

       subroutine sldsup_turnonModel_NH(a)
           import 
           implicit none
       class(SUPSolidsProblem_NH) :: a      

       end subroutine

       subroutine sldsup_GetPhysicalParameters_NH(a,ndime,sz,densi,lam,G,ielem)
          import 
          implicit none
          class(SUPSolidsProblem_NH), target :: a
          integer(ip),intent(in) :: ndime,sz
          real(rp),   intent(out):: densi,G,lam
          integer(ip),optional   :: ielem
      end subroutine
     
      subroutine sldsup_Elmope_NH(a)
         import 
         implicit none
         class(SUPSolidsProblem_NH) :: a      
      end subroutine

      subroutine sldsup_EndElmope_NH(a,task)
         import 
         implicit none
         class(SUPSolidsProblem_NH) :: a
         character(6) :: task
      end subroutine
      
  end interface

  interface SUPSolidsProblem_NH_Const
      procedure constructor
  end interface SUPSolidsProblem_NH_Const

  contains

      function constructor()
          class(SUPSolidsProblem_NH), pointer :: constructor

          allocate(constructor)

      end function constructor

end module Mod_SUPSolids_NH

