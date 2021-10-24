module Mod_UPSolids_NH
   use typre
   use Mod_UPSolids
   implicit none
   private
   public UPSolidsProblem_NH,UPSolidsProblem_NH_Const

   type, extends(UPSolidsProblem) ::  UPSolidsProblem_NH

   contains

      procedure :: SpecificBegite          => sldup_begite
      procedure :: Elmope                  => sldup_Elmope_NH
      procedure :: EndElmope               => sldup_EndElmope_NH
      procedure :: UPModelTurnon           => sldup_turnonModel_NH
      procedure :: UPGetPhysicalParameters => sldup_GetPhysicalParameters_NH

   end type

   interface

       subroutine sldup_begite(a)
           import 
           implicit none
       class(UPSolidsProblem_NH) :: a      
       end subroutine

       subroutine sldup_turnonModel_NH(a)
           import 
           implicit none
       class(UPSolidsProblem_NH) :: a      

       end subroutine

       subroutine sldup_GetPhysicalParameters_NH(a,ndime,sz,densi,lam,G,ielem)
          import 
          implicit none
          class(UPSolidsProblem_NH), target :: a
          integer(ip),intent(in) :: ndime,sz
          real(rp),   intent(out):: densi,G,lam
          integer(ip),optional   :: ielem
      end subroutine
     
      subroutine sldup_Elmope_NH(a)
         import 
         implicit none
         class(UPSolidsProblem_NH) :: a      
      end subroutine

      subroutine sldup_EndElmope_NH(a,task)
         import 
         implicit none
         class(UPSolidsProblem_NH) :: a
         character(6) :: task
      end subroutine
      
  end interface

  interface UPSolidsProblem_NH_Const
      procedure constructor
  end interface UPSolidsProblem_NH_Const

  contains

      function constructor()
          class(UPSolidsProblem_NH), pointer :: constructor

          allocate(constructor)

      end function constructor

end module Mod_UPSolids_NH

