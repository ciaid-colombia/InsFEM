module Mod_sldAlemov
   use typre   
   use Mod_Alemov
   use Mod_Solids
   implicit none
   private
   public sldAlemovProblem,sldAlemovProblem_Const

   type, extends(AlemovProblem) :: sldAlemovProblem

       !We have a pointer to linear solid statics :)
       class(SolidsProblem), pointer :: solid => NULL()         !Attached solid module (Linear)

       logical :: kfl_printPrincipalStresses= .false. !Used when coupled with SLD to deform mesh smoothly

   contains

   !ALE iteration is solved in specific BEGSTE
   procedure :: SetNdofn          => sldale_SetNdofn
   procedure :: Cvgunk            => sldale_cvgunk
   procedure :: Solite            => sldale_solite
   procedure :: SpecificReaphy    => sldale_reaphy
   procedure :: SpecificBegste    => sldale_begste
   procedure :: SpecificBegite    => sldale_begite
   procedure :: SpecificEndite    => sldale_endite
   procedure :: SpecificEndste    => sldale_endste
   procedure :: deferredTurnon    => sldale_turnon
   procedure :: deferredTurnof    => sldale_turnof
   procedure :: ParticularRestart => sldale_restar
   procedure :: ALESpecificRefine => sldale_refine
   procedure :: sldale_elmope
   
   end type

interface

   subroutine sldale_SetNdofn(a)
      import 
      implicit none
   
      class(sldAlemovProblem) :: a

   end subroutine

    subroutine sldale_reaphy(a,itask)
      use typre
      import 
      implicit none
      class(sldAlemovProblem) :: a
      integer(ip) :: itask
   end subroutine

   subroutine sldale_turnon(a)
      import 
      implicit none
      class(sldAlemovProblem) :: a
   end subroutine

   subroutine sldale_elmope(a,currentbvess)  
      import 
      implicit none
      class(sldAlemovProblem) :: a
      integer :: currentbvess
      
   end subroutine

   subroutine sldale_turnof(a)
      import 
      implicit none
      class(sldAlemovProblem) :: a
   end subroutine

   subroutine sldale_begste(a)
      import 
      implicit none
      class(sldAlemovProblem) :: a
   end subroutine

   subroutine sldale_begite(a)
      import 
      implicit none
      class(sldAlemovProblem) :: a
   end subroutine

   subroutine sldale_refine(a,itask)
      use typre
      import 
      implicit none
      class(sldAlemovProblem) :: a
      character(6) :: itask
   end subroutine

   subroutine sldale_restar(a,itask)
      import 
      implicit none
      class(sldAlemovProblem) :: a
      integer(ip)             :: itask
   end subroutine

   subroutine sldale_endite(a,itask)
      import 
      implicit none
      class(sldAlemovProblem) :: a
      integer(ip) :: itask
   end subroutine

   subroutine sldale_endste(a,itask)
      import 
      implicit none
      class(sldAlemovProblem) :: a
      integer(ip) :: itask
   end subroutine
   
   subroutine sldale_solite(a)
      import 
      implicit none
      class(sldAlemovProblem) :: a
   end subroutine

   subroutine sldale_cvgunk (a,itask)
      import 
      implicit none
      class(sldAlemovProblem) :: a
      integer(ip), intent(in) :: itask
   end subroutine
      
end interface

  interface sldAlemovProblem_Const
      procedure constructor
  end interface sldAlemovProblem_Const

contains

    function constructor()
        class(sldAlemovProblem), pointer :: constructor

        allocate(constructor)

        !Allocate solids problem
        constructor%solid => SolidsProblem_Const()
        constructor%solid%sld_type='LINEA'

    end function constructor

   subroutine sldale_reaMPI(a)
      use typre
      implicit none
      class(sldAlemovProblem) :: a

          call a%solid%SpecificReaMPI
   end subroutine


end module
