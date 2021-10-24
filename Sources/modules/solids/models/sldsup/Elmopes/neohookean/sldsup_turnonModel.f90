subroutine sldsup_turnonModel_NH(a)
   use typre
   use Mod_SUPSolids_NH
   implicit none
   class(SUPSolidsProblem_NH) :: a

       a%sld_type = 'NONLI'
       a%sld_model ='NEOHO'

end subroutine sldsup_turnonModel_NH
