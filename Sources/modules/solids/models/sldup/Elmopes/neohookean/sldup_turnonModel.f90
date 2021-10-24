subroutine sldup_turnonModel_NH(a)
   use typre
   use Mod_UPSolids_NH
   implicit none
   class(UPSolidsProblem_NH) :: a

       a%sld_type = 'NONLI'
       a%sld_model ='NEOHO'

end subroutine sldup_turnonModel_NH
