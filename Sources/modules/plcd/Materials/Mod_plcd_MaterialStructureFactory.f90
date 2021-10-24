module Mod_plcd_MaterialFactory
   use typre
   use Mod_plcd_LinearElastic
   use Mod_plcd_Material
   use Mod_plcd_DamageMaterial
   use Mod_plcd_CompositeMaterial
   use Mod_plcd_SIMPMaterial
   use Mod_plcd_SVKirchhoffMaterial
   use Mod_plcd_NeoHookeanUpdatedFormulation
   use Mod_plcd_NeoHookeanTotalFormulation
   use Mod_Listen
   use Mod_Mesh
   use Mod_plcd
   use Mod_Memor
   !use Mod_plcd_NotDef
   implicit none 
   private
   public AllocatePLCDMaterial

contains

   subroutine AllocatePLCDMaterial(string,Material)
      character(5) :: string
      class(PLCDMaterial), pointer :: Material

      select case (string)
         case ('LINEL') 
            allocate(LinearElastic::Material)
         case ('DAMAG')
            allocate(DamageMaterial:: Material)
         case ('COMPO') 
            allocate(CompositeMaterial::Material)
         case ('SIMP ') 
            allocate(SIMPMaterial::Material)
         case ('SVKIR') 
            allocate(SVKirchhoffMaterial::Material)
         case ('NEOHO')
            allocate(NeoHookeanUpdatedMaterial::Material)
         case ('TOSVK')
            allocate(LinearElastic::Material)
         case ('TONEO')
            allocate(NeoHookeanTotalMaterial::Material)
         case ('NOTDE') 
            !allocate(NotDef::Material)
         case default
            call runend('Wrong type of material')
      end select

   end subroutine
   
  
end module