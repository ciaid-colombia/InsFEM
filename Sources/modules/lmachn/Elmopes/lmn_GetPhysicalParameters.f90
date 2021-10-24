subroutine lmn_GetPhysicalParameters(a,acden,acvis,acsph,actco,actex,acgas,acpth,acrea,acsou)
   use typre
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   real(rp), optional    :: acden,acvis,acsph,actco,actex,acgas,acrea,acsou 
   real(rp), optional    :: acpth(*) 
   
   if (present(acden)) acden = a%densi
   if (present(acvis)) acvis = a%visco
   if (present(acsph)) acsph = a%cphea
   if (present(actco)) actco = a%tcond
   if (present(actex)) actex = a%texpc
   if (present(acgas)) acgas = a%sgasc
   if (present(acpth)) acpth(1:a%ncomp) = a%pther
   if (present(acrea)) acrea = a%react
   if (present(acsou)) acsou = a%sourc
end subroutine
   
