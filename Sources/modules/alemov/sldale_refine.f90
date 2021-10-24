subroutine sldale_refine(a,itask)
   use typre
   use Mod_phpRefineArrays
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem) :: a
   character(6) :: itask
   integer(ip) :: npoin
   
   call a%Mesh%GetNpoin(npoin)

   call php_RefineArrays(a,itask,a%bdisp,'bdisp') 
   call php_RefineArrays(a,itask,a%bres,'bres') 

   a%solid%kfl_adap = a%kfl_adap
   a%solid%Refiner => a%Refiner
   call a%solid%Refine(itask)

   a%solid%LinearSystem => a%LinearSystem
   a%solid%kfl_SkipLinearSystemRefinement = 1

   call a%Memor%realloc(a%solid%ndofn,npoin,a%solid%unkno,'unkno',trim(a%exmod)//'_LinearSystemMemall')
   a%solid%unkno = 0.0_rp

end subroutine
