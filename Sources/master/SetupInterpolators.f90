subroutine SetupInterpolators
   use Mod_GeneralCase
   use Mod_DC_GeneralCase
   use Mod_DistributedContainer
   use Mod_DCHashCharSize
   use def_master
   implicit none
   
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()
   character(DCHashCharSize) :: myCaseKey
   
   !Loop through all cases so that they setup the interpolators
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      
      call ExtractGeneralCase(myDC,myCase)
      call myDC%GetKey(myCaseKey)
      
      call myCase%SetupInterpolators(myCaseKey)
      
      call CaseList%GetNext(myDC)
   enddo
end subroutine

subroutine FinalizeInterpolators
   use Mod_GeneralCase
   use Mod_DC_GeneralCase
   use Mod_DistributedContainer
   use Mod_DCHashCharSize
   use def_master
   implicit none
   
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()
   character(DCHashCharSize) :: myCaseKey
   
   !Loop through all cases so that they setup the interpolators
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      
      call ExtractGeneralCase(myDC,myCase)
      call myDC%GetKey(myCaseKey)
      
      call myCase%FinalizeInterpolators
      
      call CaseList%GetNext(myDC)
   enddo
end subroutine