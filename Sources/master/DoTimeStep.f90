subroutine DoTimeStep(kfl_gotim)
   use def_master
   use Mod_DC_GeneralCase
   use Mod_caseVariables
   use Mod_DistributedContainer
   implicit none 
   integer(ip) :: kfl_gotim
   
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()
   type(masterVariables), pointer :: m => NULL()
   integer(ip) :: auxkflgotim = 0, cpiter,aux_ncases,maxiters,auxmaxiters 
   logical :: kfl_coupledconv = .false.
   logical :: kfl_iterconverged= .false.
   logical :: auxkfl_coupledconv = .false.
   character(DCHashCharSize) :: nick
   

   !Do the time step for all cases

   kfl_gotim       = 0
   kfl_coupledconv = .false.
   maxiters        = -1
   auxmaxiters     = -1
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      !Prepare time step

      call myCase%Adaptive(2)

      !Update Channels
      call myCase%ReSetupInterpolators  !Only if interpolated and necessary
      call myCase%UpdateChannels
           
      call myCase%Setste
      call myCase%Begste
      
      !check if we need to converge coupled containers
      call myCase%GetDoCoupledConv(auxkfl_coupledconv)
      if(auxkfl_coupledconv) kfl_coupledconv = .true.

      !check for max coupled iterations 
      auxmaxiters = maxiters
      call myCase%GetMaxCoupledIters(maxiters)
      maxiters = max(maxiters,auxmaxiters)

      call CaseList%GetNext(myDC)
   enddo

   cpiter = 1
   kfl_iterconverged = .false.
   case_coupling: do while (kfl_iterconverged .eqv. .false.)

      aux_ncases = 0
      kfl_iterconverged = .true.

      call CaseList%GetFirst(myDC)
      do while (associated(myDC))
         call ExtractGeneralCase(myDC,myCase)
         m => myCase%caseVars%masterVars
         m%cpiter = cpiter
         !Do the time step
         call myCase%UpdateChannels
         call myCase%DoTimeStep

         if(kfl_coupledconv) then
            call myCase%GlobalCouplingConvergence(cpiter,auxkfl_coupledconv)
            if (auxkfl_coupledconv .eqv. .true. ) aux_ncases = aux_ncases + 1
         end if

         auxkfl_coupledconv = .false.
         call CaseList%GetNext(myDC)
      enddo

      if (aux_ncases .ne. CasesTotal .and. kfl_coupledconv) then 
         kfl_iterconverged = .false.
         if (cpiter == maxiters) kfl_iterconverged = .true.
      endif

      cpiter = cpiter + 1

   end do case_coupling

   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)

      !End time step
      call myCase%Endste

      !Case output
      call myCase%Output

      !Ask if the case needs to continue
      call myCase%GetGoTimeStepping(auxkflgotim)
      if (auxkflgotim == 1) kfl_gotim = 1

      call CaseList%GetNext(myDC)
   enddo

end subroutine
