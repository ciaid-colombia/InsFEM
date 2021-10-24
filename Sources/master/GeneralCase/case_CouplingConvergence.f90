subroutine case_LocalCouplingConvergence(a)
   !This routine checks the general convergence of the run.
   use MPI
   use typre
   use Mod_Timer
   use Mod_DCHashCharSize
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_caseVariables
   use Mod_DriverInterface
   use Mod_GeneralCase
   implicit none
   class(GeneralCase), target :: a
   
   real(rp)          :: cpu_refer
   integer(ip), save :: itert=0
   integer           :: ierr
   
   type(caseVariables), pointer :: c => NULL()
   type(masterVariables), pointer :: m => NULL()
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
   character(DCHashCharSize) :: namod
   
   integer(ip) :: nmodu, imodu
   real(rp) :: glres
   
   !Saving for the internal Loop subroutines
   !save c,m
   
   c => a%caseVars
   m => a%caseVars%masterVars
   glres = 0.0_rp

   if(m%MPIrank == m%MPIroot) then
      !Initializations.
      m%kfl_gocou = 0
      call cputim(cpu_refer)
      cpu_refer = cpu_refer - m%cpu_initi
      write(m%lun_outpu,10) m%iiter ,cpu_refer

      !Write to convergence file and keep iterating...
      itert=itert+1
      if(itert==1) then
         write(m%lun_conve,11)
         write(m%lun_conve,12, advance='no')
         
         !Get drivers and write their names
         call a%DriverList%GetFirst(myDC)
         do while (associated(myDC))
            call myDC%GetKey(namod)
            write(m%lun_conve,13,advance='no') adjustl(trim(namod))
            call a%DriverList%GetNext(myDC)
         enddo
         write(m%lun_conve,30,advance='no')
      end if
      
      write(m%lun_conve,14,advance='no') m%istep,m%iiter,m%ctime
      !Loop through drivers and write their residuals
      call a%DriverList%GetFirst(myDC)
      do while (associated(myDC))
         call ExtractDriver(myDC,myDriver)
      
         call myDriver%Convergence(c,glres)
         if (glres > m%tolgl) then
            m%kfl_gocou = 1
         endif
         write(m%lun_conve,15,advance='no') glres
         
         call a%DriverList%GetNext(myDC)
      enddo
      write(m%lun_conve,30,advance='no')
      
      if(m%iiter>=m%mitgl) m%kfl_gocou = 0                                !False convergence
      if (m%kfl_flush == 1) then
         call flush(m%lun_outpu)
         call flush(m%lun_conve)
      endif
   endif
   
   CALL MPI_BCAST(m%kfl_gocou, 1, MPI_INTEGER4, m%MPIroot, m%MPIcomm, ierr)
   m%iiter = m%iiter + 1   
   
   !Formats.
   10 format(/,&
            10x,'GLOBAL ITERATION NUMBER: ',i5,/,10x, 'CURRENT CPU TIME:        ',e12.6)

   11 format('$ ','       Time','     Global','       Current','  Module residuals -->')
   12 format('$ ','       step','  iteration','          time')
   14 format(4x,i9,2x,i9,20(2x,e12.6))
   13 format(20('   ',a5))   
   20 format(&
            1x,'Current step',2x,'Maximum steps',/,&
            1x,i7,6x,i7)
   15 format(20(2x,e12.6))
   30 format(/)  
            
end subroutine 

subroutine case_GlobalCouplingConvergence(a,cpiter,kfl_coupledConv)
   !This routine checks the general convergence of the run.
   use typre
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DCHashCharSize
   use MPI
   use Mod_Timer
   use Mod_DC_Driver
   implicit none
   class(GeneralCase), target :: a
   
   type(caseVariables), pointer :: c => NULL()
   type(masterVariables), pointer :: m => NULL()
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
   logical     :: kfl_coupledConv,auxkfl_coupledConv
   integer(ip) :: ierr,aux_ndrivers,cpiter
   real(rp)    :: cpres,auxcpres,beta,maxcpiters

   c => a%caseVars
   m => a%caseVars%masterVars
   maxcpiters = m%kfl_maxCoupledIters
   beta = m%couplingBeta
   maxcpiters = m%kfl_maxCoupledIters
   auxcpres = 10
   cpres    = 1

   if(m%MPIrank == m%MPIroot) then
       aux_ndrivers = 0
       kfl_coupledConv = .false.
       call a%DriverList%GetFirst(myDC)
       do while (associated(myDC))
           call ExtractDriver(myDC,myDriver)

           call myDriver%CoupledConvergence(c,auxkfl_coupledConv,cpres)
           call myDriver%GetCoupledResidual(auxcpres)
           call myDriver%SetCoupledResidual(cpres)
           if ((auxkfl_coupledConv .eqv. .true.) .or. ((cpiter .ge. 2) .and. (auxcpres/cpres .lt. beta)) .or. (cpiter .gt. maxcpiters)) then 
               aux_ndrivers = aux_ndrivers + 1
           endif

           call a%DriverList%GetNext(myDC)
       enddo
       if (aux_ndrivers == c%ndrivers) kfl_coupledConv = .true.
   endif

   call MPI_BCAST(kfl_coupledConv, 1, MPI_LOGICAL, m%MPIroot, m%MPIcomm, ierr)

end subroutine 
