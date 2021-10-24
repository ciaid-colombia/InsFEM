subroutine case_Endste(a)
   use typre
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_DomainVariables
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_DCHashCharSize
   use Mod_DriverInterface
   implicit none
   class(GeneralCase), target :: a
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
         
   
   type(masterVariables), pointer :: m => NULL()
   type(caseVariables), pointer :: c => NULL()
   type(domainVariables), pointer :: d => NULL()
   
   integer(ip), allocatable :: processor(:)
   integer(ip) :: npoin,npoinLocal,ipoin
   
   interface
      subroutine case_PostprocessProcessor(a)
         use typre
         use Mod_GeneralCase
         implicit none
         class(GeneralCase), target :: a
      end subroutine
   
   
   end interface
   
   c => a%caseVars
   m => c%masterVars
   d => c%domainVars

   !Initializations.
   m%kfl_gotim = 0

   if (m%MPIrank == m%MPIroot) write(m%lun_outpu,10) 
   
   !Loop through drivers and do endste
   call a%DriverList%GetFirst(myDC)
   do while(associated(myDC))
   
      call ExtractDriver(myDC,myDriver)
      call myDriver%Endste(c)
      
      call a%DriverList%GetNext(myDC)
   enddo

   call case_PostprocessProcessor(a)
     
   
      ! Check if the time evolution has to be stopped or not.
   if( (m%timef-m%ctime)<=epsilon(1.0_rp) )  m%kfl_gotim = 0
   if(m%istep>=m%nsmax)              m%kfl_gotim = 0


   10 format(/,&
         10x,'Module     Total      Turnon     Getste     Begste     Doiter    LinBuild   LinSolve    Endite     Endste    Output     Turnof     Restar     Refine     Projec ')
  
end subroutine 

subroutine case_PostprocessProcessor(a)
   use typre
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_DomainVariables
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_DCHashCharSize
   use Mod_DriverInterface
   implicit none
   class(GeneralCase), target :: a
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
         
   
   type(masterVariables), pointer :: m => NULL()
   type(caseVariables), pointer :: c => NULL()
   type(domainVariables), pointer :: d => NULL()
   
   integer(ip), allocatable :: processor(:)
   integer(ip) :: npoin,npoinLocal,ipoin
   
   c => a%caseVars
   m => c%masterVars
   d => c%domainVars

   !Postprocess processor if required
   if ((m%PostprocessProcessorEvery > 0) .and. (m%istep >= m%StartPostprocessAt)) then
      if (mod(m%istep,m%PostprocessProcessorEvery) == 0) then      
         call d%Mesh%GetNpoin(npoin)
         call d%Mesh%GetNpoinLocal(npoinLocal)
         call m%MasterMemo%alloc(npoin,processor,'processor','domain')
         do ipoin = 1,npoinLocal
            processor(ipoin) = m%MPIrank
         enddo
         call d%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,processor)
         call m%FilePostpr%postpr(processor,'proc',m%istep,m%ctime,d%Mesh)
         call m%MasterMemo%dealloc(npoin,processor,'processor','domain')
      endif    
   endif    

end subroutine


subroutine case_Output(a)
    use typre
    use Mod_GeneralCase
    use Mod_caseVariables
    use Mod_DistributedContainer
    use Mod_DC_Driver
    use Mod_DCHashCharSize
    use Mod_DriverInterface
    implicit none
    class(GeneralCase), target :: a

    class(DistributedContainer), pointer :: myDC => NULL()
    class(DriverInterface), pointer :: myDriver => NULL()


    type(masterVariables), pointer :: m => NULL()
    type(caseVariables), pointer :: c => NULL()

    c => a%caseVars
    m => a%caseVars%masterVars

    !Loop through drivers and print output
    call a%DriverList%GetFirst(myDC)
    do while(associated(myDC))

        call ExtractDriver(myDC,myDriver)
        call myDriver%Output(c)

        call a%DriverList%GetNext(myDC)
    enddo

end subroutine
