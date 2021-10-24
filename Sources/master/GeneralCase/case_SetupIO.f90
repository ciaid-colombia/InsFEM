subroutine case_SetupIO(a,itask)
   use typre
   use Mod_GeneralCase
   use Mod_CaseVariables
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
   type(adaptiveVariables), pointer :: ad => NULL()
   
   integer(ip), allocatable :: processor(:)
   integer(ip) :: npoin,npoinLocal,ipoin,itask
   
   real(rp), pointer :: coord(:,:) => NULL()
   integer(ip) :: kfl_Init,kfl_rstar,kfl_inter,kfl_preli
   integer(ip) :: auxkfl_rstar,auxkfl_inter,auxkfl_preli

   interface
      subroutine case_SetupReader(a,itask)
         use typre
         use Mod_GeneralCase
         use Mod_caseVariables
         class(GeneralCase), target :: a
         integer(ip) :: itask
      end subroutine
      subroutine case_SetupWriter(a)
         use typre
         use Mod_GeneralCase
         use Mod_caseVariables
         class(GeneralCase), target :: a
         integer(ip) :: itask
      end subroutine
   end interface
   
   c  => a%caseVars
   m  => c%masterVars
   d  => c%domainVars
   ad => c%adaptiveVars

   kfl_rstar = 0
   kfl_inter = 0
   kfl_preli = 0

   call a%DriverList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractDriver(myDC,myDriver)
      call myDriver%GetRestartFlags(auxkfl_rstar,auxkfl_inter,auxkfl_preli)
      if (auxkfl_rstar == 1) kfl_rstar = 1
      if (auxkfl_inter == 1) kfl_inter = 1
      if (auxkfl_preli == 1) kfl_preli = 1
      call a%DriverList%GetNext(myDC)
   enddo

   if (kfl_inter == 1 .and. ad%UpdateFreq .ge. 0) ad%UpdateIO = .true.

   select case(itask)
   case(0)
      if (kfl_rstar==1) then 
         if (kfl_inter==1) then 
            call d%Mesh%GetCoord(coord)
            call a%OldDomain(1)
            call c%masterVars%Int_Restart%SetTol(c%masterVars%Interp_Tol)
            call c%masterVars%Int_Restart%IsInitialized(kfl_Init)
            if (kfl_Init == -1) then
               call c%masterVars%Int_Restart%Initialize_Interp(d%OldMesh,coord)
            endif
         endif
         call case_SetupReader(a,kfl_inter)
      endif
      if (kfl_preli==1) call case_SetupWriter(a)
   case(1)
      if (kfl_rstar==1) call m%Readerpr%DeallocateReader
      if (kfl_inter==1) then
         call c%masterVars%Int_Restart%IsInitialized(kfl_Init)
         if (kfl_Init /= -1) then
            call c%masterVars%Int_Restart%Finalize    
         endif
         call a%OldDomain(2)
         call a%OldDomain(3)  
      endif
      if (kfl_preli==1) call m%Writerpr%DeallocateWriter
   end select
     
end subroutine 

subroutine case_SetupWriter(a)
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
   
   integer(ip)             :: npoinLocal,gnpoin,ipoin
   integer(ip),allocatable :: Indexes(:)
   
   c => a%caseVars
   m => c%masterVars
   d => c%domainVars
      
   call d%Mesh%GetNpoinLocal(npoinLocal)
   call d%Mesh%GetGnpoin(gnpoin)

   call m%MasterMemo%alloc(npoinLocal,Indexes,'Indexes','php_restar')
   do ipoin = 1,npoinLocal
      Indexes(ipoin) = ipoin
   enddo

   call d%Mesh%Local2Global(npoinLocal,Indexes,Indexes)
   call d%Mesh%Global2Initial(npoinLocal,Indexes,Indexes)
   call m%Writerpr%SetWriter(npoinLocal,gnpoin,Indexes,m%MasterMemo)

   call m%MasterMemo%dealloc(npoinLocal,Indexes,'Indexes','php_restar')

end subroutine

subroutine case_SetupReader(a,kfl_inter)
   use typre
   use MPI
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
   
   integer(ip)             :: inipoin,endpoin,kfl_inter
   integer(ip)             :: npoinLocal,gnpoin,ipoin,ierr
   integer(ip),allocatable :: Indexes(:),Indexes2(:),ProcNumber(:),rankPoinLocal(:)
   
   c => a%caseVars
   m => c%masterVars
   d => c%domainVars

   if (kfl_inter == 1) then
      call d%OldMesh%GetNpoinLocal(npoinLocal)
      call d%OldMesh%GetGnpoin(gnpoin)
   else
      call d%Mesh%GetNpoinLocal(npoinLocal)
      call d%Mesh%GetGnpoin(gnpoin)
   endif

   call m%MasterMemo%alloc(m%MPIsize,rankPoinLocal,'rankPoinLocal','php_restar')
   rankPoinLocal = 0
   rankPoinLocal(m%MPIrank+1) = npoinLocal
   call MPI_ALLREDUCE(MPI_IN_PLACE,rankPoinLocal,m%MPIsize,MPI_INTEGER4,MPI_MAX,m%MPIcomm,ierr)

   call m%MasterMemo%alloc(npoinLocal,Indexes,'Indexes','php_restar')
   call m%MasterMemo%alloc(npoinLocal,Indexes2,'Indexes2','php_restar')
   call m%MasterMemo%alloc(npoinLocal,ProcNumber,'ProcNumber','php_restar')  

   inipoin = sum(rankPoinLocal(1:m%MPIrank))+1
   endpoin = sum(rankPoinLocal(1:m%MPIrank+1))
   call m%MasterMemo%dealloc(m%MPIsize,rankPoinLocal,'rankPoinLocal','php_restar')
   Indexes = Indexes + inipoin-1

   do ipoin = 1,npoinLocal
      Indexes(ipoin) = ipoin
   enddo
   Indexes = Indexes + inipoin-1

   if (kfl_inter == 1) then
      call d%OldMesh%Initial2ProcNumber(npoinLocal,Indexes,ProcNumber)
      call d%OldMesh%Initial2Global(npoinLocal,Indexes,Indexes)
      call m%Readerpr%SetReader(npoinLocal,gnpoin,Indexes2,Indexes,ProcNumber,m%MasterMemo)
      call d%OldMesh%Global2Local(npoinLocal,Indexes2,Indexes2)
      call m%Readerpr%SetG2L(Indexes2)
   else    
      call d%Mesh%Initial2ProcNumber(npoinLocal,Indexes,ProcNumber)
      call d%Mesh%Initial2Global(npoinLocal,Indexes,Indexes)
      call m%Readerpr%SetReader(npoinLocal,gnpoin,Indexes2,Indexes,ProcNumber,m%MasterMemo)
      call d%Mesh%Global2Local(npoinLocal,Indexes2,Indexes2)
      call m%Readerpr%SetG2L(Indexes2)
   endif

   call m%MasterMemo%dealloc(npoinLocal,ProcNumber,'ProcNumber','php_restar')
   call m%MasterMemo%dealloc(npoinLocal,Indexes2,'Indexes2','php_restar')
   call m%MasterMemo%dealloc(npoinLocal,Indexes,'Indexes','php_restar')

end subroutine
