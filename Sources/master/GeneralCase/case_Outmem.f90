module Mod_caseOutMem
   use Mod_GeneralCase
   use Mod_CaseVariables
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DC_Driver
   implicit none
   
   type(masterVariables), pointer :: m => NULL()
   type(domainVariables), pointer :: d => NULL()
   
   integer(8) :: ToMem,ToCurrent,ToMax,SolMemo
   real(rp)     :: remem,bymem
   real(rp)     :: rbyte,romem,romax,rSolMemo,rToMem
   character(6) :: lbyte
   

contains        
   subroutine LoopGetMemo(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      integer(8) :: dummy,ModMaxMemo
      
      call ExtractDriver(myDc,myDriver)
      
      call myDriver%GetMemo(dummy,ModMaxMemo,dummy,dummy)
      ToMem = ToMem + ModMaxMemo
   end subroutine
   
   subroutine LoopWriteMemo(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      character(DCHashCharSize) :: namod
      
      integer(8) :: dummy,ModMaxMemo

      
      call ExtractDriver(myDc,myDriver)
      call myDC%GetKey(namod)
      
      call myDriver%GetMemo(dummy,ModMaxMemo,dummy,dummy)
      
      bymem = ModMaxMemo/rbyte
      remem = 100.0_rp*bymem/romem
      write(m%lun_outpu,100) namod,bymem,lbyte,remem
      
      100 format(&
         10x,a6,' MODULE:              ',&
         f7.2,1x,a6,' (',f6.2, ' %)' )
   end subroutine
end module

subroutine case_outmem(a)
   !This routine writes the memory needed.
   use Mod_caseOutMem
   implicit none
   class(GeneralCase), target :: a
   
   m => a%caseVars%masterVars
   d => a%caseVars%domainVars

   if (m%MPIrank == m%MPIroot) then
      !Initialization
      write(m%lun_outpu,10)
      
      !Main memory: domain+master
      call d%Mesh%GetMemo(d%DomCurrentMemo,d%DomMaxMemo,ToCurrent,ToMax)
      ToMem = d%DomMaxMemo
      
      !Add the memory from the modules
      call a%DriverList%LoopList(LoopGetMemo)
      
      !Write total memory
      call Mem2String(Tomem,romem,rbyte,lbyte)
      bymem = d%DomMaxMemo/rbyte
      remem = 100.0_rp*bymem/romem
      write(m%lun_outpu,40) romem,lbyte,bymem,lbyte,remem
      
      !Write Memory depending on the module
      call a%DriverList%LoopList(LoopWriteMemo)

      
         
      !Maximum memory required
      call Mem2String(Tomax,romax,rbyte,lbyte)
      write(m%lun_outpu,60) romax,lbyte

      !Parallel Library maximum memory
      call m%ParallelLibrary%GetMaximumMemory(SolMemo)
      call Mem2String(SolMemo,rSolMemo,rbyte,lbyte)
      write(m%lun_outpu,80) rSolMemo, lbyte
   endif
   

   !Formats. 

10 format(//,&
        5x,'**************************',/&
        5x,'SUMMARY OF REQUIRED MEMORY',/,&
        5x,'**************************',//)
40 format(&
        5x,'     MEMORY REQUIRED:',//,&
        10x,    'TOTAL MEMORY:               ',&
        f7.2,1x,a6,/,  &
        10x,    'CONSTRUCTION OF THE DOMAIN: ',&
        f7.2,1x,a6,' (',f6.2, ' %)' )

50 format(&
        10x,'SOLVERS:                    ',&
        f7.2,1x,a6,' (',f6.2, ' %)' )
  
60 format(/,&
        10x,    'MAXIMUM MEMORY REQUIRED:    ',&
        f7.2,1x,a6,//)

70 format(&
        5x,'     WARNING, MEMORY LEAK:       ',f7.2,1x,a6,/,&
        5x,'     RUN MEMCHECK')
        
80 format(/,&
        10x,    'MAXIMUM MEMORY PAR. LIB.:   ',&
        f7.2,1x,a6,//)       



           
end subroutine 
