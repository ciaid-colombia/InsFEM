subroutine case_cputab(a)
   !This routine writes a summary of spent computing time.
   use MPI 
   use Mod_Int2Str
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_Debugging
   use Mod_DistributedContainer
   use Mod_DriverInterface
   use Mod_DC_Driver
   implicit none
   class(GeneralCase), target :: a 
   
   type(masterVariables), pointer :: m => NULL()
   type(domainVariables), pointer :: d => NULL()
   type(adaptiveVariables), pointer :: ad => NULL()
   
   real(rp) :: cpumodut, cputotal
   
   real(rp)    :: cpurefer,cpudenom,cpustart(10) !!
   real(rp)    :: cpustart2(10), cputotal2
   real(rp)    :: cpuend(10),cpuend2(10)
   integer(ip) :: imodu,ifield
   real(rp), allocatable :: cpu_modul2(:,:), auxCpuTims(:),auxCpuTims2(:)
   integer(ip) :: kfl_isSet,auxkfl_isSet,itimer,icpu
   integer :: ierr
   real(rp)    :: Meshcputimes(50),cpuadaptive=0.0_rp, cpumeshprojections=0.0_rp  !!
   
   class(DistributedContainer), pointer :: myDC => NULL()
   class(DriverInterface), pointer :: myDriver => NULL()
   real(rp) :: cpu_modul
   character(DCHashCharSize) :: namod
   
   real(rp) :: cpuminim

   
   m  => a%caseVars%masterVars
   d  => a%caseVars%domainVars
   ad => a%caseVars%adaptiveVars
   

   !Initializations.
   call m%cpu_total%Toc
   call m%cpu_total%GetValue(cputotal)
   do icpu = 1,10
      call m%cpu_start(icpu)%GetValue(cpustart(icpu))
      call m%cpu_end(icpu)%GetValue(cpuend(icpu))
   enddo   

   if (m%MPIrank == m%MPIroot) then
      
      !Total CPU and CPU for starting operations.
      write(m%lun_outpu,100)&
            cputotal,&
            sum(cpustart(2:6)),      100.0_rp* (sum(cpustart(2:6)))/cputotal,&
            cpustart(4),             100.0_rp* cpustart(2)/cputotal,&
            cpustart(5),             100.0_rp* cpustart(2)/cputotal,&
            cpustart(6),             100.0_rp* cpustart(2)/cputotal,&
            cpustart(2),             100.0_rp* cpustart(2)/cputotal,&
            cpustart(3),             100.0_rp* cpustart(3)/cputotal
      
      !Modules time 
      !Get it and write it
      cpumodut = 0.0_rp
      !Loop Drivers and write their times
      call a%DriverList%GetFirst(myDC)
      do while (associated(myDC))
         call ExtractDriver(myDC,myDriver)
         
         call myDC%GetKey(namod)
         call ExtractDriver(myDc,myDriver)
         
         call myDriver%GetTime(cpu_modul)
         cpumodut = cpumodut + cpu_modul
         cpuminim = 1.0d-6*cputotal 
            
         write(m%lun_outpu,101) namod,& 
               cpu_modul, 100.0_rp* cpu_modul/cputotal
         
         !Write the time info internally
         call myDriver%WriteTimes
         
         
         
         call a%DriverList%GetNext(myDC)
      enddo
      
      !Adaptive Information
      if (ad%kfl_AdaptiveRefinement == 1) then
         call ad%cpu_adaptive%GetValue(cpuadaptive)
         write(m%lun_outpu,104) cpuadaptive,100.0_rp*cpuadaptive/cputotal
         call ad%Refiner%WriteTimes(m%lun_outpu,cputotal)
         call d%Mesh%GetTimes(Meshcputimes)
         write(m%lun_outpu,105) Meshcputimes(22),100.0_rp*Meshcputimes(22)/cputotal
      endif
      
      !Mesh Projections
      call d%cpu_meshprojections%GetValue(cpumeshprojections)
      write(m%lun_outpu,107) cpumeshprojections,100.0_rp*cpumeshprojections/cputotal
      
      
      write(m%lun_outpu,103)&
            cpuend(1),  100.0_rp*cpuend(1)/cputotal
      
      !CPU for OUTPUT and OTHER operations.
      cpurefer = cputotal - cpustart(2) -cpustart(4) - cpustart(5) - cpustart(6) &
                 - cpumodut -cpuend(1)-cpuadaptive - cpumeshprojections
      
      write(m%lun_outpu,102)&
            cpurefer,  100.0_rp*cpurefer/cputotal
   endif
   
   
   !If AuxTimers is set, write local
   kfl_isSet = 0
   call m%MasterMemo%alloc(size(AuxTimers),auxCpuTims,'auxCpuTims','cputab')
   call m%MasterMemo%alloc(size(AuxTimers),auxCpuTims2,'auxCpuTims2','cputab')
   do itimer  = 1,size(AuxTimers)
      call AuxTimers(itimer)%isSet(auxkfl_isSet)
      if (auxkfl_isSet == 1) then
         kfl_isSet = kfl_isSet + 1
         call AuxTimers(itimer)%GetValue(auxCpuTims(kfl_isSet))
         if (m%MPIrank == m%MPIroot) write(m%lun_outpu,*) 'AuxTimer_',trim(int2str(itimer)),': ',auxCpuTims(kfl_isSet)
      endif
   enddo
   
   if (m%MPIsize > 1) then
      !AuxTimers, communicate and write
      if (kfl_isSet /= 0) then
         call MPI_REDUCE( auxCpuTims, auxCpuTims2, kfl_isSet, MPI_REAL8, MPI_MAX, m%MPIroot,m%MPIcomm, ierr )
      
         if (m%MPIrank == m%MPIroot) then
            kfl_isSet = 0
            do itimer  = 1,size(AuxTimers)
               call AuxTimers(itimer)%isSet(auxkfl_isSet)
               if (auxkfl_isSet == 1) then
                  kfl_isSet = kfl_isSet + 1
                  write(m%lun_outpu,*) 'AuxTimer_',trim(int2str(itimer)),': ',auxCpuTims2(kfl_isSet)
               endif
            enddo
         endif
      endif
      
   endif
   
   call m%MasterMemo%dealloc(size(AuxTimers),auxCpuTims,'auxCpuTims','cputab')
   call m%MasterMemo%dealloc(size(AuxTimers),auxCpuTims2,'auxCpuTims2','cputab')
   
  !Formats.
   100 format(//,&
            5x,'**************************',/&
            5x,'SUMMARY OF COMPUTING TIMES',/,&
            5x,'**************************',//,&
            10x,'TOTAL CPU TIME              :',F12.2,/,&
            10x,'STARTING OPERATIONS         :',F12.2,' (',F6.2,' % )',/,&
            10x,'     Initialize MPI         :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Reapro                 :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Initialize Par. Lib.   :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Domain Operations      :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     Output Domain          :   ',F12.2,' (',F6.2,' % )')
   
   105 format(&
            10x,'   MESH                     :',F12.2,' (',F6.2,' % )')
   
   104 format(&
            10x,'REFINEMENT OPERATIONS           :',F12.2,' (',F6.2,' % )')
   
   103 format(&
            10x,'TURNOF OPERATIONS           :',F12.2,' (',F6.2,' % )')
            
   107 format(&
            10x,'MESH PROJECTIONS            :',F12.2,' (',F6.2,' % )')         
   
      
   102 format(&
            10x,'OTHER OPERATIONS            :',F12.2,' (',F6.2,' % )')

   131 format(//,//,//,'  MAXIMUM TIME AMONGST ALL PROCESSES',//,//)   
   101 format(&
            10x,a6,' MODULE               :',F12.2,' (',F6.2,' % )')
   
  
   
end subroutine


