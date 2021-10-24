subroutine Turnof
   use def_master
   use Mod_DC_GeneralCase
   use Mod_DistributedContainer
   use def_parame
   use Mod_iofile
   use Mod_Memor
   implicit none 
    
   class(DistributedContainer), pointer :: myDC => NULL()
   type(GeneralCase), pointer :: myCase => NULL()
   
   if(MPIrank==MPIroot) write(*,*) 'Case: solving ... Done'
   
   !Turnof everything
   if(MPIrank==MPIroot) write(*,*) 'Case: Turning off'

   !Normal turnof procedures
   call CaseList%GetFirst(myDC)
   do while (associated(myDC))
      call ExtractGeneralCase(myDC,myCase)
      
      call myCase%Turnof
      !Deallocate it
      deallocate(myCase)
      
      call CaseList%GetNext(myDC)
   enddo
    
   !Free the list of cases
   call CaseList%Finalize
   
   !Write the times
   call Cputab
   
   !Stop the run.
   call runend('O.K.!')
   
   call FinalizeMPI
   
    
contains
   
   
   subroutine Cputab
   
      real(rp) :: cputotal
      real(rp)    :: cpustart(10) !!
      real(rp)    :: cpuend(10)
      
      integer(8) :: ToMem,ToCurrent,ToMax,SolMemo,CurrentMemo,maxmemo
      real(rp)     :: remem,bymem
      real(rp)     :: rbyte,romem,romax,rSolMemo,rToMem
      character(6) :: lbyte
   
      integer(ip) :: icpu
      
      type(MemoryMan) :: Memor
   
      call cpu_total%Toc
      
      call cpu_total%GetValue(cputotal)
      do icpu = 1,10
         call cpu_start(icpu)%GetValue(cpustart(icpu))
         call cpu_end(icpu)%GetValue(cpuend(icpu))
      enddo   
      
      if (MPIrank == MPIroot) then
         write(lun_outpu,100) cputotal,&
               sum(cpustart(2:6)),      100.0_rp* (sum(cpustart(2:6)))/cputotal,&
               cpustart(4),             100.0_rp* cpustart(2)/cputotal,&
               cpustart(5),             100.0_rp* cpustart(2)/cputotal
         
         write(lun_outpu,*)
         write(lun_outpu,*) 'Here the timing for the cases needs to be added'    
         write(lun_outpu,*) 'Also the timing for the MatchCases'
         
         call Memor%GetValue(CurrentMemo,MaxMemo,ToCurrent,ToMax)
         
         !Maximum memory required
         call Mem2String(Tomax,romax,rbyte,lbyte)
         write(lun_outpu,60) romax,lbyte
         if (ToCurrent /= 0_8) then
            write(lun_outpu,70) real(ToCurrent)/rbyte, lbyte
         endif
         
         !Close the file
         call iofile(two,lun_outpu,' ','CloseGlobalOutputFile','old')
      endif
      
      
      !Formats.
      100 format(//,&
            5x,'**************************',/&
            5x,'SUMMARY OF COMPUTING TIMES',/,&
            5x,'**************************',//,&
            10x,'TOTAL CPU TIME              :',F12.2,/,&
            10x,'STARTING OPERATIONS         :',F12.2,' (',F6.2,' % )',/,&
            10x,'     Initialize MPI         :   ',F12.2,' (',F6.2,' % )',/,&
            10x,'     GetArgsAndEnv          :   ',F12.2,' (',F6.2,' % )')
   
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

    

end subroutine Turnof
