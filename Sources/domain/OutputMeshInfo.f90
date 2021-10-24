subroutine OutputMeshInfo(a)
   use MPI
   use typre
   use def_parame
   use Mod_Iofile
   use Mod_Int2Str
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip)    :: optpoin,refnelem
   character(150) :: fil_outpu, fil_solve
   real(rp)       :: cputimes(25),cputimes2(25)
   
   integer :: ierr
   
   interface
      subroutine OutputMeshTimes(a,cputimes)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
         real(rp) :: cputimes(25)
      end subroutine
   end interface
   
   if (a%MPIrank == a%MPIroot) then
      fil_outpu = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//adjustl(trim(int2str(a%MPIrank)))//'.dom'
      call iofile(zero,a%lun_outpu_dom,fil_outpu,'MESH OUTPUT')
      
      write(a%lun_outpu_dom,*) '     DOMAIN DECOMPOSITION RESULTS'
      write(a%lun_outpu_dom,*) '     ------------------------------------------'
      write(a%lun_outpu_dom,100) a%gnpoin
      optpoin = a%gnpoin/a%MPIsize
      write(a%lun_outpu_dom,101) optpoin
      write(a%lun_outpu_dom,102) a%MaxNpoinPerProc
      write(a%lun_outpu_dom,103) real(a%MaxNpoinPerProc-optpoin)/real(optpoin)*100.0

      write(a%lun_outpu_dom,104) a%gnpoinGhost
      write(a%lun_outpu_dom,105) a%MaxNpoinGhostPerProc
      
      
      write(a%lun_outpu_dom,200) a%gnelem
      refnelem = a%gnelem/a%MPIsize
      write(a%lun_outpu_dom,201) refnelem
      write(a%lun_outpu_dom,202) a%MaxNelemPerProc
      write(a%lun_outpu_dom,203) real(a%MaxNelemPerProc-refnelem)/real(refnelem)*100.0
      
   endif

      
   call a%GetTimes(cputimes)
   
   if (a%MPIrank == a%MPIroot) call OutputMeshTimes(a,cputimes)
   
   if (a%MPIsize > 1) then
      call MPI_REDUCE( cputimes, cputimes2, size(cputimes), MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
      
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_outpu_dom,131) 
         call OutputMeshTimes(a,cputimes2)
      endif
   endif
   
   if (a%MPIrank == a%MPIroot) then
      call flush(a%lun_outpu_dom)
   endif
   
   
   
   100 format(10x,'GLOBAL NUMBER OF POINTS:                    ',i10)
   101 format(10x,'OPTIMUM NUMBER OF POINTS PER PROCESSOR:     ',i10)
   102 format(10x,'MAXIMUM NUMBER OF POINTS PER PROCESSOR:     ',i10)
   103 format(10x,'DIFFERENCE:                                 ',f10.2,' %'/)
   
   104 format(10x,'TOTAL NUMBER OF GHOST POINTS:               ',i10)
   105 format(10x,'MAXIMUM NUMBER OF GHOST POINTS PER PROC     ', i10,/)
   
   200 format(10x,'GLOBAL NUMBER OF ELEMENTS:                  ',i10)
   201 format(10x,'REFERENCE NUMBER OF ELEMENTS PER PROCESSOR: ',i10)
   202 format(10x,'MAXIMUM NUMBER OF ELEMENTS PER PROCESSOR:   ',i10)
   203 format(10x,'DIFFERENCE:                                 ',f10.2,' %'/)
   
   131 format(//,//,//,'  MAXIMUM TIME AMONGST ALL PROCESSES',//,//)   
   
   
end subroutine

subroutine OutputMeshTimes(a,cputime)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   real(rp) :: cputime(25)
   
   integer(ip) :: itime
   
    write(a%lun_outpu_dom,300) cputime(1), (cputime(itime),cputime(itime)/cputime(1)*100_rp,itime=2,16)
   if (a%kfl_TestCommunications == 1) then
       write(a%lun_outpu_dom,301) (cputime(itime), cputime(itime)/cputime(1)*100_rp,itime=17,21), a%Bandwidth,a%Bandwidth/a%MPIsize
   endif
   
   300 format( '     SUMMARY OF COMPUTING TIMES',/,&
               '     ------------------------------------------',//,&
               10x,'TOTAL CPU TIME:            ',f10.2,/,&
               10x,'     ReaLnods                        :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     Partitionate                    :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     LocalDimensionsAndIRenumbering  :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     GhostPoints                     :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     BuildLocalOrdering              :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     ReaCoord                        :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     LnodsToLocal                    :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     LocalGraph                      :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     BuildArrayCommunicator          :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     GhostCoord                      :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     PeriodicBoundaryConditions      :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     Dombou                          :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     ExtnorLpoty                     :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     ComputeVmass                    :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'     Blbopo                          :   ',F12.2,' (',F6.2,' % )')
   301 format( 10x,'     Test Communications: ',/,&
               10x,'        100 GhostCommunicates (Coord):   ',F12.2,' (',F6.2,' % )',/,&
               10x,'        100 AllToAll                 :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'        100 Reduce                   :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'        100 Broadcast                :   ',F12.2,' (',F6.2,' % )',/,&
               10x,'        Test Bandwidth               :   ',F12.2,' (',F6.2,' % )',//,&
               10x,'        Estimated Bandwidth (Gb/s)   :   ',F12.2, /,&
               10x,'        Estimated Bandwidth/Proc (Gb/s)   :   ',e12.4)
               

end subroutine

