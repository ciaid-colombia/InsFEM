module def_master
    use typre
    use Mod_timer
    use Mod_GeneralCase
    use Mod_DistributedContainerList
    use Mod_Memor
    implicit none

    type(DistributedContainerList) :: CaseList
    integer(ip)                    :: MPIrank, MPIsize,CasesTotal,MPIcomm
    integer(ip), parameter         :: MPIroot = 0
      
    type(Timer) :: cpu_start(10)   ! CPU for starting operations
    type(Timer) :: cpu_end(10)     ! CPU for ending operations
    type(Timer) :: cpu_total       ! Total CPU time
    integer(ip) :: lun_outpu       ! Output file

end module 
