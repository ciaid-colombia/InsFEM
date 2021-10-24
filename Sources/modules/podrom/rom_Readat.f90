subroutine rom_Readat(a)
   use typre
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   
   character(150) :: outstr

   a%kfl_Precondition = 0
   a%kfl_SubstractMean  = .false.

   a%kfl_massMatrix = .true.
   a%kfl_outBasis   = .false.

   a%cotol = 1e-6
   a%basis_energy = 1
   a%basis_number = 10
   a%kfl_basis_number = .false.
   a%numBasisFiles = 1
   a%nNonLinearIterations = 10
   a%kfl_eigentype = 'SVD'
   
   a%kfl_SaveSnap          = .false.
   a%kfl_BkpOpen           = .false.
  
   a%RefinerErrorEstimator = 'ZZ'
   a%RefinerErrorCriteria  = 'Tolerance'
   a%RefinerErrorLimits = (/ 1e-3, 1e-4 /)

   a%param = 0.0

   !Set the units for reading
   call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
   
   !Output string
   outstr = 'ROM_READAT'
   
   !Reach the section
   call a%Listener%listen(outstr)
   do while(a%Listener%words(1)/='REDUC')
      call a%Listener%listen(outstr)
   end do
   
   !Begin to read data
   do while(a%Listener%words(1)/='ENDRE')
      call a%Listener%listen(outstr)

      if(a%Listener%words(1) == 'STAGE') then
         if (a%Listener%exists('FOM  ')) then 
            a%kfl_itask = 'FOM'
            if (a%Listener%exists('SNAPS')) then
               a%kfl_SaveSnap = .true.
               if (a%Listener%exists('OLD  ')) a%kfl_BkpOpen = .true.
            end if
         elseif (a%Listener%exists('BASIS')) then
            if (a%Listener%exists('SNAPS')) then
               a%kfl_itask = 'BFS'
            elseif (a%Listener%exists('BASIS')) then
               a%kfl_itask = 'BAS'
               a%numBasisFiles = a%Listener%param(3)
            end if
         elseif (a%Listener%exists('SNAPS')) then
            a%kfl_itask = 'SNP'
            if (a%Listener%exists('OLD  ')) a%kfl_BkpOpen = .true.
         elseif (a%Listener%exists('INTER')) then
            a%kfl_itask = 'INT'
         elseif (a%Listener%exists('ROM  ')) then
            a%kfl_itask = 'ROM'
         endif
      elseif(a%Listener%words(1)=='STYPE') then
         a%kfl_eigentype = a%Listener%words(2)
      elseif (a%Listener%words(1) == 'MASSP') then
         if (a%Listener%exists('YES  ') .or. a%Listener%exists('ON   ')) a%kfl_massMatrix = .true.
         if (a%Listener%exists('NO   ') .or. a%Listener%exists('OFF  ')) a%kfl_massMatrix = .false.
      elseif(a%Listener%words(1) == 'PROJE') then
         if (a%Listener%exists('GALER')) then
            a%kfl_Precondition = 0
         endif
         if (a%Listener%exists('PETRO')) then
            a%kfl_Precondition = 1
         endif
      elseif(a%Listener%words(1)=='NSNAP') then
         a%nsnap = a%Listener%param(1)
      elseif(a%Listener%words(1)=='INTER') then
         a%SnapshotsInterval = a%Listener%param(1)
      elseif(a%Listener%words(1) == 'MEANS') then
         if (a%Listener%exists('YES  ') .or. a%Listener%exists('ON   ')) then 
            a%kfl_SubstractMean = .true.
         elseif (a%Listener%exists('NO   ') .or. a%Listener%exists('OFF  ')) then
            a%kfl_SubstractMean = .false.
         endif 
      elseif (a%Listener%words(1) == 'NONLI') then
         a%nNonLinearIterations = a%Listener%getint('NONLI',1_ip,'# of Non-linear iterations')
      elseif (a%Listener%words(1) == 'CONVE') then
         a%cotol = a%Listener%param(1)
      elseif (a%Listener%words(1) == 'ENERG') then
         a%basis_energy = a%Listener%param(1)
      elseif (a%Listener%words(1) == 'NUMBE') then
         a%basis_number = a%Listener%param(1)
         a%kfl_basis_number = .true.
      else if(a%Listener%words(1)=='ERROR') then
          a%RefinerErrorEstimator = trim(a%Listener%words(2))
          if (a%Listener%exists('TOLER')) then
              a%RefinerErrorCriteria = 'Tolerance'
          elseif (a%Listener%exists('ELEME')) then
              a%RefinerErrorCriteria = 'Elements'
          endif
          a%RefinerErrorLimits(1) = a%Listener%getrea('MAXER',1e-3_rp,'Error Estimator Maximum')
          a%RefinerErrorLimits(2) = a%Listener%getrea('MINER',1e-4_rp,'Error Estimator Minimum')
      elseif (a%Listener%words(1) == 'OUTBA') then
         if (a%Listener%exists('YES  ') .or. a%Listener%exists('ON   ')) a%kfl_outBasis = .true.
         if (a%Listener%exists('NO   ') .or. a%Listener%exists('OFF  ')) a%kfl_outBasis = .false.
      elseif (a%Listener%words(1) == 'PARAM') then
         a%param = a%Listener%param(1)
      endif  
   enddo
   if (a%kfl_itask == "BAS") a%kfl_SubstractMean = .false.
   if (a%kfl_itask == "INT") a%basis_energy = 1

end subroutine

subroutine rom_ReadatMPI(a)
   use MPI
   use typre
   use Mod_BroadCastBuffer
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   type(BroadCastBuffer) :: BBuffer
   
   call BBuffer%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call BBuffer%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call BBuffer%Initialize(100,500)
   
   call BBuffer%add(3,a%kfl_eigentype)
   call BBuffer%add(3,a%kfl_itask)
   call BBuffer%add(a%kfl_SaveSnap)
   call BBuffer%add(a%kfl_BkpOpen)
   call BBuffer%add(a%kfl_snapshotSpecific)
   
   call BBuffer%add(a%kfl_massMatrix)
   call BBuffer%add(a%kfl_outBasis)

   call BBuffer%add(a%nsnap)
   call BBuffer%add(a%SnapshotsInterval)

   call BBuffer%add(a%kfl_Precondition)
   call BBuffer%add(a%nNonLinearIterations)
   call BBuffer%add(a%cotol)
   call BBuffer%add(a%basis_energy)
   call BBuffer%add(a%basis_number)
   call BBuffer%add(a%kfl_basis_number)
   call BBuffer%add(a%numBasisFiles)

   call BBuffer%add(a%kfl_SubstractMean)
   call BBuffer%add(a%param)
   
   !Adaptive
   call BBuffer%Add(len(a%RefinerErrorEstimator),a%RefinerErrorEstimator)
   call BBuffer%Add(len(a%RefinerErrorCriteria),a%RefinerErrorCriteria)
   call BBuffer%Add(a%RefinerErrorLimits)
   
   call BBuffer%BroadCast
   call BBuffer%Dealloc

end subroutine
