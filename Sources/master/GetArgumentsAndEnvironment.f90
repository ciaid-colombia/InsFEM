subroutine GetArgumentsAndEnvironment
   use Mod_Int2Str
   use typre
   use Mod_DCHashCharSize
   use Mod_DC_GeneralCase
   use def_master
   use Mod_DistributedContainer
   use def_parame
   use Mod_iofile
   use MPI
   use iso_c_binding 
   use Mod_Listen, only : upcase
   implicit none
   
   interface 
    function mkdir(path,mode) bind(c,name="mkdir") 
      use iso_c_binding 
      integer(c_int) :: mkdir 
      character(kind=c_char,len=1) :: path(*) 
      integer(c_int16_t), value :: mode 
    end function mkdir 
   end interface 
   
   
   integer(ip) :: icase,casnam,ind
   character(150) :: chararg1,chararg2,namda,aux,auxstr,auxstr2
   character(150) :: readType, DataFolder, ResultsFolder, PostprocessFolder, RestartFolder, OldDataFolder, OldRestartFolder
   character(150) :: auxbase,auxpos,auxres,auxrst,auxolddata,auxoldrestar
   integer(ip) :: narg,iarg,multicomm = 1_ip
   logical :: MultipleCases = .false.
   integer(ip) :: numCases
   character(DCHashCharSize) :: caseNames(5) = " "
   integer(ip) :: mymulticomm,multicommColor=0,ierr
   
   type(GeneralCase), pointer :: myCase => NULL()
   class(DistributedContainer), pointer :: myDC => NULL()
   character(150) :: fil_outpu
   
   integer(ip) :: kfl_multicase


   narg = command_argument_count()

   call GETARG(1,namda)

   !Default values for optional arguments
   readType = 'SERIAL'    !Default is root reads
   DataFolder = 'data'
   ResultsFolder   = 'results'
   PostprocessFolder    = 'post'
   RestartFolder   = 'rst'
   OldDataFolder = 'data'
   OldRestartFolder = 'rst'    
   
   !default number of cases is 1
   numCases = 1

   do iarg = 2,narg-1
       call GETARG(iarg,chararg1)
       call GETARG(iarg+1,chararg2)

       if (trim(chararg1) == '-cases') then
           !write(*,*)'chararg1 :',chararg1
           !write(*,*)'chararg2 :',chararg2
           numCases = str2int(chararg2)
           !allocate(character(len=150) :: caseNames(numCases))
           do casnam = 1,numCases
               call GETARG(iarg+1+casnam,aux)
               caseNames(casnam)= aux
           end do
           MultipleCases = .true.
       elseif (trim(chararg1) == '-datafolder') then
           DataFolder = chararg2
       elseif (trim(chararg1) == '-resultsfolder') then
           ResultsFolder   = chararg2
       elseif (trim(chararg1) == '-postprocessfolder') then
           PostprocessFolder    = chararg2
       elseif (trim(chararg1) == '-restartfolder') then
           RestartFolder   = chararg2      
       elseif (trim(chararg1) == '-olddatafolder') then
           OldDataFolder  = chararg2
       elseif (trim(chararg1) == '-oldrestarfolder') then
           OldRestartFolder  = chararg2    
       elseif (trim(chararg1) == '-readtype') then
           call To_upper(66,chararg2)
           readType = trim(chararg2)
           !SERIAL
           !PARTITIONED
        elseif (trim(chararg1) == '-multicomm') then
           multicomm = str2int(chararg2)           
       endif
   enddo

   !Running several cases at the same time
   !only used in PLCD topology optimization, stochastic analysis
   if (multicomm > 1) then
     !multicommColor = mod(MPIrank,multicomm)  !old, worse performance
     multicommColor = MPIrank/(MPIsize/multicomm)
     call MPI_Comm_split(MPIcomm,multicommColor,MPIrank,mymulticomm,ierr)
     
     !Replacing original MPI_COMM_WORLD by the new one!
     MPIcomm = mymulticomm
     
     call MPI_COMM_RANK( MPIcomm, MPIrank, ierr )
     call MPI_COMM_SIZE( MPIcomm, MPIsize, ierr )
     
     ResultsFolder = adjustl(trim(ResultsFolder))//'/results'//adjustl(trim(int2str(multicommColor)))
     !Build the new results folder if not existent
     ierr = mkdir(adjustl(trim(ResultsFolder))//char(0), int(o'772',c_int16_t))
     
     PostprocessFolder = adjustl(trim(PostprocessFolder))//'/post'//adjustl(trim(int2str(multicommColor)))
     !Build the new postprocess folder if not existent
     ierr = mkdir(adjustl(trim(PostprocessFolder))//char(0), int(o'772',c_int16_t)) 
     
     RestartFolder = adjustl(trim(RestartFolder))//'/rst'//adjustl(trim(int2str(multicommColor)))
     !Build the new restart folder if not existent
     ierr = mkdir(adjustl(trim(RestartFolder))//char(0), int(o'772',c_int16_t))  
   endif
    
   !Open the Global Output file
   !Compose filenames
   fil_outpu = trim(ResultsFolder)//'/'//adjustl(trim(namda))//adjustl(trim(int2str(MPIrank)))//'.glog'
   !Open         
   if (MPIrank == MPIroot) call iofile(zero,lun_outpu,fil_outpu,'RUN EVOLUTION')

   !Initialize Case List
   call CaseList%Initialize
   if (numCases > 1) then
     kfl_multicase = 1
   else
     kfl_multicase = 0
   endif
   !Set up icase
   do icase = 1,numCases 
     ind    = index(caseNames(icase), '/')
     auxstr = caseNames(icase)(ind+1:)
     !Only put the bar if there is a case name
     if (auxstr /= "") then
        !auxstr = '/'//auxstr
     endif

     auxbase = adjustl(trim(DataFolder))//'/'//trim(auxstr)
     auxres =  adjustl(trim(ResultsFolder))//'/'//trim(auxstr)
     !mkdir if non-existent
     ierr = mkdir(adjustl(trim(auxres))//char(0), int(o'772',c_int16_t))
     auxpos  = adjustl(trim(PostprocessFolder))//'/'//trim(auxstr)
     !mkdir if non-existent
     ierr = mkdir(adjustl(trim(auxpos))//char(0), int(o'772',c_int16_t))
     auxrst = adjustl(trim(RestartFolder))//'/'//trim(auxstr)
     !mkdir if non-existent
     ierr = mkdir(adjustl(trim(auxrst))//char(0), int(o'772',c_int16_t))
     auxolddata = adjustl(trim(OldDataFolder))//'/'//trim(auxstr)
     auxoldrestar = adjustl(trim(OldRestartFolder))//'/'//trim(auxstr)

     myCase => GeneralCase_Const() 
     call myCase%SetMPI(MPIcomm,MPIsize,MPIroot,MPIrank)
     call myCase%SetArgumentsAndEnvironment(namda,readType, auxbase, auxres, auxpos, auxrst, auxolddata, auxoldrestar, kfl_multicase)
     
     call myCase%SetMulticommData(multicomm,multicommColor)
     
     myDC => DC_GeneralCase_Const(myCase)
     auxstr2 = auxstr
     call upcase(auxstr2)
     call myDC%SetKey(auxstr2)   
     call CaseList%Add(myDC)  
   enddo 

   CasesTotal = numCases

end subroutine

subroutine To_upper(n,str)
   use typre
   character(*), intent(in out) :: str
   integer(ip) :: n,i

   do i = 1, n
      select case(str(i:i))
      case("a":"z")
         str(i:i) = achar(iachar(str(i:i))-32)
      end select
   end do 
end subroutine To_upper
