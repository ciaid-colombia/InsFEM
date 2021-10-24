subroutine rom_SVD(a)
   use typre
   use Mod_Iofile
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   character(150) :: fil_outpu
   character(3)   :: exmod
   integer(ip)    :: lun_outpu

   exmod = adjustl(trim(a%exmod))
   select case(a%kfl_eigentype)
   case('EPS')
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Solve EPS system, Begin'
      call a%EigenSystem%setMatrixPointersBuild
      call a%EigenSystem%SolveSystem
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Solve EPS system, End'
   case('SVD')
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Solve SVD system, Begin'
      call a%EigenSystem%SolveSystem
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Solve SVD system, End'
   end select
      
   fil_outpu = trim(a%OutputFolder)//'/'//adjustl(trim(a%namda))//'.'//adjustl(trim(exmod))//'.pod.sol'
   if (a%MPIrank == a%MPIroot) call iofile(0,lun_outpu,fil_outpu,'podsolve','replace') ! open file
   call a%EigenSystem%WriteInfo(fil_outpu)
      
   end subroutine
