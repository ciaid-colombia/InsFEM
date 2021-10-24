subroutine php_restar(a,itask)
   use MPI
   use typre
   use def_parame
   use Mod_int2str
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip), intent(in)    :: itask
   integer(ip)             :: icomp,ierr
   integer(ip)             :: arrayaux(1)
   character(150)          :: nameV,nameG
   character(3)            :: exmod

   call a%Timer%Total%Tic
   call a%Timer%Restar%Tic

   call a%OpenFilesRestart(itask)
   exmod =adjustl(trim(a%exmod))

   select case (itask)

   case (1)

      if(a%MPIrank==a%MPIroot) write(*,*) adjustl(trim(exmod))//': Reading restart, Begin'
      write(nameG,"(A10,I1)") 'Time Step ',1
      nameV = 'TimeIntegrator'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,arrayaux,nameV,nameG)
      a%oldncomp = arrayaux(1)

      CALL MPI_BCAST(a%oldncomp,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)

      if (a%ncomp <= a%oldncomp) then
         a%neule = 0_ip
      else
         a%neule = a%ncomp - a%oldncomp
      endif
      
      call a%SpecificRestart(one)
      if(a%MPIrank==a%MPIroot) write(*,*) adjustl(trim(exmod))//': Reading restart, End'

   case (2)

      if(a%MPIrank==a%MPIroot) write(*,*) adjustl(trim(exmod))//': Writing restart, Begin'
      write(nameG,"(A10,I1)") 'Time Step ',1
      nameV = 'TimeIntegrator'
      arrayaux(1) = a%ncomp
      call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,arrayaux,nameV,nameG)
     
      call a%SpecificRestart(two)
      if(a%MPIrank==a%MPIroot) write(*,*) adjustl(trim(exmod))//': Writing restart, End'

   end select
   
   call a%CloseFilesRestart(itask)     

   call a%Timer%Total%Toc
   call a%Timer%Restar%Toc

end subroutine php_restar
