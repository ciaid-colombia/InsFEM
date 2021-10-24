subroutine php_ReadOnFunctions(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   character(150) :: outstr
   integer(ip) :: ifunc,ifunp,nfunp

   !Output
   outstr = adjustl(trim(a%exmod))//'_RootReadOnFunctions'
   
   call a%Listener%listen(outstr)
   do while(a%Listener%words(1)/='ENDFU')
      if(a%kfl_conbc==0) then
         ifunc=a%Listener%getint('FUNCT',1,'#FUNCTION NUMBER')
         if(ifunc<0.or.ifunc>10) then
            call runend(outstr//': WRONG FUNCTION NUMBER')
         else
            if(a%Listener%words(2)=='PARAB') then
               a%kfl_funty(ifunc,1)=1
               a%kfl_funty(ifunc,2)=6
            else if(a%Listener%words(2)=='PERIO') then
               a%kfl_funty(ifunc,1)=2
               a%kfl_funty(ifunc,2)=6
            else if(a%Listener%words(2)=='DISCR') then
               a%kfl_funty(ifunc,1)=3
               a%kfl_funty(ifunc,2)=2*a%Listener%getint('NUMBE',1,'#Number of data')
            else if(a%Listener%words(2)=='TSTRE') then
               a%kfl_funty(ifunc,1)=6 ! function identifier
               a%kfl_funty(ifunc,2)=6 ! number of parameters, they put six to all
            else if(a%Listener%words(2)=='ROSEN') then
               a%kfl_funty(ifunc,1)=7 ! function identifier
               a%kfl_funty(ifunc,2)=6 ! number of parameters, they put six to all
            else if(a%Listener%words(2)=='CONST') then
               a%kfl_funty(ifunc,1)=8 ! function identifier
               a%kfl_funty(ifunc,2)=6 ! number of parameters, they put six to all
            else if(a%Listener%words(2)=='VELOC') then
               a%kfl_funty(ifunc,1)=9 ! function identifier
               a%kfl_funty(ifunc,2)=6 ! number of parameters, they put six to all     
           elseif (a%Listener%words(2) == 'HYPER') then !hyperbolic tangent
               a%kfl_funty(ifunc,1) = 10
               a%kfl_funty(ifunc,2) = 1
           elseif (a%Listener%words(2) == 'SINUS') then !Sinusoidal function
               a%kfl_funty(ifunc,1) = 11
               a%kfl_funty(ifunc,2) = 1
           elseif (a%Listener%words(2) == 'FSIWA') then !Sinusoidal function
               a%kfl_funty(ifunc,1) = 12
               a%kfl_funty(ifunc,2) = 1
            else
               a%kfl_funty(ifunc,1)=0
               a%kfl_funty(ifunc,2)=6
            end if 
            
            call a%SpecificReadOnFunctions(ifunc)
            
            if(a%kfl_funty(ifunc,1)/=0) then
               call a%Memor%palloc(a%kfl_funty(ifunc,2),a%funpa(ifunc)%a,'funpa%a',outstr)
               if(a%kfl_funty(ifunc,1)==3) then
                  ifunp=0
                  nfunp=a%kfl_funty(ifunc,2)/2
                  if(nfunp<1) call runend(outstr//': WRONG DISCRETE FUNCTION PARAMETER')
                  call a%Listener%listen(outstr)
                  do while(a%Listener%words(1)/='ENDFU')
                     ifunp=ifunp+1
                     if(ifunp>nfunp) call runend(outstr//': WRONG DISCRETE FUNCTION DATA')
                     a%funpa(ifunc)%a((ifunp-1)*2+1)=a%Listener%param(1)
                     a%funpa(ifunc)%a((ifunp-1)*2+2)=a%Listener%param(2)
                     call a%Listener%listen(outstr)
                  end do
                  ! Order the function field
                  call ordena(nfunp,a%funpa(ifunc)%a)
               else
                  a%funpa(ifunc)%a(1:a%kfl_funty(ifunc,2))=a%Listener%param(3:2+a%kfl_funty(ifunc,2))
               end if
            end if
         end if
      end if
      call a%Listener%listen(outstr)
   end do


end subroutine



subroutine php_ScatterOnFunctions(a)
   use MPI
   use typre
   use Mod_PhysicalProblem
   use Mod_ToCSR
   implicit none
   class(PhysicalProblem) :: a
   
   !MPI
   integer, parameter :: mtag0 = 0, mtag1 = 1, mtag2 = 2, mtag3 = 3,mtag4 = 4, mtag5 = 5, mtag6 = 6, mtag7 = 7
   integer status(MPI_STATUS_SIZE)
   integer :: ierr,irequest1(a%MPIsize),irequest2(a%MPIsize),irequest3(a%MPIsize),irequest4(a%MPIsize),irequest5(a%MPIsize),irequest6(a%MPIsize),irequest7(a%MPIsize),irequest0(a%MPIsize)
   
   integer(ip) :: funpa_coun,ifunc,isize
   real(rp),allocatable :: aux_glfunpa(:)
   integer(ip), allocatable :: aux_gpfunpa(:)
   
   character(150) :: outstr
   
   !Output
   outstr = adjustl(trim(a%exmod))//'_ScatterOnFunctionsData'
   
   if (a%kfl_conbc/=1) then
      if (a%MPIrank == a%MPIroot) then
         funpa_coun = 0
         do ifunc = 1,10
            if (associated(a%funpa(ifunc)%a)) then
               funpa_coun = funpa_coun + size(a%funpa(ifunc)%a)
            endif   
         enddo
      endif
      CALL MPI_BCAST(a%kfl_funty, size(a%kfl_funty), MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      !Allocate auxiliar arrays
      CALL MPI_BCAST(funpa_coun, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      call a%Memor%alloc(11,aux_gpfunpa,'aux_gpfunpa',outstr)
      call a%Memor%alloc(funpa_coun,aux_glfunpa,'aux_glfunpa',outstr)
      
      if (a%MPIrank == a%MPIroot) then
         call r1p2CSR(10,a%funpa,aux_gpfunpa,aux_glfunpa,a%Memor,'funpa%a')
      endif
      CALL MPI_BCAST(aux_gpfunpa,11, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      CALL MPI_BCAST(aux_glfunpa, funpa_coun, MPI_REAL8, a%MPIroot, a%MPIcomm, ierr)
      
      !if (a%MPIrank /= a%MPIroot) then
         do ifunc = 1,10
            isize = aux_gpfunpa(ifunc+1)-aux_gpfunpa(ifunc)
            if (isize /= 0) then
               call a%Memor%palloc(isize,a%funpa(ifunc)%a,'funpa%a',outstr)
               a%funpa(ifunc)%a = aux_glfunpa(aux_gpfunpa(ifunc):aux_gpfunpa(ifunc+1)-1)
            endif
         enddo
      !endif
      call a%Memor%dealloc(11,aux_gpfunpa,'aux_gpfunpa',outstr)
      call a%Memor%dealloc(funpa_coun,aux_glfunpa,'aux_glfunpa',outstr)
      
   endif
   
   call a%SpecificScatterOnFunctions

end subroutine

