subroutine nsc_ifconf(a)
!-----------------------------------------------------------------------
! NAME 
!    nsc_ifconf
! DESCRIPTION
!    This routine checks if the flow is confined or not
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_Mesh
   use MPI
   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip), allocatable :: aux_kfl_confi(:)
   integer, parameter       :: mtag1 = 1, mtag2 = 2
   integer(ip) :: kbopo,ndime,npoin,idime,ipoin,nbopo,irank
   integer     :: ierr,irequest
   integer status(MPI_STATUS_SIZE)
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNbopo(nbopo)
   
   !Check if flow is confined
   if(a%kfl_confi/=1.and.a%kfl_confi/=-1) then
      a%kfl_confi=0
      do ipoin=1,npoin
         dimensions1: do idime=1,ndime
            if(a%kfl_fixno(idime,ipoin)==1) then
               a%kfl_confi=a%kfl_confi+1
               exit dimensions1
            end if
         end do dimensions1
      end do
      kbopo=nbopo
      if(a%kfl_confi==kbopo) then
         a%kfl_confi=1
      else
         a%kfl_confi=0
      end if
   end if
   
   call MPI_ISEND(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,irequest, ierr)
   
   if (a%MPIrank == a%MPIroot) then
      call a%Memor%alloc(a%MPIsize,aux_kfl_confi,'aux_kfl_confi','nsc_ifconf')
      do irank = 0,a%MPIsize-1
         call MPI_RECV(aux_kfl_confi(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,status,ierr)
      enddo
      a%kfl_confi = minval(aux_kfl_confi)
      call a%Memor%dealloc(a%MPIsize,aux_kfl_confi,'aux_kfl_confi','nsc_ifconf')
   endif
!   call MPI_WAIT(irequest, status, ierr)
   call MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   !If pressure node is not prescribed, find one and impose p=0 or 
   !the value given by the exact solution case.
   if(a%nodpr==-1.and.a%kfl_confi==1) then
      a%nodpr=1
      call a%Mesh%Global2Local(a%nodpr,a%nodpr)
      
!      if(a%kfl_exacs/=0) then
        
!         call a%Mesh%GetPointCoord(a%nodpr,coord)      
      
!         call exacso%nsi_ComputeSolution(ndime,coord,a)
!         call exacso%nsi_GetPressure(ndime,expre,exprg)              
               
!          a%prepr=0.0_rp!expre                   
          
!      else
         a%prepr=0.0_rp
!      end if
   end if

end subroutine nsc_ifconf
