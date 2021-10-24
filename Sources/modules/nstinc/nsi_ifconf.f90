subroutine nsi_ifconf(a)
!-----------------------------------------------------------------------
! NAME 
!    nsi_ifconf
! DESCRIPTION
!    This routine checks if the flow is confined or not and looks for a node
!    to impose pressure 
!-----------------------------------------------------------------------
   use MPI
   use typre
   use Mod_NavierStokes
   use Mod_NsiExacso   
   implicit none
   class(NavierStokesProblem) :: a
   type(NsiExacso) :: exacso
   
   integer(ip) :: kbopo,ndime,npoin,idime,ipoin,nbopo,irank
   integer(ip), allocatable :: aux_kfl_confi(:)
   real(rp)     :: dummr
   integer(ip)  :: dummi
   real(rp), pointer       :: coord(:)  => NULL() 
   !Exact Values
   real(rp)                :: expre
   real(rp), allocatable   :: exprg(:)

   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2
   integer status(MPI_STATUS_SIZE)
   integer :: ierr,irequest
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNbopo(nbopo)
   
   call a%Memor%alloc(ndime,exprg,'exprg','nsi_ifconf')   
   
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
      call a%Memor%alloc(a%MPIsize,aux_kfl_confi,'aux_kfl_confi','nsi_ifconf')
      do irank = 0,a%MPIsize-1
         call MPI_RECV(aux_kfl_confi(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,status,ierr)
      enddo
      a%kfl_confi = minval(aux_kfl_confi)
      call a%Memor%dealloc(a%MPIsize,aux_kfl_confi,'aux_kfl_confi','nsi_ifconf')
   endif
   call MPI_WAIT(irequest, status, ierr)
   call MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   !If pressure node is not prescribed, find one and impose p=0 or 
   !the value given by the exact solution case.
   if(a%nodpr==-1.and.a%kfl_confi==1) then
      a%nodpr=1
      call a%Mesh%Global2Local(a%nodpr,a%nodpr)
      
      if(a%kfl_exacs/=0) then
	 
	 call a%Mesh%GetPointCoord(a%nodpr,coord)      
      
         call exacso%nsi_ComputeSolution(ndime,coord,a%ctime,a)
         call exacso%nsi_GetPressure(ndime,expre,exprg)              
               
         a%prepr=expre
          
      else
         a%prepr=0.0_rp
      end if
   end if
   
   call a%Memor%dealloc(ndime,exprg,'exprg','nsi_ifconf')     

end subroutine nsi_ifconf
