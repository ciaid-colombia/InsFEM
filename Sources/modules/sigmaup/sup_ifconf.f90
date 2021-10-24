subroutine sup_ifconf(a)
!-----------------------------------------------------------------------
! NAME 
!    nsi_ifconf
! DESCRIPTION
!    This routine checks if the flow is confined or not and look for a node
!    to impose pressure 
!-----------------------------------------------------------------------
   use typre
   use Mod_Memor
   use Mod_Listen
   use Mod_MPIObject
   use Mod_Mesh
   use Mod_ThreeField
   use MPI
   use Mod_SupExacso     
   implicit none
   class(ThreeFieldNSProblem), target :: a
   type(SupExacso) :: exacso   
   
   integer(ip) :: kbopo,ndime,npoin,idime,ipoin,nbopo,irank,nbcu
   integer(ip), allocatable :: aux_kfl_confi(:)
   real(rp)     :: dummr
   integer(ip)  :: dummi
   real(rp), pointer       :: coord(:)   
   !Exact Values
   real(rp)                :: expre
   real(rp), allocatable   :: exprg(:)   

   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2
   integer status(MPI_STATUS_SIZE)
   integer :: ierr,irequest,aux_confi
   
   type(MemoryMan), pointer :: Memor
   
   !to do multy materials
   integer(ip) :: imat=1
   
   Memor => a%Memor
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNbopo(nbopo)
   
   call Memor%alloc(ndime,exprg,'exprg','sup_ifconf')     
   
   if(a%MatProp(imat)%lawvi>=0)then
      nbcu=0
   elseif(a%MatProp(imat)%lawvi<0)then
      nbcu=(ndime-1)*(ndime-1)+2
   end if

   !Check if flow is confined
   if(a%kfl_confi/=1.and.a%kfl_confi/=-1) then
      a%kfl_confi=0
      do ipoin=1,npoin
         dimensions1: do idime=1,ndime
            if(a%kfl_fixno(nbcu + idime,ipoin)==1) then
               a%kfl_confi=a%kfl_confi+1
            end if
         end do dimensions1
      end do
      
      kbopo=nbopo*ndime
      if(a%kfl_confi==kbopo) then
         a%kfl_confi=1
      else
         a%kfl_confi=0
      end if
   end if  
   
   call MPI_ISEND(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, mtag1, a%MPIcomm,irequest, ierr)
   
   if (a%MPIrank == a%MPIroot) then
      call Memor%alloc(a%MPIsize,aux_kfl_confi,'aux_kfl_confi','nsi_ifconf')
      do irank = 0,a%MPIsize-1
         call MPI_RECV(aux_kfl_confi(irank+1),1,MPI_INTEGER4,irank,mtag1,a%MPIcomm,status,ierr)
      enddo
      a%kfl_confi = minval(aux_kfl_confi)
      call Memor%dealloc(a%MPIsize,aux_kfl_confi,'aux_kfl_confi','nsi_ifconf')
   endif
   call MPI_WAIT(irequest, status, ierr)
   call MPI_BCAST(a%kfl_confi, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   !If pressure node is not prescribed, find one and impose p=0 or 
   !the value given by the exact solution case.
   if(a%nodpr==-1.and.a%kfl_confi==1) then
      a%nodpr=1
      call a%Mesh%Global2Local(a%nodpr,a%nodpr)
      
      call a%Mesh%GetPointCoord(a%nodpr,coord)      

      if(a%kfl_exacs/=0) then
      
         call exacso%sup_ComputeSolution(ndime,coord,a%ctime,a%LogFormulation,a)
         call exacso%sup_GetPressure(ndime,expre,exprg)   
               
          a%prepr=expre                   
          
       else
         a%prepr=0.0_rp
       end if
   end if
   
   call Memor%dealloc(ndime,exprg,'exprg','sup_ifconf')       

end subroutine sup_ifconf
