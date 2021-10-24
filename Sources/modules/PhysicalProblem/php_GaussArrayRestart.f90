subroutine GaussArrayRestart(a,array,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file 
   !    corresponding to gauss point arrays and writes them to a new restart
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_iofile
   use Mod_int2str
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip) :: ndime,npoin,oldNpoin,ipoin,idofn,ndofn,ielem,nelem
   type(r3p)    :: array(:)
   real(rp)    :: rdummy
   integer(ip) :: iaux,iaux2,iaux3
   real(rp), allocatable :: auxvesgs(:)
   
   call a%Timer%Total%Tic
   call a%Timer%Restar%Tic
  

   select case (itask)

   case (1)
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)

      iaux = 0
      do ielem = 1,nelem
         iaux = iaux+size(array(ielem)%a,3)*ndime
      enddo
      call a%Memor%alloc(iaux,auxvesgs,'auxvesgs','nsi_restar')
            
      read (a%lun_rstar) (auxvesgs(ipoin),ipoin=1,OldNpoin)
                     
      iaux = 0
      do ielem = 1,nelem
         iaux2 = size(array(ielem)%a,3)
         do iaux3 = 1,iaux2
            array(ielem)%a(1:ndime,1,iaux3)= auxvesgs(iaux+1:iaux+ndime)
            iaux = iaux+ndime
         enddo
      enddo
      call a%Memor%dealloc(iaux,auxvesgs,'auxvesgs','nsi_restar')
      
   case (2)

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNelem(nelem)

      iaux = 0
      do ielem = 1,nelem
         iaux = iaux+size(array(ielem)%a,3)*ndime
      enddo
      call a%Memor%alloc(iaux,auxvesgs,'auxvesgs','nsi_restar')
         
      iaux = 0
      do ielem = 1,nelem
         iaux2 = size(array(ielem)%a,3)
         do iaux3 = 1,iaux2
           auxvesgs(iaux+1:iaux+ndime) = array(ielem)%a(1:ndime,1,iaux3)
           iaux = iaux+ndime
         enddo
      enddo
         
      write (a%lun_rstar) (auxvesgs(ipoin),ipoin=1,npoin)
         
      call a%Memor%dealloc(iaux,auxvesgs,'auxvesgs','nsi_restar')         

   end select
   
   call a%Timer%Total%Toc
   call a%Timer%Restar%Toc

100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)

end subroutine GaussArrayRestart
