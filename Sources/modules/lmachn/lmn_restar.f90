subroutine lmn_restar(a,itask)
   use MPI
   use typre
   use def_parame
   use Mod_LowMach
   implicit none
   class(LowMachProblem)   :: a
   integer(ip), intent(in) :: itask
   integer(ip) :: ndime,npoin,npoinLocal,ierr,icomp,oldnpoin,oldnpoinLocal,ncomp
   real(rp), allocatable :: press(:,:),auxpress(:,:),auxveloc(:,:),tempe(:,:),auxtempe(:,:)
   real(rp) :: array(1)
   character(150) :: nameV,nameG

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoinLocal(npoinLocal)

   select case (itask)

   case (1)
      
      if (a%kfl_inter == 1) then
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(ndime,oldnpoin,auxveloc,'auxveloc','lmn_InterpolatedRestart')
      call a%Memor%alloc(1_ip,oldnpoin,auxpress,'auxpress','lmn_InterpolatedRestart')
      call a%Memor%alloc(1_ip,npoin,press,'press','lmn_InterpolatedRestart')
      call a%Memor%alloc(1_ip,oldnpoin,auxtempe,'auxtempe','lmn_InterpolatedRestart')
      call a%Memor%alloc(1_ip,npoin,tempe,'tempe','lmn_InterpolatedRestart')
     
      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif
      do icomp = 3,ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         nameV = 'ThermPress'
         call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,array,nameV,nameG)
         a%pther(icomp) = array(1)
         call MPI_BCAST(a%pther(icomp),1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)

         nameV = 'Velocity'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxveloc(:,1:oldnpoinLocal),nameV,nameG)

         nameV = 'Temperature'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxtempe(1,1:oldnpoinLocal),nameV,nameG)

         nameV = 'Pressure'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxpress(1,1:oldnpoinLocal),nameV,nameG)

         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxtempe)
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxpress)
            call a%Int_Restart%Interpolate(ndime,auxveloc,a%veloc(:,:,icomp))
            call a%Int_Restart%Interpolate(1,auxtempe,tempe) 
            call a%Int_Restart%Interpolate(1,auxveloc,press)
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxtempe)
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxpress)
            a%veloc(:,:,icomp) = auxveloc
            press = auxpress
            tempe = auxtempe
         endif
         a%press(:,icomp) = press(1,:)
         a%tempe(:,icomp) = tempe(1,:)
      enddo
      do icomp = a%oldncomp+1,a%ncomp
         a%press(:,icomp) = a%press(:,a%oldncomp)
         a%veloc(:,:,icomp) = a%veloc(:,:,a%oldncomp)
         a%tempe(:,icomp) = a%tempe(:,a%oldncomp)
      end do
      
      a%press(:,1) = a%press(:,3)
      a%press(:,2) = a%press(:,3)
      a%veloc(:,:,1) = a%veloc(:,:,3)
      a%veloc(:,:,2) = a%veloc(:,:,3)
      a%tempe(:,1) = a%tempe(:,3)
      a%tempe(:,2) = a%tempe(:,3)
      a%pther(1) = a%pther(3)
      a%pther(2) = a%pther(3)

      call a%Memor%dealloc(ndime,oldnpoin,auxveloc,'auxveloc','lmn_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,oldnpoin,auxpress,'auxpress','lmn_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,npoin,press,'press','lmn_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,oldnpoin,auxtempe,'auxtempe','lmn_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,npoin,tempe,'tempe','lmn_InterpolatedRestart')
     
   case (2)
      
      do icomp = 3,a%ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         nameV = 'ThermPress'
         array = a%pther(icomp)
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,array,nameV,nameG)

         nameV = 'Velocity'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%veloc(:,1:npoinLocal,icomp),nameV,nameG)

         nameV = 'Temperature'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%tempe(1:npoinLocal,icomp),nameV,nameG)

         nameV = 'Pressure'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%press(1:npoinLocal,icomp),nameV,nameG)
      enddo

   end select

end subroutine lmn_restar
