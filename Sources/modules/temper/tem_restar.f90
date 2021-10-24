subroutine tem_restar(a,itask)
   use typre
   use def_parame
   use MPI
   use Mod_Temperature
   implicit none
   class(TemperatureProblem) :: a
   integer(ip), intent(in) :: itask
   integer(ip) :: ndime,npoin,npoinLocal,ierr,icomp,oldnpoin,oldnpoinLocal,ncomp
   real(rp), allocatable :: tempe(:,:),auxtempe(:,:),auxrepro(:,:)
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

      call a%Memor%alloc(1_ip,oldnpoin,auxtempe,'auxtempe','nsc_InterpolatedRestart')
      call a%Memor%alloc(1_ip,npoin,tempe,'tempe','nsc_InterpolatedRestart')

      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif
      do icomp = 3,ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         nameV = 'Temperature'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxtempe(1,1:oldnpoinLocal),nameV,nameG)
         
         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxtempe)
            call a%Int_Restart%Interpolate(1,auxtempe,tempe) 
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxtempe)
            tempe = auxtempe
         endif
         a%tempe(:,icomp) = tempe(1,:)
      enddo
      do icomp = a%oldncomp+1,a%ncomp
         a%tempe(:,icomp) = a%tempe(:,a%oldncomp)
      end do

!      if(a%kfl_repro>=1) then
!         call a%Memor%alloc(1,oldnpoin,auxrepro,'auxrepro','tem_InterpolatedRestart')
!         write(nameG,"(A5)") 'Repro'
!         write(nameV,"(A5)") 'Repro'
!         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxrepro,nameV,nameG)
!
!         if (a%kfl_inter == 1) then
!            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxrepro)
!            call a%Int_Restart%Interpolate(1,auxrepro,a%repro)
!         else
!            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxrepro)
!            a%repro = auxrepro(1,:)
!         endif
!         call a%Memor%dealloc(1,oldnpoin,auxrepro,'auxrepro','tem_InterpolatedRestart')
!         write(nameV,"(A5)") 'Repro'
!         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,a%repro,nameV,nameG)
!         call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%repro)
!      endif

      a%tempe(:,1) = a%tempe(:,3)
      a%tempe(:,2) = a%tempe(:,3)
      
      call a%Memor%dealloc(1_ip,oldnpoin,auxtempe,'auxtempe','nsc_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,npoin,tempe,'tempe','nsc_InterpolatedRestart')

   case (2)
      
      do icomp = 3,a%ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         nameV = 'Temperature'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%tempe(1:npoinLocal,icomp),nameV,nameG)
      enddo
      
      if(a%kfl_repro>=1) then
         write(nameV,"(A5)") 'Repro'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%repro,nameV,nameG)
      endif

!      if (a%kfl_trasg >= 1 .and. a%kfl_iofor == 0) then            
!         call a%GaussArrayRestart(a%tesgs,2)
!      endif

   end select

end subroutine
