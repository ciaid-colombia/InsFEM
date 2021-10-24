subroutine nsc_restar(a,itask)
   ! This subroutine could be improved for implicit time integration high order schemes!
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file (itask=1)
   !    and writes the restart file when required.
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_NSCompressible
   use Mod_postpr
   implicit none
   class(NSCompressibleProblem) :: a
   integer(ip), intent(in)    :: itask
   integer(ip) :: ndime,npoin,npoinLocal,ierr,icomp,oldnpoin,oldnpoinLocal,ncomp
   real(rp), allocatable :: press(:,:),auxpress(:,:),auxveloc(:,:),tempe(:,:),auxtempe(:,:)
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

      call a%Memor%alloc(ndime,oldnpoin,auxveloc,'auxveloc','nsc_InterpolatedRestart')
      call a%Memor%alloc(1_ip,oldnpoin,auxpress,'auxpress','nsc_InterpolatedRestart')
      call a%Memor%alloc(1_ip,npoin,press,'press','nsc_InterpolatedRestart')
      call a%Memor%alloc(1_ip,oldnpoin,auxtempe,'auxtempe','nsc_InterpolatedRestart')
      call a%Memor%alloc(1_ip,npoin,tempe,'tempe','nsc_InterpolatedRestart')

      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif
      do icomp = 3,ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2

         nameV='Pressure'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxpress(1,1:oldnpoinLocal),nameV,nameG)

         nameV='Velocity'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxveloc(:,1:oldnpoinLocal),nameV,nameG)

         nameV='Temperature'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxtempe(1,1:oldnpoinLocal),nameV,nameG)

         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxpress)
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxtempe)
            call a%Int_Restart%Interpolate(1,auxpress,press)
            call a%Int_Restart%Interpolate(ndime,auxveloc,a%veloc(:,:,icomp))
            call a%Int_Restart%Interpolate(1,auxtempe,tempe) 
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxpress)
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxtempe)
            press = auxpress
            a%veloc(:,:,icomp) = auxveloc
            tempe = auxtempe
         endif
         a%press(:,icomp) = press(1,:)
         a%tempe(:,icomp) = tempe(1,:)
      end do
      do icomp = a%oldncomp+1,a%ncomp
         a%press(:,icomp) = a%press(:,a%oldncomp)
         a%veloc(:,:,icomp) = a%veloc(:,:,a%oldncomp)
         a%tempe(:,icomp) = a%tempe(:,a%oldncomp)
      end do

      if(a%kfl_repro>=1) then
         nameG= 'Repro'
         nameV= 'Repro'
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,a%repro,nameV,nameG)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime+2_ip,a%repro(:,:))
      endif

      a%press(:,1) = a%press(:,3)
      a%press(:,2) = a%press(:,3)
      a%veloc(:,:,1) = a%veloc(:,:,3)
      a%veloc(:,:,2) = a%veloc(:,:,3)
      a%tempe(:,1) = a%tempe(:,3)
      a%tempe(:,2) = a%tempe(:,3)

      call a%Memor%dealloc(ndime,oldnpoin,auxveloc,'auxveloc','nsc_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,oldnpoin,auxpress,'auxpress','nsc_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,npoin,press,'press','nsc_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,oldnpoin,auxtempe,'auxtempe','nsc_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,npoin,tempe,'tempe','nsc_InterpolatedRestart')

      call a%SpecificNSCompRestar(1)

   case (2)

      call a%SpecificNSCompRestar(2)

      do icomp = 3,a%ncomp

            write(nameG,"(A10,I1)") 'Time Step ',icomp-2

            nameV='Pressure'
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%press(1:npoinLocal,icomp),nameV,nameG)
   
            nameV='Velocity'
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%veloc(:,1:npoinLocal,icomp),nameV,nameG)
   
            nameV='Temperature'
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%tempe(1:npoinLocal,icomp),nameV,nameG)
      enddo


      if(a%kfl_repro>=1) then
         nameG= 'Repro'
         nameV= 'Repro'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%repro,nameV,nameG)
      endif

   end select
   
end subroutine nsc_restar
