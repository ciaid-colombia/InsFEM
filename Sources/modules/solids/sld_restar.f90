subroutine sld_restar(a,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file (itask=1)
   !    and writes the restart file when required.
   !------------------------------------------------------------------------
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem)     :: a
   integer(ip), intent(in)  :: itask
   real(rp),    allocatable :: auxdisp(:,:),auxveloc(:,:),auxaccel(:,:)
	real(rp), allocatable    :: restar_aux(:,:)
   integer(ip)              :: ndime,idime,npoin,npoinLocal,ierr,icomp,oldnpoin,oldnpoinLocal,ncomp
   real(rp)                 :: beta,omega
   integer(ip)              :: auxcoupling(1)=0,auxdynamic(1)
   character(150)           :: nameV,nameG

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoinLocal(npoinLocal)

   select case (itask)

   case (1)
      
      if (a%kfl_inter == 1) then
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoin(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(ndime,oldnpoin,restar_aux,'restar_aux','sld_restar')
      call a%Memor%alloc(ndime,oldnpoin,auxdisp,'auxdisp','sld_InterpolatedRestart')
      if (a%kfl_timei==1 ) then
         call a%Memor%alloc(ndime,oldnpoin,auxveloc,'auxveloc','sld_InterpolatedRestart')
         call a%Memor%alloc(ndime,oldnpoin,auxaccel,'auxaccel','sld_InterpolatedRestart')
      endif
     
      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A7)") 'Dynamic'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxdynamic,nameV,nameG)

      do icomp = 3,ncomp
          write(nameG,"(A10,I1)") 'Time Step ',icomp-2
          do idime=1, ndime
              write(nameV,"(A13,I1,A1,I1)") 'Displacement_',idime,'_',icomp-2
              call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxdisp(idime,1:oldnpoinLocal),nameV,nameG)
          enddo

         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxdisp)
            call a%Int_Restart%Interpolate(ndime,auxdisp,a%disp(:,:,icomp))
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxdisp)
            a%disp(:,:,icomp) = auxdisp
         endif
      enddo
      do icomp = a%oldncomp+1,a%ncomp
         a%disp(:,:,icomp) = a%disp(:,:,a%oldncomp)
      end do

      if (auxdynamic(1)==1 ) then

         do icomp = 1,2
            write(nameG,"(A10,I1)") 'Time Step ',icomp

            do idime=1, ndime
               write(nameV,"(A9,I1,A1,I1)") 'Velocity_',idime,'_',icomp
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxveloc(idime,1:oldnpoinLocal),nameV,nameG)
               write(nameV,"(A13,I1,A1,I1)") 'Acceleration_',idime,'_',icomp
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxaccel(idime,1:oldnpoinLocal),nameV,nameG)
            enddo

            if (a%kfl_inter == 1) then
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxaccel)
               call a%Int_Restart%Interpolate(ndime,auxveloc,a%veloc(:,:,icomp))
               call a%Int_Restart%Interpolate(ndime,auxaccel,a%accel(:,:,icomp))
            else
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxaccel)
               a%veloc(:,:,icomp) = auxveloc
               a%accel(:,:,icomp) = auxaccel
            endif
         end do

      endif

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Coupling'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

      if (auxcoupling(1)==1 ) then

         do idime= 1,ndime
            write(nameV,"(A9,I1)") 'Displ_cp_',idime
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(idime,:),nameV,nameG)
         end do
         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
            call a%Int_Restart%Interpolate(ndime,restar_aux,a%disp_cp)         
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
            a%disp_cp = restar_aux
         endif

      endif

      call a%Memor%dealloc(ndime,oldnpoin,restar_aux,'restar_aux','sld_restar')
      call a%Memor%dealloc(ndime,oldnpoin,auxdisp,'auxdisp','sld_InterpolatedRestart')
      if (a%kfl_timei==1 ) then
         call a%Memor%dealloc(ndime,oldnpoin,auxveloc,'auxveloc','sld_InterpolatedRestart')
         call a%Memor%dealloc(ndime,oldnpoin,auxaccel,'auxaccel','sld_InterpolatedRestart')
      endif

      !FSI does not need index update as it comes from the coupling
      if (a%kfl_timei==1 ) then
         if (auxcoupling(1)==0 ) then
         else
             a%kfl_FSI_restart = .true.
         endif
      endif

      if (a%sld_type == 'NONLI') then
         !-------------------------------------------------------------------
         !MOVE mesh so we are at initial configuration after restart
         !Dealloc and recompute ExtnorLpoty and Vmass
         call a%Mesh%DeallocExnorLpoty
         call a%Mesh%DeallocVmass

         !Recompute
         call a%Mesh%ComputeVmass
         call a%Mesh%ExtnorLpoty
      endif

   case (2)
      
      do icomp = 3,a%ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         do idime= 1,ndime
            write(nameV,"(A13,I1,A1,I1)") 'Displacement_',idime,'_',icomp-2
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%disp(idime,1:npoinLocal,icomp),nameV,nameG)
         end do
      enddo

      write(nameG,"(A10,I1)") 'Time Step ',1
      if (a%kfl_timei==1 ) then
         auxdynamic=1
         write(nameV,"(A7)") 'Dynamic'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxdynamic,nameV,nameG)

         do icomp = 1,2
            write(nameG,"(A10,I1)") 'Time Step ',icomp
            do idime= 1,ndime
               write(nameV,"(A9,I1,A1,I1)") 'Velocity_',idime,'_',icomp
               call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%veloc(idime,1:npoinLocal,icomp),nameV,nameG)

               write(nameV,"(A13,I1,A1,I1)") 'Acceleration_',idime,'_',icomp
               call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%accel(idime,1:npoinLocal,icomp),nameV,nameG)
            end do
         end do

      else
         auxdynamic=0
         write(nameV,"(A7)") 'Dynamic'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxdynamic,nameV,nameG)
      endif

      write(nameG,"(A10,I1)") 'Time Step ',1
      if (a%kfl_docoupconv) then
         auxcoupling=1
         write(nameV,"(A8)") 'Coupling'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)
         do idime= 1,ndime
            write(nameV,"(A9,I1)") 'Displ_cp_',idime
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%disp_cp(idime,:),nameV,nameG)
         end do
      else
         auxcoupling=0
         write(nameV,"(A8)") 'Coupling'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

      endif

   end select

   call a%SolidSpecificRestart(itask)

end subroutine sld_restar
