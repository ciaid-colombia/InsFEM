subroutine sldup_restar(a,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file (itask=1)
   !    and writes the restart file when required.
   !------------------------------------------------------------------------
   use typre
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem) :: a
   integer(ip), intent(in) :: itask
   real(rp), allocatable   :: auxpress(:)
   real(rp), allocatable   :: restar_aux(:)
   integer(ip)             :: ndime,idime,npoin,npoinLocal,icomp
   integer(ip)             :: oldnpoin,oldnpoinLocal,ncomp
   integer(ip)             :: auxcoupling(1)=0
   character(150)          :: nameV,nameG

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

      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif

      call a%Memor%alloc(oldnpoin,restar_aux,'restar_aux','sldsup_restar')
      call a%Memor%alloc(oldnpoin,auxpress,'auxpress','sldsup_InterpolatedRestart')

      do icomp = 3,ncomp

          write(nameG,"(A10,I1)") 'Time Step ',icomp-2
          write(nameV,"(A9,I1)") 'Pressure_',icomp-2
          call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxpress(1:oldnpoinLocal),nameV,nameG)

         if (a%kfl_inter == 1) then

            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxpress)
            call a%Int_Restart%Interpolate(1,auxpress,a%press(:,icomp))
         else

            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxpress)
            a%press(:,icomp) = auxpress
         endif
      enddo
      do icomp = a%oldncomp+1,a%ncomp
         a%press(:,icomp)   = a%press(:,a%oldncomp)
      end do
      
      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Coupling'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

      if (auxcoupling(1)==1 ) then
          if(a%kfl_docoupconv)  then

              restar_aux=0.0_rp
              write(nameV,"(A8)") 'Press_cp'
              call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(:),nameV,nameG)
              if (a%kfl_inter == 1) then
                  call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,restar_aux(:))
                  call a%Int_Restart%Interpolate(1,restar_aux(:),a%press_cp(:))         
              else
                  call a%Mesh%ArrayCommunicator%GhostCommunicate(1,restar_aux(:))
                  a%press_cp(:) = restar_aux(:)
              endif

          endif
      endif

      call a%Memor%dealloc(oldnpoin,restar_aux,'restar_aux','sldsup_restar')
      call a%Memor%dealloc(oldnpoin,auxpress,'auxpress','sldsup_InterpolatedRestart')

   case (2)

      do icomp = 3,a%ncomp
          write(nameG,"(A10,I1)") 'Time Step ',icomp-2
          write(nameV,"(A9,I1)") 'Pressure_',icomp-2
          call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%press(1:npoinLocal,icomp),nameV,nameG)
      enddo

      write(nameG,"(A10,I1)") 'Time Step ',1
      if(a%kfl_docoupconv) then
          write(nameV,"(A8)") 'Press_cp'
          call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%press_cp(:),nameV,nameG)
      endif
      
   end select

end subroutine sldup_restar
