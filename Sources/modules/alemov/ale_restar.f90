subroutine ale_restar(a,itask)
   use typre
   use Mod_Alemov
   implicit none

   class(AlemovProblem) :: a
   integer(ip), intent(in) :: itask
   integer(ip) :: npoin,oldnpoin,idime,ndime,compi,icomp
   integer(ip) :: ncomp,npoinLocal,oldnpoinLocal,lim0,aux
   real(rp), pointer :: coord(:,:) => NULL()
   character(150) :: nameV,nameG
   integer(ip)    :: auxcoupling(1)=0

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)

   select case (itask)

   case (1)

      if (a%kfl_inter == 1) then
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoin(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(ndime,oldnpoin,a%oldncomp,a%restar_displ,'restar_displ','ale_restar')
      call a%Memor%alloc(ndime,oldnpoin,a%ncomp   ,a%restar_veloc,'restar_veloc','ale_restar')

      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Coupling'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

     !ALE Displacements loop
      do icomp = 1,ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp
         do idime=1, ndime
            write(nameV,"(A9,I1,A1,I1)") 'ALEdispl_',idime,'_',icomp
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,a%restar_displ(idime,:,icomp),nameV,nameG)
         end do
     enddo

      if (auxcoupling(1)==1) then
         aux  = 1 !For FSI coupling value needs to be assigned to past index
         lim0 = 2 !Start from second index
      else
         aux  = 0
         lim0 = 1
      endif

      do icomp = lim0,ncomp
         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,a%restar_displ(:,:,icomp-aux))
            call a%Int_Restart%Interpolate(ndime,a%restar_displ(:,:,icomp-aux),a%Displacement(:,:,icomp))
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%restar_displ(:,:,icomp-aux))
            a%Displacement(:,:,icomp) = a%restar_displ(:,:,icomp-aux)
         endif
      enddo
     do icomp = a%oldncomp+1,a%ncomp
        a%Displacement(:,:,icomp) = a%restar_disp(:,:,a%oldncomp)
     end do


     !ALE Velocity loop
     do icomp = 1,ncomp
        write(nameG,"(A10,I1)") 'Time Step ',icomp
        do idime=1, ndime
           write(nameV,"(A9,I1,A1,I1)") 'ALEveloc_',idime,'_',icomp
           call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,a%restar_veloc(idime,:,icomp),nameV,nameG)
        end do

        if (a%kfl_inter == 1) then
           call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,a%restar_veloc(:,:,icomp))
           call a%Int_Restart%Interpolate(ndime,a%restar_veloc(:,:,icomp),a%Velocity(:,:,icomp))
        else
           call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%restar_veloc(:,:,icomp))
           a%Velocity(:,:,icomp) = a%restar_veloc(:,:,icomp)
        endif

     end do
     do icomp = a%oldncomp+1,a%ncomp
         a%Velocity(:,:,icomp) = a%restar_veloc(:,:,a%oldncomp)
     end do

     call a%Mesh%SetVelocities(a%Velocity)

  case (2)

      write(nameG,"(A10,I1)") 'Time Step ',1
      if (a%kfl_docoupconv) then
         auxcoupling=1
         write(nameV,"(A8)") 'Coupling'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)
      else
         auxcoupling=0
         write(nameV,"(A8)") 'Coupling'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)
      endif

      do icomp = 1,a%ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp
         do idime=1, ndime
            write(nameV,"(A9,I1,A1,I1)") 'ALEdispl_',idime,'_',icomp
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%Displacement(idime,:,icomp),nameV,nameG)
         enddo
      end do

      do icomp = 1,a%ncomp
      write(nameG,"(A10,I1)") 'Time Step ',icomp
         do idime=1, ndime
            write(nameV,"(A9,I1,A1,I1)") 'ALEveloc_',idime,'_',icomp
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%Velocity(idime,:,icomp),nameV,nameG)
         end do
      end do

  end select

  call a%ParticularRestart(itask)

end subroutine
