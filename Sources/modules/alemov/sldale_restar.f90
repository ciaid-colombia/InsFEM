subroutine sldale_restar(a,itask)
   use typre
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem) :: a
   integer(ip)             :: itask
	real(rp), allocatable   :: restar_aux(:,:)
   integer(ip)    :: npoin,oldnpoin,idime,ndime
   integer(ip)    :: npoinLocal,oldnpoinLocal
   character(150) :: nameV,nameG
   real(rp)       :: auxarray(1)

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

      call a%Memor%alloc(ndime,oldnpoin,restar_aux,'restar_aux','ale_restar')

      a%solid%disp(:,:,1) = a%Displacement(:,:,1)

      write(nameG,"(A10,I1)") 'Time Step ',1

      write(nameV,"(A8)") 'ALErelax'
      auxarray = a%sldale_dispRelax
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxarray,nameV,nameG)

      restar_aux = 0.0_rp
      do idime=1, ndime
         write(nameV,"(A8,I1)") 'ALEbres_',idime
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(idime,:),nameV,nameG)
      end do
      if (a%kfl_inter == 1) then
         call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
         call a%Int_Restart%Interpolate(ndime,restar_aux,a%bres(:,:,1))
      else
         call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
         a%bres(:,:,1) = restar_aux(:,:)
      endif

      restar_aux = 0.0_rp
      do idime=1, ndime
         write(nameV,"(A9,I1)") 'ALEbvess_',idime
         call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(idime,:),nameV,nameG)
      end do
      if (a%kfl_inter == 1) then
         call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
         call a%Int_Restart%Interpolate(ndime,restar_aux,a%bvess(:,:,1))
      else
         call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
         a%bvess(:,:,1) = restar_aux(:,:)
      endif

      call a%Memor%dealloc(ndime,oldnpoin,restar_aux,'restar_aux','ale_restar')

  case (2)

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'ALErelax'
      auxarray = a%sldale_dispRelax
      call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxarray,nameV,nameG)

      do idime=1, ndime
         write(nameV,"(A8,I1)") 'ALEbres_',idime
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%bres(idime,:,1),nameV,nameG)
      end do

      do idime=1, ndime
         write(nameV,"(A9,I1)") 'ALEbvess_',idime
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%bvess(idime,:,1),nameV,nameG)
      end do

  end select

end subroutine
