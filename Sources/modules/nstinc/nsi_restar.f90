subroutine nsi_restar(a,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file (itask=1)
   !    and writes the restart file when required.
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip), intent(in)    :: itask
   integer(ip)           :: npoin,ndime,idime,icomp,oldnpoin,npoinLocal,oldnpoinLocal,ncomp
   real(rp), allocatable :: press(:,:),auxpress(:,:),auxveloc(:,:)
   real(rp), allocatable :: restar_aux(:,:)
   character(150)        :: nameV,nameG
   integer(ip)           :: auxcoupling(1)=0

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)

   select case (itask)

   case (1)

      if (a%kfl_inter == 1) then
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(ndime,oldnpoin,restar_aux,'restar_aux','nsi_restar')
      call a%Memor%alloc(ndime,oldnpoin,auxveloc,'auxveloc','nsi_InterpolatedRestart')
      call a%Memor%alloc(1_ip,oldnpoin,auxpress,'auxpress','nsi_InterpolatedRestart')
      call a%Memor%alloc(1_ip,npoin,press,'press','nsi_InterpolatedRestart')

      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif

      do icomp = 3,ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         do idime=1, ndime
            write(nameV,"(A9,I1,A1,I1)") 'Velocity_',idime,'_',icomp-2
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxveloc(idime,1:oldnpoinLocal),nameV,nameG)
         enddo

         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
            call a%Int_Restart%Interpolate(ndime,auxveloc,a%veloc(:,:,icomp))         
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxveloc)
            a%veloc(:,:,icomp) = auxveloc
         endif
      end do
      do icomp = a%oldncomp+1,a%ncomp
         a%veloc(:,:,icomp) = a%veloc(:,:,a%oldncomp)
      end do

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Pressure'
      call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxpress(1_ip,1:oldnpoinLocal),nameV,nameG)

      if (a%kfl_inter == 1) then
         call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxpress(1,:))
         call a%Int_Restart%Interpolate(1_ip,auxpress,press)
      else
         call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxpress(1,:))
         press = auxpress
      endif
              
      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Coupling'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

      if (auxcoupling(1)==1 ) then
         if (a%kfl_docoupconv)  then

            restar_aux=0.0_rp
            do idime= 1,ndime
               write(nameV,"(A9,I1)") 'Veloc_cp_',idime
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(idime,:),nameV,nameG)
            end do
            if (a%kfl_inter == 1) then
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
               call a%Int_Restart%Interpolate(ndime,restar_aux,a%veloc_cp)         
            else
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,restar_aux)
               a%veloc_cp = restar_aux
            endif

            restar_aux=0.0_rp
            write(nameV,"(A8)") 'Press_cp'
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(1,:),nameV,nameG)
            if (a%kfl_inter == 1) then
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,restar_aux(1,:))
               call a%Int_Restart%Interpolate(1,restar_aux(1,:),a%press_cp(:))         
            else
               call a%Mesh%ArrayCommunicator%GhostCommunicate(1,restar_aux(1,:))
               a%press_cp(:) = restar_aux(1,:)
            endif

         endif
      endif

      a%press(:,3) = press(1,:)

      a%veloc(:,:,1) = a%veloc(:,:,3)
      a%veloc(:,:,2) = a%veloc(:,:,3)
      a%press(:,1) = a%press(:,3)
      a%press(:,2) = a%press(:,3)

      !if (a%kfl_trasg >= 1) then            
      !    call a%GaussArrayRestart(a%vesgs,1)
      !    if(a%kfl_repro == 2 .or. a%kfl_repro == 3) call a%GaussArrayRestart(a%vesgs2,1)
      !endif

      call a%Memor%dealloc(ndime,oldnpoin,auxveloc,'auxveloc','nsi_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,oldnpoin,auxpress,'auxpress','nsi_InterpolatedRestart')
      call a%Memor%dealloc(1_ip,npoin,press,'press','nsi_InterpolatedRestart')
      call a%Memor%dealloc(ndime,oldnpoin,restar_aux,'restar_aux','nsi_restar')

   case (2)
      
      do icomp = 3,a%ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         do idime= 1,ndime
            write(nameV,"(A9,I1,A1,I1)") 'Velocity_',idime,'_',icomp-2
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%veloc(idime,:,icomp),nameV,nameG)
         end do
      enddo

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Pressure'
      call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%press(:,3),nameV,nameG)

      if (a%kfl_docoupconv) then

         write(nameG,"(A10,I1)") 'Time Step ',1
         auxcoupling=1
         write(nameV,"(A8)") 'Coupling'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

         do idime= 1,ndime
             write(nameV,"(A9,I1)") 'Veloc_cp_',idime
             call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%veloc_cp(idime,:),nameV,nameG)
         end do

         write(nameV,"(A8)") 'Press_cp'
         call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%press_cp(:),nameV,nameG)
      else

         write(nameG,"(A10,I1)") 'Time Step ',1
         auxcoupling=0
         write(nameV,"(A8)") 'Coupling'
         call a%Writerpr%WriteAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)
      endif

   end select

   call a%NstincSpecificRestart(itask)
   call a%GaussRestart(itask)

end subroutine nsi_restar

subroutine nsi_gaussRestart(a,itask)
   use typre
   use def_parame
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip), intent(in)    :: itask
   character(150)             :: nameV, nameG 
   integer(ip)                :: npoin,ndime,idime,icomp,oldnpoin
   integer(ip)                :: npoinLocal,oldnpoinLocal,ncomp
   real(rp), allocatable      :: auxrepro(:,:)

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)

   select case (itask)

   case (1)

      if (a%kfl_inter == 1) then
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif


      if(a%kfl_repro>=1) then
         call a%Memor%alloc(a%ResidualSize,oldnpoin,auxrepro,'auxrepro','nsi_InterpolatedRestart')
         do idime = 1,a%ResidualSize
            write(nameG,"(A5)") 'Repro'
            write(nameV,"(A6,I1)") 'Repro_',idime
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxrepro(idime,:),nameV,nameG)
         end do

         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(a%ResidualSize,auxrepro)
            call a%Int_Restart%Interpolate(a%ResidualSize,auxrepro,a%repro)
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ResidualSize,auxrepro)
            a%repro = auxrepro
         endif
         call a%Memor%dealloc(a%ResidualSize,oldnpoin,auxrepro,'auxrepro','nsi_InterpolatedRestart')
      endif

   case (2)

      if(a%kfl_repro>=1) then
         do idime = 1,a%ResidualSize
            write(nameG,"(A5)") 'Repro'
            write(nameV,"(A6,I1)") 'Repro_',idime
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%repro(idime,:),nameV,nameG)
         end do
      endif

      if (a%kfl_trasg >= 1 .and. a%kfl_iofor == 0) then            
         call a%GaussArrayRestart(a%vesgs,2)
         if(a%kfl_repro == 2 .or. a%kfl_repro == 3) call a%GaussArrayRestart(a%vesgs2,2)
      endif

  end select

end subroutine nsi_gaussRestart
