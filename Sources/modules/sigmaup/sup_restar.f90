subroutine sup_restar(a,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file (itask=1)
   !    and writes the restart file when required.
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_ThreeField 
   implicit none
   class(ThreeFieldNSProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip) :: icomp   !component of the unknown vector
   integer(ip) :: ndime,npoin,idime,ipoin,ntens,itens 
   integer(ip) :: oldnpoin, oldnpoinLocal, ncomp, npoinLocal
   real(rp), allocatable :: press(:,:),auxpress(:,:),auxveloc(:,:),auxsigma(:,:)
   character(150) ::  nameV, nameG 
   real(rp), allocatable    :: restar_aux(:,:)
   integer(ip)              :: auxcoupling(1)=0
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   ntens=(ndime-1)*(ndime-1)+2
         
   select case (itask)
   
   case (1) 
   
      if (a%kfl_inter == 1) then 
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(ntens,oldnpoin,restar_aux,'restar_aux','sup_restar')
      call a%Memor%alloc(ntens,oldnpoin,auxsigma,'auxsigma','sup_InterpolatedRestart')
   
   
      if (a%ncomp <= a%oldncomp) then
         ncomp = a%ncomp
      else
         ncomp = a%oldncomp
      endif
 
      do icomp = 3,ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         do itens=1, ntens
            write(nameV,"(A9,I1,A1,I1)") 'Stress_',itens,'_',icomp-2
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxsigma(itens,1:oldnpoinLocal),nameV,nameG)
         end do

         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ntens,auxsigma)
            call a%Int_Restart%Interpolate(ntens,auxsigma,a%sigma(:,:,icomp))         
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ntens,auxsigma)
            a%sigma(:,:,icomp) = auxsigma
         endif
      end do 
         
      do icomp = a%oldncomp+1,a%ncomp
         a%sigma(:,:,icomp) = a%sigma(:,:,a%oldncomp)
      end do

      write(nameG,"(A10,I1)") 'Time Step ',1
      write(nameV,"(A8)") 'Coupling'
      call a%Readerpr%ReadAttribute(a%fil_rstar,a%lun_rstar,auxcoupling,nameV,nameG)

      if (auxcoupling(1)==1 ) then
         if (a%kfl_docoupconv)  then

            restar_aux=0.0_rp
            do idime= 1,ntens
               write(nameV,"(A9,I1)") 'Sigma_cp_',idime
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,restar_aux(idime,:),nameV,nameG)
            end do
            if (a%kfl_inter == 1) then
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ntens,restar_aux)
               call a%Int_Restart%Interpolate(ntens,restar_aux,a%sigma_cp)         
            else
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ntens,restar_aux)
               a%sigma_cp = restar_aux
            endif

         endif
      endif
      
      
      !Update Stress at previous iteration and previous time step (guesses)
      icomp = 3
      a%sigma(:,:,1) = a%sigma(:,:,icomp)
      a%sigma(:,:,2) = a%sigma(:,:,icomp)
      
      call a%Memor%dealloc(ntens,oldnpoin,auxsigma,'auxsigma','sup_InterpolatedRestart')
      call a%Memor%dealloc(ntens,oldnpoin,restar_aux,'restar_aux','sup_restar')

   case (2)

      do icomp = 3,a%ncomp
         write(nameG,"(A10,I1)") 'Time Step ',icomp-2
         do idime= 1,ntens
            write(nameV,"(A7,I1,A1,I1)") 'Stress_',idime,'_',icomp-2
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%sigma(idime,:,icomp),nameV,nameG)
         end do
      enddo
      
      if (a%kfl_docoupconv) then

         write(nameG,"(A10,I1)") 'Time Step ',1
         do idime= 1,ntens
            write(nameV,"(A9,I1)") 'Sigma_cp_',idime
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%sigma_cp(idime,:),nameV,nameG)
         end do
      endif

   end select

end subroutine sup_restar

subroutine sup_gaussRestart(a,itask)
   !------------------------------------------------------------------------
   !    This routine reads the initial values from the restart file (itask=1)
   !    and writes the restart file when required.
   !------------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_ThreeField 
   implicit none
   class(ThreeFieldNSProblem) :: a
   integer(ip), intent(in)    :: itask
   character(150) ::  nameV, nameG 
   integer(ip) :: ndime,npoin,idime,ntens
   integer(ip) :: oldnpoin,oldnpoinLocal,npoinLocal
   real(rp), allocatable :: auxrepro(:,:),auxreproSDiv(:,:), auxreproUGradU(:,:), auxreproGradP(:,:)

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   ntens=(ndime-1)*(ndime-1)+2

   select case (itask)

   case (1)

      if (a%kfl_inter == 1) then 
         call a%OldMesh%GetNpoin(oldnpoin)
         call a%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      if (a%kfl_repro>=1) then
      
         call a%Memor%alloc(ntens+ndime+1_ip,oldnpoin,auxrepro,'auxrepro','sup_InterpolatedRestart')
         do idime = 1,ntens+ndime+1_ip
            write(nameG,"(A5)") 'Repro'
            write(nameV,"(A6,I1)") 'Repro_',idime
            call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxrepro(idime,:),nameV,nameG)
         end do
         
         if (a%kfl_inter == 1) then
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(ntens+ndime+1_ip,auxrepro)
            call a%Int_Restart%Interpolate(ntens+ndime+1_ip,auxrepro,a%repro)
         else
            call a%Mesh%ArrayCommunicator%GhostCommunicate(ntens+ndime+1_ip,auxrepro)
            a%repro = auxrepro
         endif
         
         call a%Memor%dealloc(ntens+ndime+1_ip,oldnpoin,auxrepro,'auxrepro','sup_InterpolatedRestart')
         
         if (a%kfl_repro>=2) then
         
            call a%Memor%alloc(ndime,oldnpoin,auxreproSDiv,'auxreproSDiv','sup_InterpolatedRestart')
            call a%Memor%alloc(ndime,oldnpoin,auxreproUGradU,'auxreproUGradU','sup_InterpolatedRestart')
            call a%Memor%alloc(ndime,oldnpoin,auxreproGradP,'auxreproGradP','sup_InterpolatedRestart')
            
            do idime = 1,ndime
               write(nameG,"(A5)") 'ReproSDiv'
               write(nameV,"(A6,I1)") 'ReproSDiv_',idime
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxreproSDiv(idime,:),nameV,nameG)
               
               write(nameG,"(A5)") 'ReproUGradU'
               write(nameV,"(A6,I1)") 'ReproUGradU_',idime
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxreproUGradU(idime,:),nameV,nameG)
               
               write(nameG,"(A5)") 'ReproGradP'
               write(nameV,"(A6,I1)") 'ReproGradP_',idime
               call a%Readerpr%ReadArray(a%fil_rstar,a%lun_rstar,auxreproGradP(idime,:),nameV,nameG)
            end do
            
            if (a%kfl_inter == 1) then
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxreproSDiv)
               call a%Int_Restart%Interpolate(ndime,auxreproSDiv,a%reproSDiv)
               
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxreproUGradU)
               call a%Int_Restart%Interpolate(ndime,auxreproUGradU,a%reproUGradU)
               
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(ndime,auxreproGradP)
               call a%Int_Restart%Interpolate(ndime,auxreproGradP,a%reproGradP)
            else
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxreproSDiv)
               a%reproSDiv = auxreproSDiv
               
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxreproUGradU)
               a%reproUGradU = auxreproUGradU
               
               call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,auxreproGradP)
               a%reproGradP = auxreproGradP
            endif
            call a%Memor%dealloc(ndime,oldnpoin,auxreproSDiv,'auxreproSDiv','sup_InterpolatedRestart')
            call a%Memor%dealloc(ndime,oldnpoin,auxreproUGradU,'auxreproUGradU','sup_InterpolatedRestart')
            call a%Memor%dealloc(ndime,oldnpoin,auxreproGradP,'auxreproGradP','sup_InterpolatedRestart')
         end if
      endif

   case(2)


      if (a%kfl_repro>=1) then
      
         do idime = 1,ntens+ndime+1_ip
            write(nameG,"(A5)") 'Repro'
            write(nameV,"(A6,I1)") 'Repro_',idime
            call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%repro(idime,:),nameV,nameG)
         end do
         
         if(a%kfl_repro>=2) then
         
            do idime=1,ndime
               write(nameG,"(A5)") 'ReproSDiv'
               write(nameV,"(A6,I1)") 'ReproSDiv_',idime
               call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%reproSDiv(idime,:),nameV,nameG)
               
               write(nameG,"(A5)") 'ReproUGradU'
               write(nameV,"(A6,I1)") 'ReproUGradU_',idime
               call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%reproUGradU(idime,:),nameV,nameG)
               
               write(nameG,"(A5)") 'ReproGradP'
               write(nameV,"(A6,I1)") 'ReproGradP_',idime
               call a%Writerpr%WriteArray(a%fil_rstar,a%lun_rstar,a%reproGradP(idime,:),nameV,nameG)
            end do
         end if
      endif
      
      if (a%kfl_trasg >= 1 .and. a%kfl_iofor == 0) then            
         call a%GaussArrayRestart(a%vesgs,2)
      endif

  end select

end subroutine sup_gaussRestart
