module Mod_rom_IO
   use MPI
   use typre
   use def_parame
   use Mod_Int2Str
   use Mod_Iofile
   use Mod_ReadWrite
   use Mod_PodRom
   use Mod_rom_Procedure
   implicit none
   private
   public SetPointersIO

   !External Procedures
   procedure() :: NULLSUB

   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersIO(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   

      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1

      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            select case (a%kfl_itask)

            case('FOM')
               ProcHook%Output => OpenNewBasis
               call ConcatenateProcedures(ProcHook%Output,OutputBuild)
               call ConcatenateProcedures(ProcHook%Output,WriteBasis)
               if (a%kfl_eigentype .ne. "EPS") then
                  call ConcatenateProcedures(ProcHook%Output,WriteMean)
                  !call ConcatenateProcedures(ProcHook%Output,WriteBasisR)
               end if
               call ConcatenateProcedures(ProcHook%Output,CloseBasis)
               if (a%kfl_SaveSnap) then
                  if (a%kfl_BkpOpen) then
                     call ConcatenateProcedures(ProcHook%OutputBkp,OpenOldSnp)
                     call ConcatenateProcedures(ProcHook%OutputBkp,InputSnapData)
                  else
                     call ConcatenateProcedures(ProcHook%OutputBkp,OpenNewSnp)
                  endif
                  call ConcatenateProcedures(ProcHook%OutputBkp,OutputSnap)
                  call ConcatenateProcedures(ProcHook%OutputBkp,CloseSnp)
               end if

            case('SNP')
               if (a%kfl_BkpOpen) then
                  call ConcatenateProcedures(ProcHook%OutputBkp,OpenOldSnp)
                  call ConcatenateProcedures(ProcHook%OutputBkp,InputSnapData)
               else
                  call ConcatenateProcedures(ProcHook%OutputBkp,OpenNewSnp)
                  call ConcatenateProcedures(ProcHook%OutputBkp,OpenNewMsh)
               endif
               call ConcatenateProcedures(ProcHook%OutputBkp,OutputSnap)
               call ConcatenateProcedures(ProcHook%OutputBkp,CloseSnp)

               call ConcatenateProcedures(ProcHook%OutputBkp,OutputMesh)
               call ConcatenateProcedures(ProcHook%OutputBkp,CloseMsh)

            case('BFS')
               ProcHook%Output => OpenNewBasis
               call ConcatenateProcedures(ProcHook%Output,OutputBuild)
               call ConcatenateProcedures(ProcHook%Output,WriteBasis)
               if (a%kfl_eigentype .ne. "EPS") call ConcatenateProcedures(ProcHook%Output,WriteMean)
               call ConcatenateProcedures(ProcHook%Output,CloseBasis)
               call ConcatenateProcedures(ProcHook%Input,OpenOldSnp)
               call ConcatenateProcedures(ProcHook%Input,InputSnapData)
               call ConcatenateProcedures(ProcHook%Input,InputSnapBkp)
               call ConcatenateProcedures(ProcHook%Input,CloseSnp)

            case('INT')
               ProcHook%InputD => OpenOldBasis
               call ConcatenateProcedures(ProcHook%InputD,InputData)
               call ConcatenateProcedures(ProcHook%InputD,InputSigma)
               call ConcatenateProcedures(ProcHook%InputD,InputParam)
               call ConcatenateProcedures(ProcHook%Input,InputBasis)
               call ConcatenateProcedures(ProcHook%Input,CloseBasis)
               call ConcatenateProcedures(ProcHook%Input,OrthogonalizeBasis)
               ProcHook%Output => OpenNewBasis
               call ConcatenateProcedures(ProcHook%Output,OutputBuild)
               call ConcatenateProcedures(ProcHook%Output,WriteBasis)
               call ConcatenateProcedures(ProcHook%Output,OrthogonalizeDeallocate)
               call ConcatenateProcedures(ProcHook%Output,WriteMean2)
               call ConcatenateProcedures(ProcHook%Output,CloseBasis)

            case('ROM')
               ProcHook%InputD => OpenOldBasis
               call ConcatenateProcedures(ProcHook%InputD,InputData)
               call ConcatenateProcedures(ProcHook%InputD,InputSigma)
               call ConcatenateProcedures(ProcHook%InputD,InputParam)
               call ConcatenateProcedures(ProcHook%Input,InputBasis)
               call ConcatenateProcedures(ProcHook%Input,CloseBasis)

            end select

         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine
   
   subroutine OpenNewBasis

      a%fil_basis = trim(a%RestartFolder)//'/'//adjustl(trim(a%namda))//'.rom.bas'
      if (a%MPIrank == a%MPIroot) call iofile(zero,a%lun_basis,a%fil_basis,'basis','replace') ! open file

   end subroutine

   subroutine OpenOldBasis

      a%fil_basis = trim(a%OldRestartFolder)//'/'//adjustl(trim(a%oldnamda))//'.rom.bas'
      if (a%MPIrank == a%MPIroot) call iofile(zero,a%lun_basis,a%fil_basis,'basis','old') ! open file

   end subroutine

   subroutine OpenNewSnp

      a%fil_snap = trim(a%RestartFolder)//'/'//adjustl(trim(a%namda))//'.rom.snp'
      if (a%kfl_BkpOpen .eqv. .false.) then
        if (a%MPIrank == a%MPIroot)  call iofile(zero,a%lun_snap,a%fil_snap,'snaps','replace') ! open file
         a%kfl_BkpOpen = .true.
      else
        if (a%MPIrank == a%MPIroot)  call iofile(zero,a%lun_snap,a%fil_snap,'snaps','old') ! open file
      end if

   end subroutine

   subroutine OpenOldSnp

      a%fil_snap = trim(a%OldRestartFolder)//'/'//adjustl(trim(a%oldnamda))//'.rom.snp'
      if (a%MPIrank == a%MPIroot)  call iofile(zero,a%lun_snap,a%fil_snap,'snaps','old') ! open file

   end subroutine

   subroutine OpenNewMsh

      a%fil_mesh = trim(a%OldRestartFolder)//'/'//adjustl(trim(a%oldnamda))//'.msh'
      if (a%MPIrank == a%MPIroot)  call iofile(zero,a%lun_mesh,a%fil_mesh,'mesh','replace') ! open file

   end subroutine

   subroutine CloseBasis

      if (a%MPIrank == a%MPIroot) call iofile(two,a%lun_basis,a%fil_basis,'basis')

   end subroutine

   subroutine CloseSnp

      if (a%MPIrank == a%MPIroot) call iofile(two,a%lun_snap,a%fil_snap,'snaps')

   end subroutine

   subroutine CloseMsh

      if (a%MPIrank == a%MPIroot) call iofile(two,a%lun_mesh,a%fil_mesh,'mesh')

   end subroutine

   subroutine OutputBuild
      character(150) :: nameV,nameG
      real(rp), allocatable :: sigma(:)
      integer(ip) :: aux(1)

      call a%EigenSystem%GetNConvergence(a%nconv)
      nameV = 'R'
      nameG = 'Basis'
      aux(1) = a%nconv
      call a%Writerpr%WriteAttribute(a%fil_basis,a%lun_basis,aux,nameV,nameG)
      
      nameV = 'P'
      nameG = 'SnapshotMean'
      aux(1) = 1
      call a%Writerpr%WriteAttribute(a%fil_basis,a%lun_basis,aux,nameV,nameG)
      
      call a%Memor%alloc(a%nconv,sigma,'sigma','rom_IO')
      call a%EigenSystem%GetSigma(sigma)
      nameV = 'EigenValues'
      nameG = 'Basis'
      call a%Writerpr%WriteData(a%fil_basis,a%lun_basis,sigma,nameV,nameG)
      call a%Memor%dealloc(a%nconv,sigma,'sigma','rom_IO')

   end subroutine

   subroutine WriteBasis
      integer(ip)              :: npoinLocal,idofr,idofn,j
      integer(ip), allocatable :: ispos(:)
      real(rp), allocatable    :: auxbas(:,:),auxbas2(:)
      character(150)           :: nameV,nameG
      character(3)             :: exmod

      exmod =adjustl(trim(a%exmod))
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Writing basis to magnetic tape, Begin'

      call a%Mesh%GetNpoinLocal(npoinLocal)

      call a%Memor%alloc(npoinLocal,ispos,'ispos','rom_IO')
      call a%Memor%alloc(npoinLocal*a%ndofn,a%nconv,auxbas,'auxBasis','rom_IO')
      call a%Memor%alloc(npoinLocal,auxbas2,'auxBasis2','rom_IO')

      call a%EigenSystem%GetBasis(auxbas)
      nameG = 'Basis'
      do idofn = 1, a%ndofn
         ispos = (/(j,j=idofn,a%ndofn*npoinLocal,a%ndofn)/)
         do idofr=1,a%nconv
            nameV = a%nameBasis(idofn)
            auxbas2 = auxbas(ispos,idofr)
            call a%Writerpr%WriteArray(a%fil_basis,a%lun_basis,auxbas2,nameV,nameG,idofr-1)
         end do
      enddo

      call a%Memor%dealloc(npoinLocal*a%ndofn,a%nconv,auxbas,'auxBasis','rom_IO')
      call a%Memor%dealloc(npoinLocal,auxbas2,'auxBasis2','rom_IO')
      call a%Memor%dealloc(npoinLocal,ispos,'ispos','rom_IO')

      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Writing basis to magnetic tape, End'

   end subroutine

   subroutine WriteBasisR
      integer(ip)              :: npoinLocal,idofr,idofn
      real(rp), allocatable    :: auxbas(:,:)
      character(150)           :: nameV,nameG
      character(3)             :: exmod

      exmod =adjustl(trim(a%exmod))
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Writing R-basis to magnetic tape, Begin'

      call a%Mesh%GetNpoinLocal(npoinLocal)

      call a%Memor%alloc(a%nconv,a%nconv,auxbas,'auxBasis','rom_IO')

      call a%EigenSystem%GetBasisR(auxbas)
      nameG = 'Basis'
      nameV = 'R-Basis'
      !call a%Writerpr%WriteAttribute(a%fil_basis,a%lun_basis,auxbas(1:a%ndofr,1:a%ndofr),nameV,nameG)

      call a%Memor%dealloc(npoinLocal*a%ndofn,a%nconv,auxbas,'auxBasis','rom_IO')

      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Writing R-basis to magnetic tape, End'

   end subroutine

   subroutine WriteMean
      integer(ip)              :: npoinLocal,idofn
      character(150)           :: nameV,nameG
      real(rp)                 :: param(1)

      call a%Mesh%GetNpoinLocal(npoinLocal)

      call a%EigenSystem%GetMean(a%SnapMean)
      nameG = 'SnapshotMean'
      do idofn = 1, a%ndofn
         nameV = a%nameBasis(idofn)
         call a%Writerpr%WriteArray(a%fil_basis,a%lun_basis,a%SnapMean(idofn,1:npoinLocal),nameV,nameG,0)
      enddo
      nameV = 'Parameters'
      param = 1
      call a%Writerpr%WriteData(a%fil_basis,a%lun_basis,param,nameV,nameG)

   end subroutine

   subroutine WriteMean2
      integer(ip)              :: npoinLocal,idofn
      character(150)           :: nameV,nameG
      real(rp)                 :: param(1)

      call a%Mesh%GetNpoinLocal(npoinLocal)
      
      nameG = 'SnapshotMean'
      do idofn = 1, a%ndofn
         nameV = a%nameBasis(idofn)
         call a%Writerpr%WriteArray(a%fil_basis,a%lun_basis,a%SnapMean(idofn,1:npoinLocal),nameV,nameG,0)
      enddo
      nameV = 'Parameters'
      param = 1
      call a%Writerpr%WriteData(a%fil_basis,a%lun_basis,param,nameV,nameG)

   end subroutine
   
   subroutine InputSnapData
      character(150) :: nameV,nameG
      integer(ip)    :: ierr,aux(1)

      nameV = 'NSnaps'
      nameG = 'Snapshots'
      call a%Readerpr%ReadAttribute(a%fil_snap,a%lun_snap,aux,nameV,nameG)
      a%wsnap = aux(1)
      call MPI_BCAST(a%wsnap,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)

      a%wsnap = floor(real(a%wsnap)/real(a%SnapshotsInterval))

   end subroutine

   subroutine OutputSnap
      integer(ip)    :: npoinLocal,idofn
      character(150) :: nameV,nameG
      integer(ip)    :: aux(1)

      call a%Mesh%GetNpoinLocal(npoinLocal)
      
      a%wsnap = a%wsnap+1
      nameV = 'NSnaps'
      nameG = 'Snapshots'
      aux = a%wsnap
      call a%Writerpr%WriteAttribute(a%fil_snap,a%lun_snap,aux,nameV,nameG)

      do idofn = 1, a%ndofn
         nameV = a%nameBasis(idofn)
         call a%Writerpr%WriteArray(a%fil_snap,a%lun_snap,a%Snapshots(idofn,1:npoinLocal),nameV,nameG,a%wsnap-1)
      end do

   end subroutine

   subroutine OutputMesh
      integer(ip)          :: npoinLocal,idime,ndime,ielem,nelem,pnode
      character(150)       :: nameV,nameG
      character (len=150),dimension(3):: coords
      real(rp), pointer    :: coord(:,:)
      integer(ip), pointer :: lnode(:)

      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNelem(nelem)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetCoord(coord)
      coords(1) = 'X'
      coords(2) = 'Y'
      coords(3) = 'Z'

      nameG = 'Coordinates'
      do idime = 1, ndime
         write(nameV,"(A1)") coords(idime)
         call a%Writerpr%WriteArray(a%fil_mesh,a%lun_mesh,coord(idime,1:npoinLocal),nameV,nameG)
      end do

      nameV = 'Mass'
      nameG = ''
      call a%Writerpr%WriteArray(a%fil_mesh,a%lun_mesh,a%Mesh%vmass(1:npoinLocal),nameV,nameG)
   end subroutine

   subroutine InputData
      character(150) :: nameV,nameG
      character(3)   :: exmod
      integer(ip)    :: ierr, aux(1)

      exmod =adjustl(trim(a%exmod))
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Reading input, Begin'

      nameV = 'P'
      nameG = 'SnapshotMean'
      call a%Readerpr%ReadAttribute(a%fil_basis,a%lun_basis,aux,nameV,nameG)
      a%npara = aux(1) 

      nameV = 'R'
      nameG = 'Basis'
      call a%Readerpr%ReadAttribute(a%fil_basis,a%lun_basis,aux,nameV,nameG)
      a%nconv = aux(1)
      a%ndofr = a%nconv
      
      call a%EigenSystem%SetNConvergence(a%nconv)

      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Reading input, End'
      
   end subroutine

   subroutine InputSigma
      character(150) :: nameV,nameG
      integer(ip)    :: ierr

      call a%Memor%alloc(a%nconv,a%sigma,'sigma','rom_IO')
      nameV = 'EigenValues'
      nameG = 'Basis'
      call a%Readerpr%ReadData(a%fil_basis,a%lun_basis,a%sigma,nameV,nameG)
      call a%EigenSystem%SetSigma(a%sigma)
      call a%BasisFilter
      call a%EigenSystem%SetBasisSize(a%ndofr)
      call a%Memor%dealloc(a%nconv,a%sigma,'sigma','rom_IO')

   end subroutine

   subroutine InputParam
      character(150)        :: nameV,nameG
      integer(ip)           :: ierr, imean,aux(1)
      real(rp), allocatable :: auxpa(:)
      logical, allocatable  :: auxmk(:)

      a%parfun(1,:) = 1
      a%parfun(2,:) = 0
      if (a%npara > 1) then
         call a%Memor%alloc(a%npara,auxpa,'param','rom_IO')
         call a%Memor%alloc(a%npara,auxmk,'param','rom_IO')
         nameV = 'Parameters'
         nameG = 'SnapshotMean'
         call a%Readerpr%ReadData(a%fil_basis,a%lun_basis,auxpa,nameV,nameG)
         auxmk = .true.
         auxpa(:) = abs(auxpa(:) - a%param)
         do imean = 1, 2
            aux = minloc(auxpa,auxmk)
            a%parfun(1,imean) = aux(1)
            a%parfun(2,imean) = minval(auxpa,auxmk)
            auxmk(int(a%parfun(1,imean))) = .false.
         end do
         a%parfun(2,:) = a%parfun(2,:)/sum(a%parfun(2,:))
         call a%Memor%dealloc(a%npara,auxpa,'param','rom_IO')
         call a%Memor%dealloc(a%npara,auxmk,'param','rom_IO')
         a%npara = 2
      end if

   end subroutine

   subroutine InputBasis
      real(rp), allocatable   :: auxBasis(:,:),auxSnapMeanT(:),auxSnapMean(:,:)
      character(150)          :: nameV,nameG
      integer(ip)             :: npoinLocal,npoin,ndime,idofr,idofn,ipdeg
      integer(ip)             :: oldnpoin,oldnpoinLocal,imean
      character(3)            :: exmod

      exmod =adjustl(trim(a%exmod))
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Reading basis, Begin'

      call a%EigenSystem%GetSigma(a%sigma)
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      if (a%kfl_inter == 1) then
         call a%Problem%OldMesh%GetNpoin(oldnpoin)
         call a%Problem%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif
      call a%Memor%alloc(oldnpoin,a%ndofr,auxBasis,'auxBasis','rom_IO')
      call a%Memor%alloc(oldnpoin,2,auxSnapMean,'auxSnapMean','rom_IO')
      call a%Memor%alloc(oldnpoin,auxSnapMeanT,'auxSnapMeanT','rom_IO')

      do idofn = 1, a%ndofn
         auxSnapMeanT = 0
         nameG = 'Basis'
         nameV = a%nameBasis(idofn)
         do idofr=1,a%ndofr
            call a%Readerpr%ReadArray(a%fil_basis,a%lun_basis,auxBasis(1:oldnpoinLocal,idofr),nameV,nameG,idofr-1)
         end do
         if (a%kfl_eigentype .ne. "EPS") then
            nameG = 'SnapshotMean'
            do imean = 1, a%npara
               call a%Readerpr%ReadArray(a%fil_basis,a%lun_basis,auxSnapMean(1:oldnpoinLocal,imean),nameV,nameG,int(a%parfun(1,imean)-1))
            enddo
         endif
         if (a%kfl_inter == 1) then
            do idofr = 1, a%ndofr
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxBasis(:,idofr))
               call a%Int_Restart%Interpolate(1,auxBasis(:,idofr),a%Basis(idofn,:,idofr))
            enddo
            do imean = 1, a%npara
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxSnapMean(:,imean))
               auxSnapMeanT = auxSnapMeanT + (1-a%parfun(2,imean))*auxSnapMean(:,imean)
            enddo
            call a%Int_Restart%Interpolate(1,auxSnapMeanT,a%SnapMean(idofn,:)) 

            if (a%kfl_updateBasis) then
               a%OldMeshBasis(idofn,:,:) = auxBasis
               a%OldMeshSnapMean(idofn,:) = auxSnapMeanT
            end if

         else
            do idofr = 1, a%ndofr
               call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxBasis(:,idofr))
               a%Basis(idofn,:,idofr) = auxBasis(:,idofr)
            enddo
            do imean = 1, a%npara
               call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxSnapMean(:,imean))
               auxSnapMeanT = auxSnapMeanT + (1-a%parfun(2,imean))*auxSnapMean(:,imean)
            enddo
            a%SnapMean(idofn,:) = auxSnapMeanT
         endif
      enddo

      call a%EigenSystem%SetMean(a%SnapMean(:,1:npoinLocal))
      do idofr=1,a%ndofr
         call a%EigenSystem%SetBasis(idofr,a%Basis(:,1:npoinLocal,idofr))
      end do

      call a%Memor%dealloc(oldnpoin,a%ndofr,auxBasis,'auxBasis','rom_IO')
      call a%Memor%dealloc(oldnpoin,2,auxSnapMean,'auxSnapMean','rom_IO')
      call a%Memor%dealloc(oldnpoin,auxSnapMeanT,'auxSnapMeanT','rom_IO')

      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Reading basis, End'
   end subroutine

   subroutine InputSnapBkp
      real(rp), allocatable   :: auxSnapshots(:)
      character(150)          :: nameV,nameG
      integer(ip)             :: npoinLocal,npoin,ndime,isnap,idofn
      integer(ip)             :: oldnpoin,oldnpoinLocal,isnapaux
      character(3)            :: exmod
      real(rp) :: times(4)


      exmod =adjustl(trim(a%exmod))
      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Reading snapshots, Begin'
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      if (a%kfl_inter == 1) then
         call a%Problem%OldMesh%GetNpoin(oldnpoin)
         call a%Problem%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(oldnpoin,auxSnapshots,'auxSnapshots','rom_IO')

      isnapaux = 1
      do isnap=1,a%wsnap
         a%isnap = isnap
         nameG = 'Snapshots'
         do idofn = 1, a%ndofn
            nameV = a%nameBasis(idofn)
            call a%Timer%ReadSnapshots%Tic
            call a%Readerpr%ReadArray(a%fil_snap,a%lun_snap,auxSnapshots(1:oldnpoinLocal),nameV,nameG,isnap-1)
            call a%Timer%ReadSnapshots%Toc
            call a%Timer%ReadSnapshots%GetValue(times(1))

            if (a%kfl_inter == 1) then
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxSnapshots)
               call a%Int_Restart%Interpolate(1,auxSnapshots,a%Snapshots(idofn,:))
            else
               call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxSnapshots)
               a%Snapshots(idofn,:) = auxSnapshots(:)
            endif
         end do

         call a%Timer%SetSnapshots%Tic
         call a%EigenSystem%SetSnapshots(a%Snapshots(:,1:npoinLocal))
         call a%Timer%SetSnapshots%Toc
         call a%Timer%SetSnapshots%GetValue(times(2))

         isnapaux = isnapaux + a%SnapshotsInterval

      enddo

      call a%Memor%dealloc(oldnpoin,auxSnapshots,'auxSnapshots','rom_IO')

      if(a%MPIrank==a%MPIroot) write(*,*) 'ROM '//adjustl(trim(exmod))//': Reading snapshots, End'

   end subroutine

   subroutine InputBasisSnaps
      real(rp), allocatable   :: auxBasis(:,:),auxSnapMean(:)
      real(rp), allocatable   :: auxBasis2(:,:,:),auxSnapMean2(:,:)
      character(150)          :: nameV,nameG
      integer(ip)             :: npoinLocal,npoin,ndime,idofr,idofn
      integer(ip)             :: oldnpoin,oldnpoinLocal
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      if (a%kfl_inter == 1) then
         call a%Problem%OldMesh%GetNpoin(oldnpoin)
         call a%Problem%OldMesh%GetNpoinLocal(oldnpoinLocal)
      else
         oldnpoin = npoin
         oldnpoinLocal = npoinLocal
      endif

      call a%Memor%alloc(oldnpoin,a%ndofr,auxBasis,'auxBasis','rom_IO')
      call a%Memor%alloc(oldnpoin,auxSnapMean,'auxSnapMean','rom_IO')
      call a%Memor%alloc(a%ndofn,oldnpoin,a%ndofr,auxBasis2,'auxBasis2','rom_IO')
      call a%Memor%alloc(a%ndofn,oldnpoin,auxSnapMean2,'auxSnapMean2','rom_IO')

      do idofn = 1, a%ndofn
         nameG = a%nameBasis(idofn)
         do idofr=1,a%ndofr
            write(nameV,"(I4)") idofr
            nameV = 'Basis '// adjustl(trim(namev))
            call a%Readerpr%ReadArray(a%fil_basis,a%lun_basis,auxBasis(1:oldnpoinLocal,idofr),nameV,nameG)
         end do
         nameG = 'SnapshotMean'
         nameV = a%nameBasis(idofn)
         call a%Readerpr%ReadArray(a%fil_basis,a%lun_basis,auxSnapMean(1:oldnpoinLocal),nameV,nameG)
         if (a%kfl_inter == 1) then
            do idofr = 1, a%ndofr
               call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxBasis(:,idofr))
               call a%Int_Restart%Interpolate(1,auxBasis(:,idofr),auxBasis2(idofn,:,idofr))
            enddo
            call a%OldMesh%ArrayCommunicator%GhostCommunicate(1,auxSnapMean)
            call a%Int_Restart%Interpolate(1,auxSnapMean,auxSnapMean2(idofn,:)) 
         else
            do idofr = 1, a%ndofr
               call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxBasis(:,idofr))
               auxBasis2(idofn,:,idofr) = auxBasis(:,idofr)
            enddo
            call a%Mesh%ArrayCommunicator%GhostCommunicate(1,auxSnapMean)
            auxSnapMean2(idofn,:) = auxSnapMean
         endif
      enddo
      a%SnapMean = a%SnapMean + auxSnapMean2
         
      do idofr = 1, a%ndofr
         a%inumBasis = a%inumBasis + 1
         a%isnap = a%inumBasis
         call a%EigenSystem%SetSnapshots(auxBasis2(:,1:npoinLocal,idofr))
      enddo

      call a%Memor%dealloc(oldnpoin,a%ndofr,auxBasis,'auxBasis','rom_IO')
      call a%Memor%dealloc(oldnpoin,auxSnapMean,'auxSnapMean','rom_IO')
      call a%Memor%dealloc(a%ndofn,oldnpoin,a%ndofr,auxBasis2,'auxBasis2','rom_IO')
      call a%Memor%dealloc(a%ndofn,oldnpoin,auxSnapMean2,'auxSnapMean2','rom_IO')

   end subroutine

   subroutine OrthogonalizeBasis
      call a%EigenSystem%OrthogonalizeBasis
   end subroutine

   subroutine OrthogonalizeDeallocate
      call a%EigenSystem%OrthoDeallocate
   end subroutine

end module

subroutine rom_Output(PodRom)
   use typre
   use def_parame
   use Mod_Int2Str
   use Mod_Iofile
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_IO
   implicit none
   class(PodRomProblem), target :: PodRom
   
   a => PodRom
  
   call a%Timer%Output%Tic
   call ProcHook%Output
   call a%Timer%Output%Toc

end subroutine

subroutine rom_InputData(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_IO
   implicit none
   class(PodRomProblem), target :: PodRom

   a => PodRom
   
   call ProcHook%InputD

end subroutine

subroutine rom_Input(PodRom)
   use typre
   use Mod_PodRom
   use Mod_rom_Procedure
   use Mod_rom_IO
   implicit none
   class(PodRomProblem), target :: PodRom

   a => PodRom
   
   call a%Timer%Input%Tic
   call ProcHook%Input
   call a%Timer%Input%Toc

end subroutine
