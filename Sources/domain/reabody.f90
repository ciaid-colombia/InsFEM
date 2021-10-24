subroutine reabody(a)
   use MPI
   use typre
   use Mod_Int2Str
   use Mod_Mesh
   use Mod_OneToAllBuffer   
   use Mod_OneToAllLoop 
   use Mod_OneToAllReadOnBodies   
   implicit none
   class(FemMesh) :: a
   
   character(150) :: outstr
   
   integer(ip), allocatable :: knodb(:),aux_knodb(:)
   
   integer(ip) :: kfl_Reading, kfl_doreceiver
   integer(ip) :: icount,counter
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2,mtag3 = 3
   integer            :: ierr
   integer            :: irequest01(a%MPIsize), irequest02(a%MPIsize), irequest03(a%MPIsize)
   integer            :: status(MPI_STATUS_SIZE)
   
   !Buffer
   integer(ip) :: SBMaxNSends,SBMaxSize, SBSafetySend
   type(OneToAllBuffer) :: Buffer
   !Aux Buffer
   integer(ip) :: RBufferNSends
   integer(ip) :: isend,nbody
   
   integer(ip) :: idummy,ioerr
   integer(ip) :: irank, sendFlag, SendFlag_Buffer(a%MPIsize)
   
   type(OneToAllReadOnBodies) :: OTAReadOnBodies
   
   
   interface
   
   end interface   
   
   !List of body allocation
   call a%Memor%alloc(a%nboun,a%lbody,'lbody','lbody')  
   
   !Reach the Bodies Section
   if (a%MPIrank == a%MPIroot) then
      !Reach the ExnorSystems section.
      call a%Listener%listen(outstr)
      
      !We count how many lines we read so that we can rewind them if there were no BODYD 
      !(only for old cases, new ones always have BODYD)
      !(in old cases it only reads 4 lines so it is not computationally costly)
      counter = 0
      do while(a%Listener%words(1)/='BODYD' .and. a%Listener%words(1)/='ENDGE')
         call a%Listener%listen(outstr)
         counter = counter +1
      end do
   endif

   CALL MPI_BCAST(a%Listener%words(1), 5, MPI_CHARACTER, a%MPIroot, a%MPIcomm, ierr)
   
   !If I arrived to the end of the geometry, the rewind just the counted ones
   if (a%Listener%words(1) == 'ENDGE') then
      if (a%MPIrank == a%MPIroot) then
         do icount = 1,counter
            backspace(a%Listener%nunit)
         enddo
      endif
      
   !Else do what we need to do   
   elseif (a%Listener%words(1)=='BODYD') then
         
      !Bodies
      call OTAReadOnBodies%SetType(a%kfl_ReadType)
      call OTAReadOnBodies%Initialize(a)
      call OTAReadOnBodies%Loop
      call OTAReadOnBodies%Finalize
      kfl_Reading = 0   
   endif

end subroutine

   
   
