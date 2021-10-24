subroutine ale_reaMPI(a)
   use typre
   use Mod_BroadCastBuffer
   use Mod_ALEmov
   use Mod_sldAlemov
   implicit none
   class(AlemovProblem) :: a
   type(BroadCastBuffer) :: BBuffer
   
   !Initialize Buffer
   call BBuffer%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call BBuffer%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call BBuffer%Initialize(100,500)

   !Communicate nsi_reaphy
   call BBuffer%Add(a%kfl_isFixedMesh)
   call BBuffer%Add(a%kfl_RemeshingCriteria)
   call BBuffer%Add(a%kfl_doAitken)
   
   !If we are working with solid ALE and want to print strains
   select type (a)
   type is (sldAlemovProblem)
       call BBuffer%Add(a%kfl_printPrincipalStresses)
       call BBuffer%Add(a%sldale_dispRelax)
       call BBuffer%Add(a%sldale_dispRelax_max)
       call a%solid%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
       call a%solid%SpecificReaMPI
   end select
   
   !BroadCast and DeallocateBuffer
   call BBuffer%BroadCast
   call BBuffer%Dealloc

end subroutine
