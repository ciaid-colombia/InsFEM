subroutine sldsup_reampi(a)
   use typre
   use Mod_BroadCastBuffer
   use Mod_SUPSolids
   implicit none
   class(SUPSolidsProblem)          :: a
   type(BroadCastBuffer),   target  :: BBuffer
   class(BroadCastBuffer) , pointer :: bb => NULL()
   integer(ip) :: ierr

   bb => BBuffer

   !Communicate variables to other cpu's

   call bb%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call bb%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call bb%Initialize(100,100)

   !-----SLDSUP stabilization constants-----
   call bb%Add(a%tau_u)
   call bb%Add(a%tau_s)
   call bb%Add(a%tau_p)
   call bb%Add(a%kfl_PrbCharL)
   call bb%Add(a%kfl_repro)
   call bb%Add(a%kfl_trasg)
   call bb%Add(a%kfl_nolsg)
   call bb%Add(a%kfl_tacsg)
   call bb%Add(a%kfl_printJ2Stresses)
   call bb%Add(a%kfl_printDevStrain)
   call bb%Add(a%kfl_printResiduals)
   call bb%Add(a%kfl_useSecantModulus)

   call bb%BroadCast

   call bb%Dealloc

end subroutine
