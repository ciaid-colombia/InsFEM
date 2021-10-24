subroutine sld_reampi(a)
   use typre
   use Mod_BroadCastBuffer
   use Mod_Solids
   implicit none
   class(SolidsProblem)             :: a
   type(BroadCastBuffer),   target  :: BBuffer
   class(BroadCastBuffer) , pointer :: bb => NULL()
   integer(ip) :: ierr

   bb => BBuffer

   !Communicate variables to other cpu's

   call bb%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call bb%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call bb%Initialize(100,100)

   call bb%Add(a%beta)
   call bb%Add(a%omega)
   call bb%Add(a%densi)
   call bb%Add(a%young)
   call bb%Add(a%poisson)
   call bb%Add(a%lambda)
   call bb%Add(a%mu)
   call bb%Add(a%bulkMod)
   call bb%Add(a%grnor)
   call bb%Add(a%gravi)
   call bb%Add(a%traction)
   call bb%Add(a%kfl_traction)
   call bb%Add(a%kfl_outfm)
   call bb%Add(a%kfl_constPointForce)
   call bb%Add(a%kfl_PointForceTime)
   call bb%Add(a%force_factor)
   call bb%Add(a%delta_force)
   call bb%Add(a%sld_dispRelax)

   call bb%Add(a%kfl_eigenter)

   call bb%Add(5,a%sld_type)
   call bb%Add(5,a%sld_model)
   call bb%Add(5,a%sld_strain)

   !For correct printing of stress and strain
   call bb%Add(a%kfl_GaussStress)
   call bb%Add(a%kfl_GaussStrain)
   call bb%Add(a%kfl_NodalStress)
   call bb%Add(a%kfl_NodalStrain)
   call bb%Add(a%kfl_printStrain)
   call bb%Add(a%kfl_printStress)
   call bb%Add(a%kfl_printSigma)
   call bb%Add(a%kfl_printGaussSigma)
   call bb%Add(a%kfl_printNodalSigma)
   call bb%Add(a%kfl_printPress)
   call bb%Add(a%kfl_printGaussPress)
   call bb%Add(a%kfl_printNodalPress)

   call bb%BroadCast

   call bb%Dealloc

   call a%SolidSpecificReampi

end subroutine
