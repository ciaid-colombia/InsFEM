subroutine lev_memall(a)
   use typre
  
   use Mod_Memor
   use Mod_MPIObject
   use Mod_PhysicalProblem 
   use Mod_Mesh
   use Mod_LevelSet
   use def_parame
   use Mod_TimeIntegrator
   use Mod_CutMesh
   implicit none
   class(LevelSetProblem) :: a
   
   integer(ip) :: ndime,npoin,nelem

   type(TimeIntegratorDt1) :: Integrator
   integer(ip)             :: nsteps
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   
   !Set the Number of components necessary for the arrays
   call Integrator%Init(a%kfl_tsche_1st_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)
   a%ncomp = 1 + nsteps
   
   !level set variable
   call a%Memor%alloc(npoin,a%ncomp,a%level,'level','lev_memall')
   !Auxiliar list of elements
   call a%Memor%alloc(nelem,a%ElementListByLayers, 'ElementListByLayers','lev_memall')
   
   !Cut element type: compute element type, point status, subelements...   
   allocate(a%CutMesh)
   call a%Memor%allocObj(0,'CutMesh','lev_memall',1_ip)
   call a%CutMesh%allocCutMesh(a%Memor,a%Mesh)
   
   a%kfl_SmoothGradient=1
   if (a%kfl_SmoothGradient == 1) then
      call a%Memor%alloc(ndime,npoin,a%SmoothGradient,'SmoothGradient','lev_memall')
   endif   
   
end subroutine
