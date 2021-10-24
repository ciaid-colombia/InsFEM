subroutine ale_refine(a,itask)
   use typre
   use Mod_phpRefineArrays
   use Mod_Alemov
   implicit none
   class(AlemovProblem) :: a
   character(6) :: itask
   integer(ip), allocatable :: auxkfl_fixno(:,:)
   
   integer(ip) :: nelem,npoin,ndime,oldnpoin
   
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   
   call php_RefineArrays(a,itask,a%displacement,'Displacement')
   call php_RefineArrays(a,itask,a%velocity,'Velocity') 

   !kfl_Fixno
   oldnpoin = size(a%kfl_fixno0,2)
   call a%Memor%alloc(a%ndofbc,npoin,auxkfl_fixno,'kfl_fixno0','php_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(a%ndofbc,a%kfl_fixno0,auxkfl_fixno,'minim')
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(a%ndofbc,a%kfl_fixno0,auxkfl_fixno)
   endif
   call move_alloc(auxkfl_fixno,a%kfl_fixno0)
   call a%Memor%deallocObj(0,'kfl_fixno0','php_Refine',ip*a%ndofbc*oldnpoin)
   
   !Setting the displacements to the mesh
   call a%Mesh%SetALE(1_ip)
   call a%Mesh%SetDisplacements(a%Displacement)
   call a%Mesh%SetVelocities(a%Velocity)
   
   call a%ALESpecificRefine(itask)
   
end subroutine
