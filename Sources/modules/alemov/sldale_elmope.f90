subroutine sldale_elmope(a,externalcurrentbvess)
   !-----------------------------------------------------------------------
   !> This routine performs an element loop and computes the elemental matrix 
   !! and RHS for the mesh displacements, applying boundary conditions. The free
   !! surface condition (kfl_fixno=3) is first changed to 1 so it can enter the 
   !! php_elmdir, where the Dirichlet boundary conditions are applied. Later on, 
   !! they are reset to 3.
   !-----------------------------------------------------------------------
   use typre
   use Mod_sldAlemov
   use Mod_Solids
   
   implicit none
   class(sldAlemovProblem), target :: a
   class(SolidsProblem), pointer :: sld => NULL()
   integer(ip) :: externalcurrentbvess,currentbvess,npoin,ipoin
   integer(ip) :: kfl_HangingNodes

   call a%Mesh%GetNpoin(npoin)

   sld => a%solid

   !For a purely Lagrangian Free Surface
   sld%kfl_fixno=a%kfl_fixno
   do ipoin=1,npoin                              ! 4 is used in FSI where disp are prescribed
       do currentbvess=1,a%ndofbc
           if (a%kfl_fixno(currentbvess,ipoin)==3 .or. a%kfl_fixno(currentbvess,ipoin)==4) then 
               sld%kfl_fixno(currentbvess,ipoin)=1
           end if
       end do
   end do

   sld%bvess = a%bvess

   !First solve system
   call sld%LinearSystem%ToZero
   call sld%Elmope
   !call sld%Bouope
   !Hanging nodes
   !call sld%Mesh%GetHanging(kfl_HangingNodes)
   !if (kfl_HangingNodes == 1) call sld%Mesh%AssemblyHangingNodesDiag(sld%ndofn,sld%LinearSystem,sld%Memor)
   call sld%LinearSystem%Solve(a%solid%unkno)
   call sld%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%solid%unkno)
   !HangingNodes
   !call sld%Mesh%GetHanging(kfl_HangingNodes)
   !if (kfl_HangingNodes == 1) call sld%Mesh%InterpolateHangingValues(sld%ndofn,sld%unkno)
   call sld%SpecificEndite(0)
   call sld%SpecificEndite(1_ip)
   a%solid%unkno=0.0_rp !just in case we cleanup

   sld%kfl_saveStrainsALE = .true.
   !Calculate stress and strains
   call sld%EndElmope('Endste')
   sld%kfl_saveStrainsALE = .false.

   !Second allow module to modify properties
   sld%kfl_modProperties = .true.

   !third solve system again for smooth strain field
   call sld%LinearSystem%ToZero
   call sld%Elmope
   !call sld%Bouope
   !Hanging nodes
   !call sld%Mesh%GetHanging(kfl_HangingNodes)
   !if (kfl_HangingNodes == 1) call sld%Mesh%AssemblyHangingNodesDiag(sld%ndofn,sld%LinearSystem,sld%Memor)
   call sld%LinearSystem%Solve(a%solid%unkno)
   call sld%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%solid%unkno)
   !HangingNodes
   !call sld%Mesh%GetHanging(kfl_HangingNodes)
   !if (kfl_HangingNodes == 1) call sld%Mesh%InterpolateHangingValues(sld%ndofn,sld%unkno)
   call sld%SpecificEndite(0)
   call sld%SpecificEndite(1_ip)

   !Calculate NEW stress and strains
   !Just for visualization purposes
   call sld%EndElmope('Endste')

   !reset flag for next step
   sld%kfl_modProperties = .false.

end subroutine
