subroutine sldale_turnon(a)
   !-----------------------------------------------------------------------
   !> This routine sets the mesh to ALE. The displacement and velocity pointers
   !! of the mesh are also initialized. 
   !-----------------------------------------------------------------------
   use typre  
   use Mod_sldAlemov
   implicit none

   class(sldAlemovProblem)  :: a
   integer(ip):: npoin,ipoin,idime,ndime,nelem,ielem,nboun,auxsize

   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   !We trick the module to believe its normal
   call a%solid%SetMesh(a%Mesh)

   call a%solid%SetParallelLibrary(a%ParallelLibrary)
   call a%solid%Mesh%SetParallelLibrary(a%ParallelLibrary)

   a%solid%LinearSystem => a%LinearSystem

   !set up solid as linear
   a%solid%sld_model='LINEA'
   a%solid%sld_type ='LINEA'
   a%solid%kfl_linea= 0

   a%solid%kfl_printPrincipalStresses=a%kfl_printPrincipalStresses
   a%solid%kfl_isSlave     = .true.
   a%solid%kfl_printStress = .true.
   a%solid%kfl_printStrain = .true.

   !set up DOF!
   call a%solid%SetNdofn
   call a%solid%SetNdofbc

   a%solid%kfl_tsche_2nd_datafile= 'NONE '
   a%solid%kfl_tsche_2nd_current = 'NONE '
   call a%solid%memall

   call a%Mesh%GetNpoin(npoin)

   !allocate basic vectors
   call a%solid%Memor%alloc(nboun,a%solid%kfl_bours,'kfl_bours','sldale_turnon')
   call a%solid%Memor%alloc(npoin,a%solid%kfl_fixrs,'kfl_fixrs','sldale_turnon')
   call a%solid%Memor%alloc(nboun,a%solid%kfl_fixbo,'kfl_fixbo','sldale_turnon')
   call a%solid%Memor%alloc(a%solid%ndofbc,npoin,a%solid%kfl_fixno,'kfl_fixno','sldale_turnon')
   call a%solid%Memor%alloc(a%solid%ndofbc,npoin,2,a%solid%bvess,'bvess','sldale_turnon')
   call a%solid%Memor%alloc(a%solid%ndofn,npoin,a%solid%unkno,'unkno','sldale_turnon')

   auxsize=size(a%bvnat,1)
   call a%solid%Memor%alloc(auxsize,a%solid%bvnat                  ,'bvnat'    ,'sldale_turnon')
   call a%solid%Memor%alloc(npoin  ,a%solid%kfl_funno              ,'kfl_funno','sldale_turnon')
   call a%solid%Memor%alloc(nboun  ,a%solid%kfl_funbo              ,'kfl_funbo','sldale_turnon')
   call a%solid%Memor%alloc(10     ,2            ,a%solid%kfl_funty,'kfl_funty','sldale_turnon')
   call a%solid%Memor%alloc(10     ,a%solid%funpa                  ,'funpa'    ,'sldale_turnon')

   call a%solid%SpecificTurnon

end subroutine sldale_turnon
