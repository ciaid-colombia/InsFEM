subroutine sldale_turnof(a)
   !-----------------------------------------------------------------------
   !> This routine sets the mesh to ALE. The displacement and velocity pointers
   !! of the mesh are also initialized. 
   !-----------------------------------------------------------------------
   use typre  
   use Mod_sldAlemov
   implicit none

   class(sldAlemovProblem)  :: a
   integer(ip):: npoin,ipoin,idime,ndime,ielem,nelem,nboun,auxsize

   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNdime(ndime)

   call a%solid%SpecificTurnof

   !deallocate basic vectors
   call a%solid%Memor%dealloc(a%solid%ndofbc,npoin,a%solid%kfl_fixno,'kfl_fixno','sldale_turnof')
   call a%solid%Memor%dealloc(nboun,a%solid%kfl_fixbo,'kfl_fixbo','sldale_turnof')
   call a%solid%Memor%dealloc(a%solid%ndofbc,npoin,2,a%solid%bvess,'bvess','sldale_turnof')
   call a%solid%Memor%dealloc(a%solid%ndofn,npoin,a%solid%unkno,'unkno','sldale_turnof')
   auxsize=size(a%bvnat,1)
   call a%solid%Memor%dealloc(auxsize,a%solid%bvnat                  ,'bvnat'    ,'sldale_turnof')
   call a%solid%Memor%dealloc(npoin  ,a%solid%kfl_funno              ,'kfl_funno','sldale_turnof')
   call a%solid%Memor%dealloc(nboun  ,a%solid%kfl_funbo              ,'kfl_funbo','sldale_turnof')
   call a%solid%Memor%dealloc(10     ,2            ,a%solid%kfl_funty,'kfl_funty','sldale_turnof')
   call a%solid%Memor%dealloc(10     ,a%solid%funpa                  ,'funpa'    ,'sldale_turnof')

end subroutine sldale_turnof
