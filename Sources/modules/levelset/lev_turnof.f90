subroutine lev_turnof(a)
   use typre
   use Mod_LevelSet
   use def_parame
   implicit none
   class(LevelSetProblem) :: a

   integer(ip) :: ndime,npoin,nelem,pnode,pgaus,ielem,nboun,ncomp

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   ncomp = a%ncomp

   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outph,'(//,a,/)') '     * * * END OF LEVELSET RUN * * *'
      write(a%lun_solve,'(//,a,/)') '     * * * END OF LEVELSET RUN * * *'
   endif
   call a%Memor%dealloc(npoin,a%ncomp,a%level,'level','lev_turnof')
   !Auxiliar list of elements
   call a%Memor%dealloc(nelem,a%ElementListByLayers, 'ElementListByLayers','lev_turnof')

   !Deallocate the cut elements
   call a%CutMesh%deallocCutElement

   !Deallocate the cut mesh
   call a%CutMesh%deallocCutMesh
   call a%Memor%deallocObj(0,'CutMesh','lev_turnof',1_ip)

   if (a%kfl_SmoothGradient == 1) then
      call a%Memor%dealloc(ndime,npoin,a%SmoothGradient,'SmoothGradient','lev_turnof')
   endif

   !HeightGauges
   if (a%nHeightGauges > 0) then
      deallocate(a%HeightGauges)
      call a%Memor%deallocObj(0,'HeightGauges','lev_reaous',a%nHeightGauges)
   endif
end subroutine
