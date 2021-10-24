subroutine MeshCloseDataFiles(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   if (  ((a%kfl_ReadType == 0) .and. (a%MPIrank == a%MPIroot)) .or. (a%kfl_ReadType == 1) ) then
         
         call a%Listener%CloseIncludeFiles
         call CloseFileDom(a%lun_pdata)
  
   endif
end subroutine