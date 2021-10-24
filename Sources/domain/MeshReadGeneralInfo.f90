subroutine MeshOpenFilesAndReadGeneralInfo(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   interface
      subroutine countv(a)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine reastr(a)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
      end subroutine
      
      subroutine ReadatMPI(a)
         use typre
         use Mod_Mesh
         implicit none
         class(FemMesh) :: a
      end subroutine
   end interface
   
   if (a%MPIrank == a%MPIroot .or. a%kfl_Readtype == 1) then
      call OpenFileDom(a%InputFolder,a%namda,a%lun_pdata)
      call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
   endif
   
   if (a%MPIrank == a%MPIroot) then
      !Read general domaindata
      call countv(a)           ! Count variables
      call reastr(a)           ! Read strategy
   endif
   call ReadatMPI(a)
end subroutine



