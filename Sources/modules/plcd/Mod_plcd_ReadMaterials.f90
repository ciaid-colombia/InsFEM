   subroutine ReadMaterials_Root(Listener,NumberOfMaterials,Materials,AuxRead_MaterialTypeList,MPIcomm,MPIsize,MPIroot,MPIrank,ndime,LargeStrainsflag,Memor)
      use typre
      use Mod_Listen
      use Mod_Mesh
      use Mod_plcd
      use Mod_Memor
      use Mod_plcd_MaterialFactory
      !use Mod_plcd_NotDef
      implicit none
      type(ListenFile) :: Listener
      integer(ip) :: NumberOfMaterials
      type(MatArray), allocatable :: Materials(:)
      character(5), allocatable :: AuxRead_MaterialTypeList(:)
      integer(ip) :: MPIsize,MPIroot,MPIrank,ndime,MPIcomm,LargeStrainsflag
      type(MemoryMan) :: Memor
      
      integer(ip) :: idnum
      character(5) :: MaterialName
      
      if(Listener%words(1) == 'NUMBE' .and. Listener%words(3) == 'MATER') then
         NumberOfMaterials = Listener%param(3)
         allocate(Materials(NumberOfMaterials))
         call Memor%allocObj(0,'Materials','plcd_memall',1*NumberOfMaterials)
         allocate(AuxRead_MaterialTypeList(NumberOfMaterials))
         AuxRead_MaterialTypeList = 'NOTDE'
         call Memor%allocObj(0,'AuxRead_MaterialTypeList','plcd_memall',1*NumberOfMaterials)
      elseif (Listener%words(1) == 'MATER') then
         do while(Listener%words(1)/='ENDMA')
            call Listener%Listen('plcd_reaphy')
            
            if (Listener%words(1) == 'BEGIN') then
               idnum = -1
               do while(Listener%words(1)/='ENDMA')
                  call Listener%Listen('plcd_reaphy')
                  
                  if(Listener%words(1) == 'IDNUM') then
                     idnum = Listener%param(1)
                     if(idnum > NumberOfMaterials) call runend('The ID material number can not be greater than Number of Materials')
                     
                  elseif(Listener%words(1) == 'MATNA') then
                    MaterialName = Listener%words(2)
            
                  elseif (Listener%words(1) == 'TYPEL') then
                     AuxRead_MaterialTypeList(idnum) = Listener%words(2)
                     call AllocatePLCDMaterial(AuxRead_MaterialTypeList(idnum),Materials(idnum)%p)
                     call Materials(idnum)%p%SetMPI(MPIcomm,MPIsize,MPIroot,MPIrank)
                     call Materials(idnum)%p%SetNdime(ndime)
                     call Materials(idnum)%p%SetFlagLargeStrains(LargeStrainsflag)
                     call Materials(idnum)%p%ReadData(Listener)
                     call Materials(idnum)%p%SetIdNum(idnum)   !This is necessary for load rebalancing
                     call Materials(idnum)%p%SetMaterialName(MaterialName)   !This is necessary for load rebalancing
                  endif   
                  
               enddo              
            endif
 
         enddo
      
      endif
   
   end subroutine
   
   subroutine ReadMaterials_Scatter(NumberOfMaterials,Materials,AuxRead_MaterialTypeList,MPIcomm,MPIsize,MPIroot,MPIrank,ndime,LargeStrainsflag,Memor)
      use MPI
      use typre
      use Mod_Listen
      use Mod_Mesh
      use Mod_plcd
      use Mod_Memor
      use Mod_plcd_MaterialFactory
      !use Mod_plcd_NotDef
      implicit none
      type(ListenFile) :: Listener
      integer(ip) :: NumberOfMaterials
      type(MatArray), allocatable :: Materials(:)
      character(5), allocatable :: AuxRead_MaterialTypeList(:)
      integer(ip) :: MPIcomm,MPIsize,MPIroot,MPIrank,ndime,LargeStrainsflag
      type(MemoryMan) :: Memor
      
      integer(ip) :: imaterial,ierr
   
      call MPI_BCAST(NumberOfMaterials,1,MPI_INTEGER4,MPIroot,MPIcomm,ierr)
      if (MPIrank /= MPIroot) then
         allocate(Materials(NumberOfMaterials))
         call Memor%allocObj(0,'Materials','plcd_memall',1*NumberOfMaterials)
         allocate(AuxRead_MaterialTypeList(NumberOfMaterials))
         AuxRead_MaterialTypeList = 'NOTDE'
         call Memor%allocObj(0,'AuxRead_MaterialTypeList','plcd_memall',1*NumberOfMaterials)
      endif
      call MPI_BCAST(AuxRead_MaterialTypeList,len(AuxRead_MaterialTypeList)*size(AuxRead_MaterialTypeList),MPI_CHARACTER,MPIroot,MPIcomm,ierr)
      if (MPIrank /= MPIroot) then
         do imaterial = 1, NumberOfMaterials 
            
            call AllocatePLCDMaterial(AuxRead_MaterialTypeList(imaterial),Materials(imaterial)%p)
            call Materials(imaterial)%p%SetMPI(MPIcomm,MPIsize,MPIroot,MPIrank)
            call Materials(imaterial)%p%SetNdime(ndime)
            !call Materials(idnum)%p%SetIdNum(idnum)   !This is necessary for load rebalancing
            call Materials(imaterial)%p%SetFlagLargeStrains(LargeStrainsflag)
         end do
      endif
      !Material type list is no longer needed
      deallocate(AuxRead_MaterialTypeList)
      call Memor%deallocObj(0,'AuxRead_MaterialTypeList','plcd_turnof',1*NumberOfMaterials)
   
      !Scatter Material Data
      do imaterial = 1,NumberOfMaterials
         call Materials(imaterial)%p%ScatterData
      enddo
      
      
      
   end subroutine
