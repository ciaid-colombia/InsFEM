subroutine geovtk(a,Mesh)
   !This routine dumps the geometry in VTK format
   use typre
   use def_parame
   use Mod_Memor
   use Mod_Listen
   use Mod_Element
   use Mod_Mesh
   use Mod_int2str
   use Mod_postpr_VTK
   use Mod_vtkPost
   implicit none
   
   class(PostprFile_VTK), intent(inout) :: a
   class(FemMesh)   , intent(inout) :: Mesh
   class(FiniteElement) , pointer     :: e => NULL()
   
   integer(ip)    :: ielty,ielty2,idime,ielem,inode,ipoin,ispos,aux_lnods(27),elemi,poini,nblock
   character(13)  :: elemt
   character(150) :: dumml
   real(rp)       :: rz,gpblock(1)
   real(rp), allocatable :: elemblocks(:)
   
   integer(ip), pointer :: nnode(:) => NULL(),ngaus(:) => NULL(),ltopo(:) => NULL(),lnode(:) => NULL()
   integer(ip)          :: ndime,npoin,nelty,nelem,pnode,npoin0,minpoin,npoinLocal,ghostLevel
   
   real(rp), pointer    :: coord(:) => NULL()


#ifdef VTK
   call Mesh%GetNpoin(npoin)
   call Mesh%GetNpoinLocal(npoinLocal)
   call Mesh%GetNdime(ndime)
   call Mesh%GetNelem(nelem)
   call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)

   if (npoin>0) then
      call Mesh%Local2Global(1,npoin0)
   else
      npoin0 = 0
   endif

   !Point data
   call VTK_BeginPoints(a%VTKInterface,npoin)
   rz = 0.0_rp
   do ipoin=1,npoin
      call Mesh%GetPointCoord(ipoin,coord)
      if (ndime == 3) rz = coord(3)
      call Mesh%Local2Global(ipoin,poini)
      !ipoin-1 vtk counter starts in 0
      call VTK_WritePoints(a%VTKInterface,ipoin-1,coord(1),coord(2),rz)
   end do
   call VTK_EndPoints(a%VTKInterface)

   !Loop over element types
   do ielty=1,nelty
      select case(ndime)
         case (2)
         if(nnode(ielty)==3.or.nnode(ielty)==6.or.nnode(ielty)==7.or.nnode(ielty)==10.or.nnode(ielty)==15) then
            elemt='Triangle' 
            if(nnode(ielty) == 6) then
               call VTK_SetElementType(a%VTKInterface,VTK_QUADRATIC_TRIANGLE)
            else
               call VTK_SetElementType(a%VTKInterface,VTK_TRIANGLE)
            end if
         else
            elemt='Quadrilateral'
            if(nnode(ielty) == 9) then
               call VTK_SetElementType(a%VTKInterface,VTK_BIQUADRATIC_QUAD)
            else
               call VTK_SetElementType(a%VTKInterface,VTK_QUAD)
            end if
         end if
      case (3)
         if(nnode(ielty)==4.or.nnode(ielty)==10) then 
            elemt='Tetrahedra'
            call VTK_SetElementType(a%VTKInterface,VTK_TETRA)
         else if(nnode(ielty)==8.or.nnode(ielty)==20.or.nnode(ielty)==27) then 
            elemt='Hexahedra'
            call VTK_SetElementType(a%VTKInterface,VTK_HEXAHEDRON)
         else if(nnode(ielty)==6.or.nnode(ielty)==15) then 
            elemt='Wedge'
            call VTK_SetElementType(a%VTKInterface,VTK_WEDGE)
         end if
      end select

      !Connectivity
      call VTK_BeginMesh(a%VTKInterface,nelem)
      do ielem=1,nelem
         ghostLevel = 1
         call Mesh%GetIelty(ielem,ielty2)
         if(ielty2==ielty) then
            call Mesh%GetLnode(ielem,pnode,lnode)
            call Mesh%Local2Global(pnode,lnode,aux_lnods)
            call Mesh%ElementLocal2Global(ielem,elemi)

            minpoin = minval(aux_lnods(1:pnode))
            if (minpoin >= npoin0 .and. minpoin <= npoin0+npoinLocal-1) ghostLevel = 0
            call VTK_BeginElements(a%VTKInterface,pnode)
            do inode=1,pnode
               !ipoin-1 vtk counter starts in 0
               call VTK_WriteElements(a%VTKInterface,lnode(inode)-1)
            enddo
            call VTK_EndElements(a%VTKInterface,ielem-1,ghostLevel)
         endif
      end do

   end do

   call VTK_WriteMesh(a%VTKInterface)

   !----------------Multiblock------------------

   if(a%kfl_writeType== 1) then
      call Mesh%GetNumBlocks(nblock)
      call VTK_BeginBlocks(a%VTKInterface,nblock)
      call VTK_BeginScalarGP(a%VTKInterface,len_trim("blocks"),trim("blocks"))
      !call VTK_BeginVectorGP(a%VTKInterface,len_trim("blocks"),trim("blocks"))

      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','geovtk')
      call a%Memor%alloc(e%pnode,elemblocks,'elemblocks','geovtk')

      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)  
         call e%gather(1,elemblocks,a%Mesh%geoblock)
         call e%interpc(1,elemblocks,gpblock(1))
         call VTK_WriteScalarGP(a%VTKInterface,ielem-1,gpblock(1))
         !call VTK_WriteVectorGP(a%VTKInterface,ielem-1,elemblocks(1),elemblocks(2),elemblocks(3))
      enddo

      call VTK_EndScalarGP(a%VTKInterface)
      !call VTK_EndVectorGP(a%VTKInterface)
      call a%Memor%dealloc(e%pnode,elemblocks,'elemblocks','geovtk')
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','geovtk')
   endif

#endif  

end subroutine geovtk

