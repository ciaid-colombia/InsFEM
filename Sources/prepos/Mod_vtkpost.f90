module Mod_vtkPost
   use   typre
   use, intrinsic :: ISO_C_Binding 
   implicit none
   integer, parameter :: VTK_EMPTY_CELL = 0, VTK_VERTEX = 1, VTK_POLY_VERTEX = 2, VTK_LINE = 3, VTK_POLY_LINE = 4, VTK_TRIANGLE = 5, VTK_TRIANGLE_STRIP = 6, VTK_POLYGON = 7, VTK_PIXEL = 8, VTK_QUAD = 9, VTK_TETRA = 10, VTK_VOXEL = 11, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13, VTK_PYRAMID = 14, VTK_PENTAGONAL_PRISM = 15, VTK_HEXAGONAL_PRISM = 16, VTK_QUADRATIC_EDGE = 21, VTK_QUADRATIC_TRIANGLE = 22, VTK_QUADRATIC_QUAD = 23, VTK_QUADRATIC_POLYGON = 36, VTK_QUADRATIC_TETRA = 24, VTK_QUADRATIC_HEXAHEDRON = 25, VTK_QUADRATIC_WEDGE = 26, VTK_QUADRATIC_PYRAMID = 27, VTK_BIQUADRATIC_QUAD = 28, VTK_TRIQUADRATIC_HEXAHEDRON = 29, VTK_QUADRATIC_LINEAR_QUAD = 30, VTK_QUADRATIC_LINEAR_WEDGE = 31, VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32, VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33, VTK_BIQUADRATIC_TRIANGLE = 34, VTK_CUBIC_LINE = 35, VTK_CONVEX_POINT_SET = 41, VTK_POLYHEDRON = 42, VTK_PARAMETRIC_CURVE = 51, VTK_PARAMETRIC_SURFACE = 52, VTK_PARAMETRIC_TRI_SURFACE = 53, VTK_PARAMETRIC_QUAD_SURFACE = 54, VTK_PARAMETRIC_TETRA_REGION = 55, VTK_PARAMETRIC_HEX_REGION = 56, VTK_HIGHER_ORDER_EDGE = 60, VTK_HIGHER_ORDER_TRIANGLE = 61, VTK_HIGHER_ORDER_QUAD = 62, VTK_HIGHER_ORDER_POLYGON = 63, VTK_HIGHER_ORDER_TETRAHEDRON = 64, VTK_HIGHER_ORDER_WEDGE = 65, VTK_HIGHER_ORDER_PYRAMID = 66, VTK_HIGHER_ORDER_HEXAHEDRON = 67
   integer, parameter :: Ascii = 0, Binary = 1, Appended = 2
   character(1,KIND=C_CHAR), PARAMETER :: C_CHAR_ARRAY(1) = (/ C_NULL_CHAR /)
   character(LEN=*), PARAMETER, PUBLIC :: VTK_NULL = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
   integer,                  PARAMETER :: VTK_MAX_CHAR_LEN   = 64
   integer(C_INT), PARAMETER :: MAX_STRING_LEN = 1024
   type CVTKInterface
      private
      type(C_ptr) :: object = C_NULL_ptr
   end type CVTKInterface

   public :: CVTKInterface
#ifdef VTK
   public :: VTK_Interface,VTK_Interface_des
   public :: VTK_BeginBlocks,VTK_EndBlocks,VTK_SetBlockName
   public :: VTK_BeginPoints,VTK_WritePoints,VTK_EndPoints
   public :: VTK_BeginMesh,VTK_SetElementType,VTK_BeginElements,VTK_WriteElements,VTK_EndElements
   public :: VTK_BeginWriter,VTK_WriteMesh,VTK_SetWriteFile,VTK_Flush

contains

   subroutine VTK_Interface(this)
      type(CVTKInterface), intent(out) :: this
      interface
         function I_VTK_Interface() result(this) bind(C,name="C_VTKInterface_new")
            import
            type(C_ptr) :: this
         end function I_VTK_Interface
      end interface
      this%object = I_VTK_Interface()
   end subroutine VTK_Interface
   subroutine VTK_Interface_des(this)
      type(CVTKInterface), intent(in) :: this
      interface
         subroutine I_VTK_Interface_des(this) bind(C,name="C_VTKInterface_des")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_Interface_des
      end interface
      call I_VTK_Interface_des(this%object)
   end subroutine VTK_Interface_des
   subroutine VTK_SetBlockName(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_SetBlockName(this,FORSTR_LEN,filename) bind(C,name="C_VTK_SetBlockName")
            import
            type(C_ptr),value                     :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_SetBlockName
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_SetBlockName(this%object,n,C_filename)
   end subroutine VTK_SetBlockName
   subroutine VTK_BeginBlocks(this,numB)
      type(CVTKInterface), intent(inout) :: this
      integer :: numB
      interface
         subroutine I_VTK_BeginBlocks(this,a) bind(C,name="C_VTK_BeginBlocks")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: a
         end subroutine I_VTK_BeginBlocks
      end interface
      call I_VTK_BeginBlocks(this%object,int(numB,C_int))
   end subroutine VTK_BeginBlocks
   subroutine VTK_EndBlocks(this)
      type(CVTKInterface), intent(inout) :: this
      interface
          subroutine I_VTK_EndBlocks(this) bind(C,name="C_VTK_EndBlocks")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndBlocks
      end interface
      call I_VTK_EndBlocks(this%object)
   end subroutine VTK_EndBlocks
   subroutine VTK_BeginPoints(this,numpoints)
      type(CVTKInterface), intent(inout) :: this
      integer :: numpoints
      interface
         subroutine I_VTK_BeginPoints(this,a) bind(C,name="C_VTK_BeginPoints")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: a
         end subroutine I_VTK_BeginPoints
      end interface
      call I_VTK_BeginPoints(this%object,int(numpoints,C_int))
   end subroutine VTK_BeginPoints
   subroutine VTK_WritePoints(this,id,x,y,z)
      type(CVTKInterface), intent(inout) :: this
      integer :: id
      real(rp)    :: x, y, z
      interface
         subroutine I_VTK_WritePoints(this,id,x,y,z) bind(C,name="C_VTK_WritePoints")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: id
            real(C_double), value :: x, y, z
         end subroutine I_VTK_WritePoints
      end interface
      call I_VTK_WritePoints(this%object,int(id,C_int),real(x,C_double),real(y,C_double),real(z,C_double))
   end subroutine VTK_WritePoints
   subroutine VTK_EndPoints(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndPoints(this) bind(C,name="C_VTK_EndPoints")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndPoints
      end interface
      call I_VTK_EndPoints(this%object)
   end subroutine VTK_EndPoints
   subroutine VTK_SetElementType(this,CellType)
      type(CVTKInterface), intent(inout) :: this
      integer :: CellType
      interface
         subroutine I_VTK_SetCellType(this,CellType) bind(C,name="C_VTK_SetCellType")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: CellType
         end subroutine I_VTK_SetCellType
      end interface
      call I_VTK_SetCellType(this%object,int(CellType,C_int))
   end subroutine VTK_SetElementType
   subroutine VTK_BeginMesh(this,numelem)
      type(CVTKInterface), intent(inout) :: this
      integer :: numelem
      interface
         subroutine I_VTK_BeginMesh(this,numelem) bind(C,name="C_VTK_BeginMesh")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: numelem
         end subroutine I_VTK_BeginMesh
      end interface
      call I_VTK_BeginMesh(this%object,int(numelem,C_int))
   end subroutine VTK_BeginMesh
   subroutine VTK_BeginElements(this,numpoints)
      type(CVTKInterface), intent(inout) :: this
      integer :: numpoints
      interface
         subroutine I_VTK_BeginCells(this,numnodes) bind(C,name="C_VTK_BeginCells")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: numnodes
         end subroutine I_VTK_BeginCells
      end interface
      call I_VTK_BeginCells(this%object,int(numpoints,C_int))
   end subroutine VTK_BeginElements
   subroutine VTK_WriteElements(this,globalids)
      type(CVTKInterface), intent(inout) :: this
      integer :: globalids
      interface
         subroutine I_VTK_WriteCells(this,globalids) bind(C,name="C_VTK_WriteCells")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: globalids
         end subroutine I_VTK_WriteCells
      end interface
      call I_VTK_WriteCells(this%object,int(globalids,C_int))
   end subroutine VTK_WriteElements
   subroutine VTK_EndElements(this,cellId,ghostLevel)
      type(CVTKInterface), intent(inout) :: this
      integer :: cellId,ghostLevel
      interface
         subroutine I_VTK_EndCells(this,cellId,ghostLevel) bind(C,name="C_VTK_EndCells")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: cellId,ghostLevel
         end subroutine I_VTK_EndCells
      end interface
      call I_VTK_EndCells(this%object,int(cellId,C_int),int(ghostLevel,C_int))
   end subroutine VTK_EndElements
   subroutine VTK_BeginWriter(this,VTKformat,wType)
      type(CVTKInterface), intent(inout) :: this
      integer :: VTKformat,wType
      interface
         subroutine I_VTK_BeginWriter(this,VTKformat,wType) bind(C,name="C_VTK_BeginWriter")
            import
            type(C_ptr),value :: this
            integer(C_int), value :: VTKformat,wType
         end subroutine I_VTK_BeginWriter
      end interface
      call I_VTK_BeginWriter(this%object,VTKformat,wType)
   end subroutine VTK_BeginWriter
   subroutine VTK_SetWriteFile(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_SetWriteFile(this,FORSTR_LEN,filename) bind(C,name="C_VTK_SetWriteFile")
            import
            type(C_ptr),value                     :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_SetWriteFile
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_SetWriteFile(this%object,n,C_filename)
   end subroutine VTK_SetWriteFile
   subroutine VTK_WriteMesh(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_WriteMesh(this) bind(C,name="C_VTK_WriteMesh")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_WriteMesh
      end interface
      call I_VTK_WriteMesh(this%object)
   end subroutine VTK_WriteMesh
   subroutine VTK_Flush(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_Flush(this) bind(C,name="C_VTK_Flush")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_Flush
      end interface
      call I_VTK_Flush(this%object)
   end subroutine VTK_Flush

   subroutine VTK_SetMPI(this,mpisize,mpirank)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: mpisize, mpirank
      interface
         subroutine I_VTK_SetMPI(this,mpisize,mpirank) bind(C,name="C_VTK_SetMPI")
            import
            type(C_ptr),value :: this
            integer(C_int),value :: mpisize, mpirank
         end subroutine I_VTK_SetMPI
      end interface
      call I_VTK_SetMPI(this%object,mpisize,mpirank)
   end subroutine VTK_SetMPI
   
   subroutine VTK_SetMulticomm(this,kfl_multicomm,MulticommColor)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: kfl_multicomm, MulticommColor
      interface
         subroutine I_VTK_SetMulticomm(this,kfl_multicomm,MulticommColor) bind(C,name="C_VTK_SetMulticomm")
            import
            type(C_ptr),value :: this
            integer(C_int),value :: kfl_multicomm, MulticommColor
         end subroutine I_VTK_SetMulticomm
      end interface
      call I_VTK_SetMulticomm(this%object,kfl_multicomm,MulticommColor)
   end subroutine VTK_SetMulticomm
   
   subroutine VTK_Reset(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_Reset(this) bind(C,name="C_VTK_Reset")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_Reset
      end interface
      call I_VTK_Reset(this%object)
   end subroutine VTK_Reset

   subroutine VTK_SetTimeWriter(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_SetTimeWriter(this,FORSTR_LEN,filename) bind(C,name="C_VTK_SetTimeWriter")
            import
            type(C_ptr),value                     :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_SetTimeWriter
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_SetTimeWriter(this%object,n,C_filename)
   end subroutine VTK_SetTimeWriter
   subroutine VTK_AddTime(this,FORSTR_LEN,F_filename,time,istep)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      real(rp), intent(in)               :: time
      integer(ip), intent(in)            :: istep
      interface
         subroutine I_VTK_AddTime(this,FORSTR_LEN,filename,time,istep) bind(C,name="C_VTK_AddTime")
            import
            type(C_ptr),value     :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
            real(C_double), value :: time
            integer(C_int), value :: istep
         end subroutine I_VTK_AddTime
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_AddTime(this%object,n,C_filename,time,istep)
   end subroutine VTK_AddTime
   subroutine VTK_EndTimeWriter(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_EndTimeWriter(this,FORSTR_LEN,filename) bind(C,name="C_VTK_EndTimeWriter")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_EndTimeWriter
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_EndTimeWriter(this%object,n,C_filename)
   end subroutine VTK_EndTimeWriter
   subroutine VTK_AddCycle(this,cycle)
      type(CVTKInterface), intent(inout) :: this
      integer(ip), intent(in)               :: cycle
      interface
         subroutine I_VTK_AddCycle(this,cycle) bind(C,name="C_VTK_AddCycle")
            import
            type(C_ptr),value     :: this
            integer(C_int), value :: cycle
         end subroutine I_VTK_AddCycle
      end interface
      call I_VTK_AddCycle(this%object,cycle)
   end subroutine VTK_AddCycle

   subroutine VTK_BeginScalarGP(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginScalarGP(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginScalarGP")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginScalarGP
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginScalarGP(this%object,n,C_filename)
   end subroutine VTK_BeginScalarGP
   subroutine VTK_WriteScalarGP(this,pointid,scalar)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: scalar
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteScalarGP(this,pointid,scalar) bind(C,name="C_VTK_WriteScalarGP")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: scalar
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteScalarGP
      end interface
      call I_VTK_WriteScalarGP(this%object,pointid,scalar)
   end subroutine VTK_WriteScalarGP
   subroutine VTK_EndScalarGP(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndScalarGP(this) bind(C,name="C_VTK_EndScalarGP")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndScalarGP
      end interface
      call I_VTK_EndScalarGP(this%object)
   end subroutine VTK_EndScalarGP

   subroutine VTK_BeginScalar(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginScalar(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginScalar")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginScalar
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginScalar(this%object,n,C_filename)
   end subroutine VTK_BeginScalar
   subroutine VTK_WriteScalar(this,pointid,scalar)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: scalar
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteScalar(this,pointid,scalar) bind(C,name="C_VTK_WriteScalar")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: scalar
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteScalar
      end interface
      call I_VTK_WriteScalar(this%object,pointid,scalar)
   end subroutine VTK_WriteScalar
   subroutine VTK_EndScalar(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndScalar(this) bind(C,name="C_VTK_EndScalar")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndScalar
      end interface
      call I_VTK_EndScalar(this%object)
   end subroutine VTK_EndScalar

   subroutine VTK_BeginVectorGP(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginVectorGP(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginVectorGP")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginVectorGP
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginVectorGP(this%object,n,C_filename)
   end subroutine VTK_BeginVectorGP
   subroutine VTK_WriteVectorGP(this,pointid,x,y,z)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x,y,z
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteVectorGP(this,pointid,x,y,z) bind(C,name="C_VTK_WriteVectorGP")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x,y,z
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteVectorGP
      end interface
      call I_VTK_WriteVectorGP(this%object,pointid,x,y,z)
   end subroutine VTK_WriteVectorGP
   subroutine VTK_EndVectorGP(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndVectorGP(this) bind(C,name="C_VTK_EndVectorGP")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndVectorGP
      end interface
      call I_VTK_EndVectorGP(this%object)
   end subroutine VTK_EndVectorGP

   subroutine VTK_BeginVector(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginVector(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginVector")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginVector
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginVector(this%object,n,C_filename)
   end subroutine VTK_BeginVector
   subroutine VTK_WriteVector(this,pointid,x,y,z)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x,y,z
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteVector(this,pointid,x,y,z) bind(C,name="C_VTK_WriteVector")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x,y,z
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteVector
      end interface
      call I_VTK_WriteVector(this%object,pointid,x,y,z)
   end subroutine VTK_WriteVector
   subroutine VTK_EndVector(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndVector(this) bind(C,name="C_VTK_EndVector")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndVector
      end interface
      call I_VTK_EndVector(this%object)
   end subroutine VTK_EndVector

   subroutine VTK_BeginArray4GP(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginArray4GP(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginArray4GP")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginArray4GP
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginArray4GP(this%object,n,C_filename)
   end subroutine VTK_BeginArray4GP
   subroutine VTK_WriteArray4GP(this,pointid,x,y,z,w)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x,y,z,w
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteArray4GP(this,pointid,x,y,z,w) bind(C,name="C_VTK_WriteArray4GP")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x,y,z,w
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteArray4GP
      end interface
      call I_VTK_WriteArray4GP(this%object,pointid,x,y,z,w)
   end subroutine VTK_WriteArray4GP
   subroutine VTK_EndArray4GP(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndArray4GP(this) bind(C,name="C_VTK_EndArray4GP")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndArray4GP
      end interface
      call I_VTK_EndArray4GP(this%object)
   end subroutine VTK_EndArray4GP

   subroutine VTK_BeginArray4(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginArray4(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginArray4")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginArray4
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginArray4(this%object,n,C_filename)
   end subroutine VTK_BeginArray4
   subroutine VTK_WriteArray4(this,pointid,x,y,z,w)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x,y,z,w
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteArray4(this,pointid,x,y,z,w) bind(C,name="C_VTK_WriteArray4")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x,y,z,w
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteArray4
      end interface
      call I_VTK_WriteArray4(this%object,pointid,x,y,z,w)
   end subroutine VTK_WriteArray4
   subroutine VTK_EndArray4(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndArray4(this) bind(C,name="C_VTK_EndArray4")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndArray4
      end interface
      call I_VTK_EndArray4(this%object)
   end subroutine VTK_EndArray4

   subroutine VTK_BeginArray6GP(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginArray6GP(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginArray6GP")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginArray6GP
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginArray6GP(this%object,n,C_filename)
   end subroutine VTK_BeginArray6GP
   subroutine VTK_WriteArray6GP(this,pointid,x,y,z,w,v,t)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x,y,z,w,v,t
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteArray6GP(this,pointid,x,y,z,w,v,t) bind(C,name="C_VTK_WriteArray6GP")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x,y,z,w,v,t
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteArray6GP
      end interface
      call I_VTK_WriteArray6GP(this%object,pointid,x,y,z,w,v,t)
   end subroutine VTK_WriteArray6GP
   subroutine VTK_EndArray6GP(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndArray6GP(this) bind(C,name="C_VTK_EndArray6GP")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndArray6GP
      end interface
      call I_VTK_EndArray6GP(this%object)
   end subroutine VTK_EndArray6GP

   subroutine VTK_BeginArray6(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginArray6(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginArray6")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginArray6
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginArray6(this%object,n,C_filename)
   end subroutine VTK_BeginArray6
   subroutine VTK_WriteArray6(this,pointid,x,y,z,w,v,t)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x,y,z,w,v,t
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteArray6(this,pointid,x,y,z,w,v,t) bind(C,name="C_VTK_WriteArray6")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x,y,z,w,v,t
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteArray6
      end interface
      call I_VTK_WriteArray6(this%object,pointid,x,y,z,w,v,t)
   end subroutine VTK_WriteArray6
   subroutine VTK_EndArray6(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndArray6(this) bind(C,name="C_VTK_EndArray6")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndArray6
      end interface
      call I_VTK_EndArray6(this%object)
   end subroutine VTK_EndArray6

   subroutine VTK_BeginArray9GP(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginArray9GP(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginArray9GP")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginArray9GP
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginArray9GP(this%object,n,C_filename)
   end subroutine VTK_BeginArray9GP
   subroutine VTK_WriteArray9GP(this,pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x1,x2,x3,x4,x5,x6,x7,x8,x9
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteArray9GP(this,pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9) bind(C,name="C_VTK_WriteArray9GP")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x1,x2,x3,x4,x5,x6,x7,x8,x9
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteArray9GP
      end interface
      call I_VTK_WriteArray9GP(this%object,pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9)
   end subroutine VTK_WriteArray9GP
   subroutine VTK_EndArray9GP(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndArray9GP(this) bind(C,name="C_VTK_EndArray9GP")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndArray9GP
      end interface
      call I_VTK_EndArray9GP(this%object)
   end subroutine VTK_EndArray9GP

   subroutine VTK_BeginArray9(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginArray9(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginArray9")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginArray9
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginArray9(this%object,n,C_filename)
   end subroutine VTK_BeginArray9
   subroutine VTK_WriteArray9(this,pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: x1,x2,x3,x4,x5,x6,x7,x8,x9
      integer, intent(in)                :: pointid
      interface
         subroutine I_VTK_WriteArray9(this,pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9) bind(C,name="C_VTK_WriteArray9")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: x1,x2,x3,x4,x5,x6,x7,x8,x9
            integer(C_int), value :: pointid
         end subroutine I_VTK_WriteArray9
      end interface
      call I_VTK_WriteArray9(this%object,pointid,x1,x2,x3,x4,x5,x6,x7,x8,x9)
   end subroutine VTK_WriteArray9
   subroutine VTK_EndArray9(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndArray9(this) bind(C,name="C_VTK_EndArray9")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndArray9
      end interface
      call I_VTK_EndArray9(this%object)
   end subroutine VTK_EndArray9

   subroutine VTK_BeginField(this,FORSTR_LEN,F_filename)
      type(CVTKInterface), intent(inout) :: this
      integer, intent(in)                :: FORSTR_LEN
      integer                            :: i,n
      character(FORSTR_LEN), intent(in)  :: F_filename
      character(1,kind=C_char)           :: C_filename(LEN_TRIM(F_filename)+1)
      interface
         subroutine I_VTK_BeginField(this,FORSTR_LEN,filename) bind(C,name="C_VTK_BeginField")
            import
            type(C_ptr),value :: this
            integer(C_int),value                  :: FORSTR_LEN
            character(1,kind=C_char), intent(in)  :: filename(FORSTR_LEN)
         end subroutine I_VTK_BeginField
      end interface
      n = LEN_TRIM(F_filename)
      do i = 1, n
         C_filename(i) = F_filename(i:i)
      end do
      C_filename(n+1) = C_NULL_char
      call I_VTK_BeginField(this%object,n,C_filename)
   end subroutine VTK_BeginField
   subroutine VTK_WriteField(this,scalar)
      type(CVTKInterface), intent(inout) :: this
      real(rp), intent(in)               :: scalar
      interface
         subroutine I_VTK_WriteField(this,scalar) bind(C,name="C_VTK_WriteField")
            import
            type(C_ptr),value     :: this
            real(C_double), value :: scalar
         end subroutine I_VTK_WriteField
      end interface
      call I_VTK_WriteField(this%object,scalar)
   end subroutine VTK_WriteField
   subroutine VTK_EndField(this)
      type(CVTKInterface), intent(inout) :: this
      interface
         subroutine I_VTK_EndField(this) bind(C,name="C_VTK_EndField")
            import
            type(C_ptr),value :: this
         end subroutine I_VTK_EndField
      end interface
      call I_VTK_EndField(this%object)
   end subroutine VTK_EndField

#endif
end module Mod_vtkPost
