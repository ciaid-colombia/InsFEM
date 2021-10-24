module Mod_ReadWrite_F90
   use typre
   use Mod_ReadWrite
   implicit none
   private
   public Reader_F90, Writer_F90, Reader_F90_Const, Writer_F90_Const

   type, extends(Writer) :: Writer_F90

   contains
      procedure :: WriteAttributeReal
      procedure :: WriteAttributeInteger
      procedure :: WriteDataReal1
      procedure :: WriteDataReal2
      procedure :: WriteArray1
      procedure :: WriteArray2
      procedure :: WriteArray3
      procedure :: InitWriter
      procedure :: EndWriter

   end type

   type, extends(Reader) :: Reader_F90

   contains
      procedure :: ReadAttributeReal
      procedure :: ReadAttributeInteger
      procedure :: ReadDataReal1
      procedure :: ReadDataReal2
      procedure :: ReadArray1
      procedure :: ReadArray2
      procedure :: ReadArray3
      procedure :: InitReader
      procedure :: EndReader

   end type

   interface Writer_F90_Const
      procedure constructorW
   end interface Writer_F90_Const

   interface Reader_F90_Const
      procedure constructorR
   end interface Reader_F90_Const

contains

   function constructorW()
        class(Writer_F90), pointer :: constructorW

        allocate(constructorW)

    end function constructorW

   function constructorR()
        class(Reader_F90), pointer :: constructorR

        allocate(constructorR)

    end function constructorR

   subroutine InitWriter(a)
      implicit none
      class(Writer_F90) :: a
   end subroutine
   
   subroutine EndWriter(a)
      implicit none
      class(Writer_F90) :: a
   end subroutine
   
   subroutine WriteAttributeInteger(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,s
      integer(ip)      :: Array(:)

      !write(fil_unit) Array
      call runend('Attributes not supported')
   end subroutine

   subroutine WriteAttributeReal(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i
      real(rp)         :: Array(:)

      !write(fil_unit) Array
      call runend('Attributes not supported')
   end subroutine

   subroutine WriteDataReal2(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,j,s(2)
      real(rp)         :: Array(:,:)

      s(1) = size(Array,1)
      s(2) = size(Array,2)
      write(fil_unit) ((Array(i,j),i=1,s(1)),j=1,s(2))
   end subroutine

   subroutine WriteDataReal1(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,s
      real(rp)         :: Array(:)

      s = size(Array,1)
      write(fil_unit) (Array(i),i=1,s)
   end subroutine

   subroutine WriteArray1(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i
      real(rp)         :: Array(:)
      integer(ip), optional :: ncol

      write(fil_unit) (Array(i),i=1,a%sizeLocal)
   end subroutine

   subroutine WriteArray2(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,s,i,j
      real(rp)         :: Array(:,:)
      integer(ip), optional :: ncol

      s = size(Array,2)
      write(fil_unit) ((Array(i,j),i=1,a%sizeLocal),j=1,s)
   end subroutine

   subroutine WriteArray3(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      use typre
      implicit none
      class(Writer_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,s(2),i,j,k
      real(rp)         :: Array(:,:,:)
      integer(ip), optional :: ncol

      s(1) = size(Array,2)
      s(2) = size(Array,3)
      write(fil_unit) (((Array(i,j,k),i=1,a%sizeLocal),j=1,s(1)),k=1,s(2))
   end subroutine

   subroutine InitReader(a)
      implicit none
      class(Reader_F90) :: a
   end subroutine
   
   subroutine EndReader(a)
      implicit none
      class(Reader_F90) :: a
   end subroutine

   subroutine ReadAttributeReal(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,s
      real(rp)         :: Array(:)

      !read(fil_unit) Array(i)
      call runend('Attributes not supported')
   end subroutine

   subroutine ReadAttributeInteger(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i
      integer(ip)         :: Array(:)

      !read(fil_unit) Array
      call runend('Attributes not supported')
   end subroutine

   subroutine ReadDataReal2(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,j,s(2)
      real(rp)         :: Array(:,:)

      s(1) = size(Array,1)
      s(2) = size(Array,2)
      read(fil_unit) ((Array(i,j),i=1,s(1)),j=1,s(2))
   end subroutine

   subroutine ReadDataReal1(a,fil_name,fil_unit,Array,nameV,nameG)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,s
      real(rp)         :: Array(:)

      s = size(Array,1)
      read(fil_unit) (Array(i),i=1,s)
   end subroutine

   subroutine ReadArray1(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i
      real(rp)         :: Array(:)
      integer(ip), optional :: ncol

      read(fil_unit) (Array(i),i=1,a%sizeLocal)
   end subroutine

   subroutine ReadArray2(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,i,j,s
      real(rp)         :: Array(:,:)
      integer(ip), optional :: ncol

      s = size(Array,2)
      read(fil_unit) ((Array(i,j),i=1,a%sizeLocal),j=1,s)
   end subroutine

   subroutine ReadArray3(a,fil_name,fil_unit,Array,nameV,nameG,ncol)
      use typre
      implicit none
      class(Reader_F90) :: a
      character(150)   :: fil_name,nameV,nameG
      integer(ip)      :: fil_unit,s(2),i,j,k
      real(rp)         :: Array(:,:,:)
      integer(ip), optional :: ncol

      s(1) = size(Array,2)
      s(2) = size(Array,3)
      read(fil_unit) (((Array(i,j,k),i=1,a%sizeLocal),j=1,s(1)),k=1,s(2))
   end subroutine

end module Mod_ReadWrite_F90
