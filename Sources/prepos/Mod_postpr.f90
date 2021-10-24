module Mod_postpr
  use typre
  use def_parame
  use Mod_iofile
  use Mod_Memor
  use Mod_Mesh 
  use Mod_int2str
  implicit none
  private
  public PostprFile, opengrafFileTrue
  
  type,abstract :: PostprFile 
      type(MemoryMan), pointer   :: Memor => NULL()
      integer(ip) :: lun_postp, lun_outpu_dom, lun_outpu_graf
      integer(ip) :: kfl_outfo, kfl_outco, kfl_flush,kfl_DuplicateInterfaceElements = 0
      integer(ip)    :: kfl_writeType= 0
      
      logical :: DidIOpenGrafFile = .false.
      character(150) :: fil_outpu_graf
      
      integer(ip) :: MPIrank, MPISize, MPIroot
      
      contains
      procedure :: SetMPI
      procedure :: SetOutputFormat
      procedure :: SetOutputCompression
      procedure :: SetFlush
      procedure :: SetDuplicateInterfaceElements
      procedure :: SetMulticomm
            
      procedure :: OpenPostFile
      procedure :: ClosePostFile
      procedure :: OpenGrafFile
      procedure :: CloseGrafFile
      procedure :: Initialize
      procedure :: Finalize

      procedure :: posfiel
      procedure :: possca
      procedure :: posvec
      procedure :: posvecbopo
      procedure :: posvec3
      procedure :: poisca
      procedure :: poisvec
      procedure :: posmat
      procedure :: MeshPostpr
      procedure :: bpossca
      procedure :: bposvec
      procedure :: bposvec3
      
      procedure :: pgprpsca
      procedure :: pgprpvec
      procedure :: pgpipsca
      procedure :: pgpr1p
      procedure :: pgpr2p
      procedure :: pgpr3p
      
      generic :: postpr  => posfiel,possca,posvec,posvecbopo,poisca,poisvec,posvec3,MeshPostpr
      generic :: bpostpr => bpossca,bposvec,bposvec3    !Boundary vectors
      generic :: postgp  => pgprpsca,pgprpvec,pgpr1p,pgpr2p,pgpr3p,pgpipsca
      
      procedure :: AverageToOneDimension
      procedure :: AverageToOneDimension1D
      
      generic :: postAvg1D => AverageToOneDimension, AverageToOneDimension1D
      
  end type
  
    interface
      
      subroutine AverageToOneDimension(a,AvgDime,array,wopos,istep,ttime,Mesh,Memor)
         use typre
         use Mod_Element
         use Mod_Memor
         use Mod_Mesh
         use MPI
         import PostprFile
         implicit none
         class(PostprFile) :: a
         integer(ip) :: AvgDime
         real(rp) :: array(:,:)
         integer(ip) :: istep
         character(*), intent(in) :: wopos
         real(rp) :: ttime
         class(FemMesh)   :: Mesh
         type(MemoryMan) :: Memor
      end subroutine
      
      subroutine AverageToOneDimension1D(a,AvgDime,array,wopos,istep,ttime,Mesh,Memor)
         use typre
         use Mod_Memor
         use Mod_Mesh
         import PostprFile
         implicit none
         class(PostprFile) :: a
         integer(ip) :: AvgDime
         real(rp), target :: array(:)
         character(*), intent(in) :: wopos
         integer(ip) :: istep
         real(rp) :: ttime
         class(FemMesh)   :: Mesh
         type(MemoryMan) :: Memor
      end subroutine
      
   end interface
  
contains

   subroutine SetMPI(a,MPIsize,MPIroot,MPIrank,Memor)
      implicit none
      class(PostprFile) :: a
      integer(ip) :: MPIsize,MPIrank,MPIroot
      type(MemoryMan), target   :: Memor
      
      a%Memor => Memor
      a%MPIsize = MPIsize
      a%MPIroot = MPIroot
      a%MPIrank = MPIrank
   end subroutine
   
   subroutine SetMulticomm(a,kfl_multicomm,MulticommColor)
      implicit none
      class(PostprFile) :: a
      integer(ip) :: kfl_multicomm, MulticommColor
      !Default is nothing to be done
   end subroutine

   subroutine Initialize(a,PostProcessFolder,namda)
      implicit none
      class(PostprFile) :: a
      character(150)    :: PostProcessFolder
      character(150)    :: namda
      
      call runend('Initialize not defined')
   end subroutine
   
   subroutine Finalize(a)
      implicit none
      class(PostprFile) :: a
      
      call runend('Finalize not defined')
   end subroutine
  
   subroutine SetOutputFormat(a,kfl_outfo)
      implicit none
      class(PostprFile) :: a
      integer(ip)        :: kfl_outfo

      a%kfl_outfo = kfl_outfo
   end subroutine 
  
   subroutine SetOutputCompression(a,kfl_outco,kfl_writeType)
      implicit none
      class(PostprFile) :: a
      integer(ip)        :: kfl_outco,kfl_writeType

      a%kfl_writeType= kfl_writeType
      a%kfl_outco = kfl_outco
   end subroutine 
   
   subroutine SetDuplicateInterfaceElements(a,kfl_DuplicateInterfaceElements)
      implicit none
      class(PostprFile) :: a
      integer(ip)        :: kfl_DuplicateInterfaceElements

      a%kfl_DuplicateInterfaceElements = kfl_DuplicateInterfaceElements
   end subroutine 
  
   subroutine OpenGrafFile(a,fil_outpu_graf)
      implicit none
      class(PostprFile) :: a
      character(150)     :: fil_outpu_graf
      
      a%fil_outpu_graf = fil_outpu_graf
      a%DidIOpenGrafFile = .false.
   end subroutine   
   
   subroutine OpenGrafFileTrue(a)
      implicit none
      class(PostprFile) :: a
      
      if (a%DidIOpenGrafFile .eqv. .false.) then
         a%DidIOpenGrafFile = .true.
         
         call iofile(zero,a%lun_outpu_graf,a%fil_outpu_graf,'OpenGrafFile')
      endif
   end subroutine   
  
   subroutine CloseGrafFile(a)
      implicit none
      class(PostprFile) :: a
      
      if (a%DidIOpenGrafFile .eqv. .true.) then
         call iofile(two,a%lun_outpu_graf,' ','CloseGrafFile','old')
      endif
   end subroutine
  
   subroutine SetFlush(a,kfl_flush)
      implicit none
      class(PostprFile) :: a
      integer(ip) :: kfl_flush
      
      a%kfl_flush = kfl_flush
   end subroutine
   
   subroutine OpenPostFile(a,PostProcessFolder,namda,meshCounter)
      implicit none
      class(PostprFile) :: a
      character(150)    :: PostProcessFolder
      character(150)    :: namda
      integer(ip),optional  :: meshCounter

      call runend('OpenPostFile not defined')
   end subroutine

   subroutine ClosePostFile(a)
      implicit none
      class(PostprFile) :: a

      call runend('ClosePostFile not defined')
   end subroutine    

   subroutine MeshPostpr(a,Mesh)
      implicit none
      class(PostprFile), intent(inout) :: a
      class(FemMesh), target, intent(inout)    :: Mesh

      call runend('MeshPostpr not defined')
   end subroutine

   subroutine inipos(a,MPIrank,MPIroot)
      use typre
      implicit none
      class(PostprFile) :: a
      integer(ip), intent(in) :: MPIrank, MPIroot

      call runend('inipos not defined')
   end subroutine

   subroutine endpos(a,istep,ctime,MPIrank,MPIroot,MPIsize)
      use typre
      implicit none
      class(PostprFile) :: a
      real(rp), intent(in)    :: ctime
      integer(ip), intent(in) :: istep, MPIrank, MPIsize, MPIroot

      call runend('endpos not defined')
   end subroutine

   subroutine closepos(a,MPIrank,MPIroot)
      implicit none
      class(PostprFile)   :: a
      integer(ip), intent(in) :: MPIrank, MPIroot

      call runend('closepos not defined')
   end subroutine
 
   subroutine posfiel(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)

      call runend('posfiel not defined')
   end subroutine posfiel

   subroutine possca(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)

      call runend('possca not defined')
   end subroutine possca

   subroutine posvec(a,bridge,wopos,istep,ttime,Mesh,auxString)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      character(6), optional    :: auxString

      call runend('posvec not defined')
   end subroutine posvec

   subroutine posvecbopo(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  intent(in)  :: kpoin,lpoty(:)

      call runend('posvecbopo not defined')
   end subroutine posvecbopo

   subroutine posvec3(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      use typre
      use Mod_int2str, only : int2str
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)

      call runend('posvec3 not defined')
   end subroutine posvec3

   subroutine posmat(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh) :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)

      call runend('posmat not defined')
   end subroutine posmat

   subroutine poisca(a,bridge,wopos,istep,ttime,Mesh)
      use typre
      implicit none
      class(FemMesh)    :: Mesh
      class(PostprFile) :: a
      character(*), intent(in)  :: wopos
      integer(ip),  intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: istep
      real(rp)    , intent(in)  :: ttime

      call runend('poisca not defined')
   end subroutine poisca

   subroutine poisvec(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      integer(ip),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)

      call runend('poisvec not defined')
   end subroutine poisvec

   subroutine bpossca(a,bridge,wopos,istep,ttime,Mesh)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime

      call runend('bpossca not defined')
   end subroutine

   subroutine bposvec(a,bridge,wopos,istep,ttime,Mesh)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime

      call runend('bposvec not defined')
   end subroutine

   subroutine bposvec3(a,bridge,wopos,istep,ttime,Mesh)
      use typre
      use Mod_int2str, only : int2str
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime

      call runend('bposvec3 not defined')
   end subroutine

   subroutine pgprpsca(a,bridge,wopos,itste,ttime,Mesh)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh  
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime

      call runend('pgprpsca not defined')
   end subroutine pgprpsca

   subroutine pgprpvec(a,bridge,wopos,itste,ttime,Mesh)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime

      call runend('pgprpvec not defined')
   end subroutine pgprpvec

   subroutine pgpipsca(a,bridge,wopos,itste,ttime,Mesh)
      use typre
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh  
      character(*), intent(in)  :: wopos
      integer(ip),     intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime

      call runend('pgpipsca not defined')
   end subroutine pgpipsca 
   subroutine pgpr1p(a,bridge,wopos,itste,ttime,Mesh)
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      type(r1p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime

      call runend('pgr1p not defined')
   end subroutine pgpr1p

   subroutine pgpr2p(a,bridge,wopos,itste,ttime,Mesh,auxString)
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh  
      character(*), intent(in)  :: wopos
      type(r2p),    intent(in)  :: bridge(:) 
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(6), optional    :: auxString

      call runend('pgr2p not defined')
   end subroutine pgpr2p

   subroutine pgpr3p(a,bridge,wopos,itste,ttime,Mesh,auxString)
      implicit none
      class(PostprFile) :: a
      class(FemMesh)    :: Mesh  
      character(*), intent(in)  :: wopos
      type(r3p),    intent(in)  :: bridge(:) 
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(6), optional    :: auxString

      call runend('pgr3p not defined')
   end subroutine pgpr3p


end module Mod_postpr

