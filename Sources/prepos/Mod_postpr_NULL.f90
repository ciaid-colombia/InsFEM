module Mod_postpr_NULL
   use typre
   use def_parame
   use Mod_postpr
   use Mod_iofile
   use Mod_Memor
   use Mod_Mesh
   use Mod_int2str
   implicit none
   private
   public PostprFile_NULL, PostprFile_NULL_const
 
   type, extends(PostprFile) ::  PostprFile_NULL
      integer(ip) :: ipagp = 0
      integer(ip) :: ipagp0 = 0

      contains
      procedure :: Initialize => Initialize
      procedure :: OpenPostFile   => OpenPostFile
      procedure :: ClosePostFile  => ClosePostFile
      procedure :: OpenGrafFile  => OpenGrafFileNULL
      procedure :: CloseGrafFile => CloseGrafFileNULL
      procedure :: Finalize => Finalize

      procedure :: MeshPostpr => geonull

      procedure :: posfiel    => posfiel
      procedure :: possca     => possca
      procedure :: posvec     => posvec
      procedure :: posvecbopo => posvecbopo
      procedure :: posvec3    => posvec3
      procedure :: poisca     => poisca
      procedure :: poisvec    => poisvec
      procedure :: posmat     => posmat

      procedure :: bpossca    => bpossca
      procedure :: bposvec    => bposvec
      procedure :: bposvec3   => bposvec3

      procedure :: pgprpsca   => pgprpsca
      procedure :: pgprpvec   => pgprpvec
      procedure :: pgpipsca   => pgpipsca
      procedure :: pgpr1p     => pgpr1p_NULL
      procedure :: pgpr2p     => pgpr2p_NULL
      procedure :: pgpr3p     => pgpr3p_NULL
   end type


   
   interface PostprFile_NULL_Const
       procedure constructor
   end interface PostprFile_NULL_Const

contains

   function constructor()
      class(PostprFile_NULL), pointer :: constructor
      allocate(constructor)
   end function constructor
   
   subroutine geonull(a,Mesh)
      class(PostprFile_NULL), intent(inout)  :: a
      class(FemMesh), target, intent(inout) :: Mesh
   end subroutine

   subroutine Initialize(a,PostProcessFolder,namda)
      implicit none
      class(PostprFile_NULL) :: a
      character(150)    :: PostProcessFolder
      character(150)    :: namda
     
   end subroutine
   
   subroutine OpenPostFile(a,PostProcessFolder,namda,meshCounter)
      implicit none
      class(PostprFile_NULL) :: a
      character(150)    :: PostProcessFolder
      character(150)    :: namda
      integer(ip),optional  :: meshCounter
     
   end subroutine

   subroutine ClosePostFile(a)
      implicit none
      class(PostprFile_NULL) :: a
 
   end subroutine
 
   subroutine Finalize(a)
      implicit none
      class(PostprFile_NULL) :: a
     
   end subroutine

   subroutine OpenGrafFileNULL(a,fil_outpu_graf)
      implicit none
      class(PostprFile_NULL) :: a
      character(150)     :: fil_outpu_graf

   end subroutine   
  
   subroutine CloseGrafFileNULL(a)
      implicit none
      class(PostprFile_NULL) :: a

   end subroutine

  
 
   subroutine inipos(a,MPIrank,MPIroot)
      use typre
      implicit none
      class(PostprFile_NULL) :: a
      integer(ip), intent(in) :: MPIrank, MPIroot

   end subroutine

   subroutine endpos(a,istep,ctime,MPIrank,MPIroot,MPIsize)
      use typre
      implicit none
      class(PostprFile_NULL) :: a
      real(rp), intent(in)    :: ctime
      integer(ip), intent(in) :: istep, MPIrank, MPIsize, MPIroot

   end subroutine

   subroutine closepos(a,MPIrank,MPIroot)
      implicit none
      class(PostprFile_NULL)   :: a
      integer(ip), intent(in) :: MPIrank, MPIroot

   end subroutine
 
   subroutine posfiel(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo,poini
      real(rp)                  :: dummr
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
 
     
   end subroutine posfiel
 
   subroutine possca(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo,poini
      real(rp)                  :: dummr
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
 
     
   end subroutine possca

   subroutine posvec(a,bridge,wopos,istep,ttime,Mesh,auxString)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      character(6), optional    :: auxString
      integer(ip)               :: npoin,ndime,ipoin,idime,poini
      real(rp)                  :: rz
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
 
   end subroutine posvec
 
   subroutine posvecbopo(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  intent(in)  :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo,poini
      real(rp)                  :: dummr,rz
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
 
   end subroutine posvecbopo
 
   subroutine posvec3(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      !Write a vector in postprocess file
      use Mod_int2str, only : int2str
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo,idime2,ndime2,poini
      real(rp)                  :: dummr,rz
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)

     
   end subroutine posvec3
 
   subroutine posmat(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh) :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ipoin,idime,kbopo,ibopo,auxdim,poini
      real(rp)                  :: dummr
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)

   end subroutine posmat

   subroutine poisca(a,bridge,wopos,istep,ttime,Mesh)
    
      !-----------------------------------------------------------------------
      !
      ! Write a scalar in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FemMesh)    :: Mesh
      class(PostprFile_NULL) :: a
      character(*), intent(in)  :: wopos
      integer(ip),  intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: istep
      real(rp)    , intent(in)  :: ttime
      integer(ip)               :: npoin,ipoin,poini
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)

   end subroutine poisca

   subroutine poisvec(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      !-----------------------------------------------------------------------
      !
      ! Write a vector in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      integer(ip),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo,poini
      real(rp)                  :: dummr,rz
      character(8)              :: state
      character(4)              :: NULL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
 
     
   end subroutine poisvec

   subroutine bpossca(a,bridge,wopos,istep,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a boundary scalar in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
 
   
   end subroutine
 
   subroutine bposvec(a,bridge,wopos,istep,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a boundary vector in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime

   end subroutine
  
   subroutine bposvec3(a,bridge,wopos,istep,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a boundary vector in postprocess file
      !
      !-----------------------------------------------------------------------
      use Mod_int2str, only : int2str
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
 
  end subroutine
  
   subroutine pgprpsca(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a scalar at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(8)              :: state
      character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,elemi,ielem,nelem,igaus,ndime
      character(13)             :: elemt
 
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
 
      real(rp), pointer :: posgp(:,:,:) => NULL()
      real(rp) :: rz
      integer(ip) :: npoin0
      integer(ip) :: pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: aux_lnods(27),minpoin,npoinLocal
 
      logical :: checkit
 
   
   end subroutine pgprpsca

   subroutine pgprpvec(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a vector at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(8)              :: state
      character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,elemi,ielem,nelem,igaus,ndime
      character(13)             :: elemt
 
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
 
      real(rp), pointer :: posgp(:,:,:) => NULL()
      real(rp) :: rz
      integer(ip) :: npoin0
      integer(ip) :: pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: aux_lnods(27),minpoin,npoinLocal
 
      logical :: checkit
 
   end subroutine pgprpvec

   subroutine pgpr1p_NULL(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a matrix at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      type(r1p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(8)              :: state
      character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,elemi,ielem,nelem,igaus,ndime
      character(13)             :: elemt

      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
 
      real(rp), pointer :: posgp(:,:,:) => NULL()

      real(rp) :: rz
      integer(ip) :: npoin0
      integer(ip) :: pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: aux_lnods(27),minpoin,npoinLocal
 
      logical :: checkit
 
   end subroutine pgpr1p_NULL

   subroutine pgpr2p_NULL(a,bridge,wopos,itste,ttime,Mesh,auxString)
      !-----------------------------------------------------------------------
      !
      ! Write a matrix at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      type(r2p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(8)              :: state
      character(6), optional    :: auxString
      character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,elemi,ielem,nelem,igaus,ndime
      character(13)             :: elemt
 
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
 
      real(rp), pointer :: posgp(:,:,:) => NULL()
 
       real(rp) :: rz
      integer(ip) :: npoin0
      integer(ip) :: pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: aux_lnods(27),minpoin,npoinLocal
 
      logical :: checkit
 
   
   end subroutine pgpr2p_NULL

   subroutine pgpr3p_NULL(a,bridge,wopos,itste,ttime,Mesh,auxString)
      !-----------------------------------------------------------------------
      !
      ! Write a matrix at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      type(r3p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(8)              :: state
      
      character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,elemi,ielem,nelem,igaus,ndime
      character(13)             :: elemt
 
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
      character(6), optional    :: auxString
 
      real(rp), pointer :: posgp(:,:,:) => NULL()
 
      real(rp) :: rz
      integer(ip) :: npoin0
      integer(ip) :: pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: aux_lnods(27),minpoin,npoinLocal
 
      logical :: checkit

    
   end subroutine pgpr3p_NULL

   subroutine pgpipsca(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a scalar at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_NULL) :: a
      class(FemMesh)    :: Mesh  
      character(*), intent(in)  :: wopos
      integer(ip),     intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(8)              :: state
      character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,elemi,ielem,nelem,igaus,ndime
      character(13)             :: elemt
      
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
      
      real(rp), pointer :: posgp(:,:,:) => NULL()
      real(rp) :: rz
      integer(ip) :: npoin0
      integer(ip) :: pnode
      integer(ip), pointer :: lnode(:) => NULL()
      integer(ip) :: aux_lnods(27),minpoin,npoinLocal
      
      logical :: checkit
    
   end subroutine pgpipsca

   subroutine CheckIfElementHasToBeWritten(Mesh,ielem,npoin0,npoinLocal,checkit)
      use typre
      implicit none
      class(FemMesh) :: Mesh
      integer(ip) :: ielem,npoin0,npoinLocal
      logical :: checkit
 
      integer(ip) :: minpoin,aux_lnods(27),pnode
      integer(ip), pointer :: lnode(:)
 
       
   end subroutine

   subroutine WriteGidGaussPoints(a,nelty,nnode,ngaus,ndime,posgp)
      use typre
      implicit none
      class(PostprFile_NULL) :: a
      integer(ip), pointer      :: nnode(:), ngaus(:)
      integer(ip)               :: nelty,ndime
      real(rp), pointer         :: posgp(:,:,:)
 
   end subroutine
 
end module Mod_postpr_NULL
