module Mod_postpr_VTK 
   use typre
   use def_parame
   use Mod_iofile
   use Mod_Memor
   use Mod_Mesh
   use Mod_Element
   use Mod_int2str
   use Mod_postpr
   use Mod_vtkPost

   implicit none
   private
   public PostprFile_VTK, PostprFile_VTK_Const

   type, extends(PostprFile) ::  PostprFile_VTK

      type(CVTKInterface) :: VTKInterface

      integer(ip)    :: ipagp = 0
      integer(ip)    :: kfl_written = 0
      integer(ip)    :: isteppost = 0
      character(150) :: fil_postp
      character(150) :: namda
      
      class(FemMesh), pointer :: Mesh => NULL()
      
      integer(ip) :: istep = -20
      real(rp) :: ctime = -1e5
      
   contains

#ifdef VTK

      procedure :: Initialize => Initialize
      procedure :: MeshPostpr => MeshPostprVTK
      procedure :: OpenPostFile   => OpenPostFile
      procedure :: ClosePostFile  => ClosePostFile
      procedure :: Finalize => Finalize
      procedure :: SetMulticomm

      procedure :: posfiel    => posfiel
      procedure :: possca     => possca
      procedure :: posvec     => posvec
      procedure :: posvecbopo => posvecbopo
      procedure :: poisca     => poisca
      procedure :: poisvec    => poisvec

      procedure :: bpossca    => bpossca
      procedure :: bposvec    => bposvec
      procedure :: bposvec3   => bposvec3

      procedure :: pgprpsca   => pgprpsca
      procedure :: pgprpvec   => pgprpvec
      procedure :: pgpipsca   => pgpipsca
      procedure :: pgpr1p     => pgpr1p
      procedure :: pgpr2p     => pgpr2p
      procedure :: pgpr3p     => pgpr3p
      
#endif

   end type


   interface
      subroutine geovtk(a,Mesh)
            import FemMesh
            import PostprFile_VTK
            class(PostprFile_VTK), intent(inout) :: a
            class(FemMesh), intent(inout) :: Mesh
        end subroutine
    end interface

    interface PostprFile_VTK_Const
        procedure constructor
    end interface PostprFile_VTK_Const

contains

   function constructor()
        class(PostprFile_VTK), pointer :: constructor

        allocate(constructor)

    end function constructor
 
#ifdef VTK

   subroutine MeshPostprVTK(a,Mesh)
      implicit none
      class(POstprFile_VTK), intent(inout) :: a
      class(FemMesh), target, intent(inout) :: Mesh
      
      call PointMesh(a,Mesh)
      call geovtk(a,a%Mesh)
      
   end subroutine

   subroutine Initialize(a,PostProcessFolder,namda)
      implicit none
      class(PostprFile_VTK) :: a
      character(150)        :: PostProcessFolder
      character(150)        :: namda
      integer(ip)           :: strlen
      !logical, save :: notInitialized = .true.
      

      a%fil_postp = trim(PostProcessFolder)//'/'//adjustl(trim(namda))
      a%namda = namda
         
      call VTK_Interface(a%VTKInterface)
      !if (notInitialized) then
      !TODO, this needs to be fixed, apparently the pointer is not being recognized
      !correctly when cases are used and it tries to initialize the same MPI more than
      !once
      call VTK_SetMPI(a%VTKInterface,a%MPIsize,a%MPIrank)
          !notInitialized = .false.
      !endif 

      if (a%MPIrank == a%MPIroot) then 
         strlen = len(trim(adjustl(a%namda)))
         call VTK_SetTimeWriter(a%VTKInterface,strlen,trim(adjustl(a%namda)))
      end if
         
   end subroutine
   
   subroutine SetMulticomm(a,kfl_multicomm,MulticommColor)
      implicit none
      class(PostprFile_VTK) :: a
      integer(ip) :: kfl_multicomm, MulticommColor
      
      call VTK_SetMulticomm(a%VTKInterface,kfl_multicomm,MulticommColor)
   end subroutine

   subroutine OpenPostFile(a,PostProcessFolder,namda,meshCounter)
      implicit none
      class(PostprFile_VTK) :: a
      character(150)        :: PostProcessFolder
      character(150)        :: namda
      integer(ip)           :: strlen
      integer(ip),optional  :: meshCounter

      select case (a%kfl_outco)
      case (Ascii)
         call VTK_BeginWriter(a%VTKInterface,Ascii,a%kfl_writeType)
      case (Binary)
         call VTK_BeginWriter(a%VTKInterface,Binary,a%kfl_writeType)
      case (Appended)
         call VTK_BeginWriter(a%VTKInterface,Appended,a%kfl_writeType)
      end select
      
   end subroutine
   
   subroutine EndPrePos(a,istep,ctime)
      implicit none
      class(PostprFile_VTK) :: a
      integer(ip) :: istep
      real(rp)    :: ctime
      
      if (a%istep /= istep .or. a%ctime /= ctime) then
         a%istep = istep
         a%ctime = ctime
         call endpos(a)
      endif

   end subroutine

   subroutine ClosePostFile(a)
      implicit none
      class(PostprFile_VTK) :: a

      call endpos(a)
      call VTK_Reset(a%VTKInterface)

   end subroutine

   subroutine Finalize(a)
      implicit none
      class(PostprFile_VTK) :: a
      integer(ip)             :: strlen

      if (a%MPIrank == a%MPIroot) then 
         strlen = len(trim(adjustl(a%namda)))
         call VTK_EndTimeWriter(a%VTKInterface,strlen,trim(adjustl(a%namda)))
      end if
      call VTK_Interface_des(a%VTKInterface)

   end subroutine

   subroutine endpos(a)
      implicit none
      class(PostprFile_VTK) :: a
      integer(ip)             :: strlen,npoin
      strlen = len(trim(adjustl(a%fil_postp))//"_"//trim(int2str(a%isteppost)))
      call VTK_SetWriteFile(a%VTKInterface,strlen,trim(adjustl(a%fil_postp))//"_"//trim(int2str(a%isteppost)))
      if (a%kfl_written /= 0) then
         if (a%MPIrank == a%MPIroot) then 
            strlen = len(trim(adjustl(a%namda)))
            call VTK_AddTime(a%VTKInterface,strlen,trim(adjustl(a%namda)),a%ctime,a%isteppost) ! Set name for VTK time step
         end if
         if(a%kfl_writeType== 1) call VTK_EndBlocks(a%VTKInterface)
         call VTK_Flush(a%VTKInterface)
         a%isteppost = a%isteppost +1
      end if
      a%kfl_written = 0

   end subroutine
      
   subroutine PointMesh(a,Mesh)
      implicit none
      class(POstprFile_VTK), intent(inout) :: a
      class(FemMesh), intent(inout), target :: Mesh
      
      a%Mesh => Mesh
      
   end subroutine

   subroutine posfiel(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo
      real(rp)                  :: dummr
 
      !If not open, open the postfile
      call EndPrePos(a,istep,ttime)
 
      call VTK_BeginField(a%VTKInterface,len_trim(wopos),trim(wopos))
      call VTK_WriteField(a%VTKInterface,bridge)
      call VTK_EndField(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1
 
   end subroutine posfiel
 
   subroutine possca(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo
      real(rp)                  :: dummr
 
      !If not open, open the postfile
      call EndPrePos(a,istep,ttime)
 
      !Dimensions
      if(present(kpoin)) then
         npoin=kpoin
      else
         npoin=size(bridge)
      end if
 
      !Check if there is arguments
      if(present(lpoty)) then
         kbopo=1               ! bridge only defined on boundary nodes
         dummr=0.0_rp
      else
         kbopo=0               ! bridge is defined on all nodes
      end if
 
      call VTK_BeginScalar(a%VTKInterface,len_trim(wopos),trim(wopos))
      if(kbopo==0) then
         do ipoin=1,npoin
            call VTK_WriteScalar(a%VTKInterface,ipoin-1,bridge(ipoin))
         end do
      else
         do ipoin = 1,npoin
            ibopo=lpoty(ipoin)
            if(ibopo/=0) then
               call VTK_WriteScalar(a%VTKInterface,ipoin-1,bridge(ibopo))
            else
               call VTK_WriteScalar(a%VTKInterface,ipoin-1,dummr)
            end if
         end do
      end if
      call VTK_EndScalar(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1
 
   end subroutine possca

   subroutine posvec(a,bridge,wopos,istep,ttime,Mesh,auxString)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_VTK)     :: a
      class(FemMesh)            :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      character(6), optional    :: auxString
      integer(ip)               :: npoin,ndime,ipoin,idime,sz
      real(rp)                  :: rz,rzt1,rzt2,rzt3

      !If not open, open the postfile
      call EndPrePos(a,istep,ttime)
 
      call Mesh%GetNdime(ndime)
      npoin=size(bridge,2)

      if (present(auxString)) then

          if (auxString == 'voigtT') then
              sz = (ndime*(ndime+1))/2
              rz   = 0.0_rp
              rzt1 = 0.0_rp
              rzt2 = 0.0_rp
              rzt3 = 0.0_rp

              call VTK_BeginArray6(a%VTKInterface,len_trim(wopos),trim(wopos))

              do ipoin=1,npoin
                  if (ndime == 3) then 
                      rz   = bridge(3,ipoin)
                      rzt1 = bridge(4,ipoin)
                      rzt2 = bridge(5,ipoin)
                      rzt3 = bridge(6,ipoin)
                  else
                      rzt1 = bridge(3,ipoin)
                  end if
                  call VTK_WriteArray6(a%VTKInterface,ipoin-1,bridge(1,ipoin),bridge(2,ipoin),rz,rzt1,rzt2,rzt3)
              end do

              call VTK_EndArray6(a%VTKInterface)

          end if
      else
          call VTK_BeginVector(a%VTKInterface,len_trim(wopos),trim(wopos))

          rz = 0.0_rp
          do ipoin=1,npoin
              if (ndime == 3) rz = bridge(3,ipoin)
              call VTK_WriteVector(a%VTKInterface,ipoin-1,bridge(1,ipoin),bridge(2,ipoin),rz)
          end do

          call VTK_EndVector(a%VTKInterface)
      end if
 
      a%kfl_written = a%kfl_written + 1
 
   end subroutine posvec
 
   subroutine posvecbopo(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      ! Write a vector in postprocess file
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  intent(in)  :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo
      real(rp)                  :: dummr,rz

      !If not open, open the postfile
      call EndPrePos(a,istep,ttime)
 
      !Dimensions
      npoin=kpoin
      ndime=size(bridge,1)

      kbopo=1               ! bridge only defined on boundary nodes
      dummr=0.0_rp

      call VTK_BeginVector(a%VTKInterface,len_trim(wopos),trim(wopos))

      rz = 0.0_rp
      do ipoin = 1,npoin
          ibopo=lpoty(ipoin)
          if(ibopo/=0) then
              if (ndime == 3) rz = bridge(3,ibopo)
              call VTK_WriteVector(a%VTKInterface,ipoin-1,bridge(1,ibopo),bridge(2,ibopo),rz)
          else
              call VTK_WriteVector(a%VTKInterface,ipoin-1,dummr,dummr,dummr)
          end if
      end do

      call VTK_EndVector(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1
 
   end subroutine posvecbopo

   subroutine poisca(a,bridge,wopos,istep,ttime,Mesh)
    
      !-----------------------------------------------------------------------
      !
      ! Write a scalar in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FemMesh)    :: Mesh
      class(PostprFile_VTK) :: a
      character(*), intent(in)  :: wopos
      integer(ip),  intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: istep
      real(rp)    , intent(in)  :: ttime
      integer(ip)               :: npoin,ipoin

      !If not open, open the postfile
      call EndPrePos(a,istep,ttime)
      
      npoin=size(bridge)
      call VTK_BeginScalar(a%VTKInterface,len_trim(wopos),trim(wopos))
      do ipoin=1,npoin
         call VTK_WriteScalar(a%VTKInterface,ipoin-1,real(bridge(ipoin),8))
      end do
      call VTK_EndScalar(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1

   end subroutine poisca

   subroutine poisvec(a,bridge,wopos,istep,ttime,Mesh,kpoin,lpoty)
      !-----------------------------------------------------------------------
      !
      ! Write a vector in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      integer(ip),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
      integer(ip),  optional    :: kpoin,lpoty(:)
      integer(ip)               :: npoin,ndime,ipoin,idime,kbopo,ibopo
      real(rp)                  :: dummr,rz
 
      !If not open, open the postfile
      call EndPrePos(a,istep,ttime)
 
      !Dimensions
      if(present(kpoin)) then
         npoin=kpoin
      else
         npoin=size(bridge,2)
      end if
      ndime=size(bridge,1)
 
      !Check if there is arguments
      if(present(lpoty)) then
         kbopo=1               ! bridge only defined on boundary nodes
         dummr=0.0_rp
      else
         kbopo=0               ! bridge is defined on all nodes
      end if
      rz = 0.0_rp
      call VTK_BeginVector(a%VTKInterface,len_trim(wopos),trim(wopos))
      if(kbopo==0) then
         do ipoin=1,npoin
            if (ndime == 3) rz = real(bridge(3,ipoin),8)
            call VTK_WriteVector(a%VTKInterface,ipoin-1,real(bridge(1,ipoin),8),real(bridge(2,ipoin),8),rz)
         end do
      else
         do ipoin = 1,npoin
            ibopo=lpoty(ipoin)
            if(ibopo/=0) then
               if (ndime == 3) rz = real(bridge(3,ipoin),8)
               call VTK_WriteVector(a%VTKInterface,ipoin-1,real(bridge(1,ibopo),8),real(bridge(2,ibopo),8),rz)
            else
               call VTK_WriteVector(a%VTKInterface,ipoin-1,dummr,dummr,dummr)
            end if
         end do
      end if
      call VTK_EndVector(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1

   end subroutine poisvec

   subroutine bpossca(a,bridge,wopos,istep,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a boundary scalar in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
 
      call a%possca(bridge,adjustl(trim(wopos)),istep,ttime,Mesh,Mesh%npoin,Mesh%lpoty)
 
   end subroutine
 
   subroutine bposvec(a,bridge,wopos,istep,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a boundary vector in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime

      call a%posvecbopo(bridge,adjustl(trim(wopos)),istep,ttime,Mesh,Mesh%npoin,Mesh%lpoty)

   end subroutine
  
   subroutine bposvec3(a,bridge,wopos,istep,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a boundary vector in postprocess file
      !
      !-----------------------------------------------------------------------
      use Mod_int2str, only : int2str
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:,:)
      integer(ip),  intent(in)  :: istep
      real(rp),     intent(in)  :: ttime
 
      call a%posvec3(bridge,adjustl(trim(wopos)),istep,ttime,Mesh,Mesh%npoin,Mesh%lpoty)
 
  end subroutine
  
   subroutine pgprpsca(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a scalar at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      integer(ip)               :: ielty,jelty,ielem,nelem
 
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
 
      !If not open, open the postfile
      call EndPrePos(a,itste,ttime)
      call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
      call Mesh%GetNelem(nelem)

      call VTK_BeginScalarGP(a%VTKInterface,len_trim(wopos),trim(wopos))
 
      jelty = 1
      do ielty=1,nelty
         do ielem = 1,nelem
            call Mesh%GetIelty(ielem,jelty)
            if(jelty==ielty) then
               call VTK_WriteScalarGP(a%VTKInterface,ielem-1,bridge(ielem))
            end if
         enddo
      enddo
      call VTK_EndScalarGP(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1
   end subroutine pgprpsca

   subroutine pgprpvec(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a vector at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      real(rp),     intent(in)  :: bridge(:,:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,ielem,nelem,igaus,ndime
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
 
      real(rp) :: rz
      integer(ip) :: pnode
 
      !If not open, open the postfile
      call EndPrePos(a,itste,ttime)
 
      call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
      call Mesh%GetNdime(ndime)
      call Mesh%GetNelem(nelem)
 
      !Gauss point headers
      call VTK_BeginVectorGP(a%VTKInterface,len_trim(wopos),trim(wopos))

      rz = 0.0_rp
 
      jelty = 1
      do ielty=1,nelty
         do ielem = 1,nelem
            call Mesh%GetIelty(ielem,jelty)
            if(jelty==ielty) then
               if (ndime == 3) rz = bridge(3,ielem)
               call VTK_WriteVectorGP(a%VTKInterface,ielem-1,bridge(1,ielem),bridge(2,ielem),rz)
            end if
         enddo
      enddo
      call VTK_EndVectorGP(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1
   end subroutine pgprpvec

   subroutine pgpr1p(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a matrix at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK)     :: a
      class(FemMesh)            :: Mesh
      character(*), intent(in)  :: wopos
      type(r1p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      integer(ip),  save        :: ipass=0
      integer(ip)               :: iesta,iesto,ielty,icomp,ncomp,jelty,ielem,nelem,igaus,ndime
      character(13)             :: elemt
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty,pgaus,pnode
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable     :: elvar(:,:), gaussArray(:,:)
      real(rp)                  :: gpvar(1)
 
      !If not open, open the postfile
      call EndPrePos(a,itste,ttime)
 
      call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
      call Mesh%GetNdime(ndime)
      call Mesh%GetNelem(nelem)
 
      call VTK_BeginScalarGP(a%VTKInterface,len_trim(wopos),trim(wopos))
      !Things to be done prior to refinement
      call Mesh%ElementAlloc(e,a%Memor,'DefaultRule','postpr_VTK')
      call a%Memor%alloc(1,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
      call a%Memor%alloc(1,e%mnode,elvar,'elvar','postpr_VTK') 
 
      jelty = 1
      do ielty=1,nelty
         do ielem = 1,nelem
            call Mesh%GetIelty(ielem,jelty)
            if(jelty==ielty) then
                call Mesh%ElementLoad(ielem,e)
                call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
                gaussArray(1,1:pgaus) = bridge(ielem)%a(1:pgaus)
                call e%GaussToNodes(gaussArray,elvar)
                call e%interpc(1,elvar,gpvar(1))
                call VTK_WriteScalarGP(a%VTKInterface,ielem-1,gpvar(1))
            end if
         enddo
      enddo

      call a%Memor%dealloc(1,e%pgaus,gaussArray,'gaussArray','postpr_VTK') 
      call a%Memor%dealloc(1,e%mnode,elvar,'elvar','postpr_VTK') 
      call Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','postpr_VTK')

      call VTK_EndScalarGP(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1

   end subroutine pgpr1p

   subroutine pgpr2p(a,bridge,wopos,itste,ttime,Mesh,auxString)
      !-----------------------------------------------------------------------
      !
      ! Write a matrix at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      type(r2p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      character(6), optional    :: auxString
      integer(ip)               :: ielty,icomp,ncomp,jelty,ielem,nelem,igaus,ndime,iboun,nboun
      character(13)             :: elemt
      real(rp)                  :: rz,rzt1,rzt2,rzt3
      integer(ip)               :: nelty,pgaus,pnode,sz
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable     :: elvar(:,:),gpvar(:),gaussArray(:,:)
 
      !If not open, open the postfile
      call EndPrePos(a,itste,ttime)
 
      call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
      call Mesh%GetNdime(ndime)
      call Mesh%GetNelem(nelem)
 
      call Mesh%ElementAlloc(e,a%Memor,'DefaultRule','postpr_VTK')
 
      rz = 0.0_rp
      jelty = 1
      if (present(auxString)) then
         if (auxString == 'SCALAR') then
            call a%Memor%alloc(1,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
            call a%Memor%alloc(1,e%mnode,elvar,'elvar','postpr_VTK') 
            call a%Memor%alloc(1,gpvar,'gpvar','postpr_VTK') 
            call VTK_BeginScalarGP(a%VTKInterface,len_trim(wopos),trim(wopos))
            do ielty=1,nelty
               do ielem = 1,nelem
                  call Mesh%GetIelty(ielem,jelty)
                  if(jelty==ielty) then
                     call Mesh%ElementLoad(ielem,e)
                     call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
                     gaussArray(:,1:pgaus) = bridge(ielem)%a(:,1:pgaus)
                     call e%GaussToNodes(gaussArray,elvar)
                     call e%interpc(1,elvar,gpvar(1))
                     call VTK_WriteScalarGP(a%VTKInterface,ielem-1,real(gpvar(1),8))
                  end if
               enddo
            enddo
            call a%Memor%dealloc(1,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
            call a%Memor%dealloc(1,e%mnode,elvar,'elvar','postpr_VTK') 
            call a%Memor%dealloc(1,gpvar,'gpvar','postpr_VTK') 
            call VTK_EndScalarGP(a%VTKInterface)
         elseif (auxString == 'voigtT') then
            sz = (ndime*(ndime+1))/2
            rzt1 = 0.0_rp
            rzt2 = 0.0_rp
            rzt3 = 0.0_rp
            call a%Memor%alloc(sz,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
            call a%Memor%alloc(sz,e%mnode,elvar,'elvar','postpr_VTK') 
            call a%Memor%alloc(sz,gpvar,'gpvar','postpr_VTK') 
            call VTK_BeginArray6GP(a%VTKInterface,len_trim(wopos),trim(wopos))
            do ielty=1,nelty
               do ielem = 1,nelem
                  call Mesh%GetIelty(ielem,jelty)
                  if(jelty==ielty) then
                     call Mesh%ElementLoad(ielem,e)
                     call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
                     gaussArray(:,1:pgaus) = bridge(ielem)%a(:,1:pgaus)
                     call e%GaussToNodes(gaussArray,elvar)
                     call e%interpc(sz,elvar,gpvar)
                     if (ndime == 3) then 
                        rz = gpvar(3)
                        rzt1 = gpvar(4)
                        rzt2 = gpvar(5)
                        rzt3 = gpvar(6)
                     else
                        rzt1 = gpvar(3)
                     end if
                     call VTK_WriteArray6GP(a%VTKInterface,ielem-1,gpvar(1),gpvar(2),rz,rzt1,rzt2,rzt3)
                  end if
               enddo
            enddo
            call a%Memor%dealloc(sz,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
            call a%Memor%dealloc(sz,e%mnode,elvar,'elvar','postpr_VTK') 
            call a%Memor%dealloc(sz,gpvar,'gpvar','postpr_VTK') 
            call VTK_EndArray6GP(a%VTKInterface)
         elseif (auxString == 'BOUNDA') then
            call a%Mesh%GetNboun(nboun)   
            call a%Memor%alloc(ndime,Mesh%mgaub,gaussArray,'gaussArray','postpr_VTK') 
            call a%Memor%alloc(ndime,e%mnodb,elvar,'elvar','postpr_VTK') 
            call a%Memor%alloc(ndime,gpvar,'gpvar','postpr_VTK') 
            call VTK_BeginVectorGP(a%VTKInterface,len_trim(wopos),trim(wopos))
            do ielty=1,nelty
               do ielem = 1,nelem
                  call Mesh%GetIelty(ielem,jelty)
                  if(jelty==ielty) then
                     call Mesh%ElementLoad(ielem,e)
                     call VTK_WriteVectorGP(a%VTKInterface,ielem-1,0.0_rp,0.0_rp,0.0_rp)
                  end if
               enddo
               do iboun = 1,nboun
                  call Mesh%BoundaryLoad(iboun,e)
                  call Mesh%GetBoundaryIelem(iboun,ielem)
                  call Mesh%GetIelty(ielem,jelty)
                  if (jelty==ielty) then
                     call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
                     gaussArray(:,1:pgaus) = bridge(ielem)%a(:,1:pgaus)
                     call e%GaussToNodesB(gaussArray,elvar)
                     call e%interpc(1,elvar,gpvar(1))
                     if (ndime == 3) rz = gpvar(3)
                     call VTK_WriteVectorGP(a%VTKInterface,ielem-1,gpvar(1),gpvar(2),rz)
                  end if
               enddo
            enddo
            call a%Memor%dealloc(ndime,Mesh%mgaub,gaussArray,'gaussArray','postpr_VTK') 
            call a%Memor%dealloc(ndime,e%mnodb,elvar,'elvar','postpr_VTK') 
            call a%Memor%dealloc(ndime,gpvar,'gpvar','postpr_VTK') 
            call VTK_EndVectorGP(a%VTKInterface)
         end if
      else
         call a%Memor%alloc(ndime,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
         call a%Memor%alloc(ndime,e%mnode,elvar,'elvar','postpr_VTK') 
         call a%Memor%alloc(ndime,gpvar,'gpvar','postpr_VTK') 
         call VTK_BeginVectorGP(a%VTKInterface,len_trim(wopos),trim(wopos))
            do ielty=1,nelty
               do ielem = 1,nelem
                  call Mesh%GetIelty(ielem,jelty)
                  if(jelty==ielty) then
                     call Mesh%ElementLoad(ielem,e)
                     call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
                     gaussArray(:,1:pgaus) = bridge(ielem)%a(:,1:pgaus)
                     call e%GaussToNodes(gaussArray,elvar)
                     call e%interpc(ndime,elvar,gpvar)
                     if (ndime == 3) rz = gpvar(3)
                     call VTK_WriteVectorGP(a%VTKInterface,ielem-1,gpvar(1),gpvar(2),rz)
                  end if
               enddo
            enddo
         call a%Memor%dealloc(ndime,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
         call a%Memor%dealloc(ndime,e%mnode,elvar,'elvar','postpr_VTK') 
         call a%Memor%dealloc(ndime,gpvar,'gpvar','postpr_VTK') 
         call VTK_EndVectorGP(a%VTKInterface)
      end if
      call Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','postpr_VTK')
      a%kfl_written = a%kfl_written + 1

   end subroutine pgpr2p

   subroutine pgpr3p(a,bridge,wopos,itste,ttime,Mesh,auxString)
      !-----------------------------------------------------------------------
      !
      ! Write a matrix at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh
      character(*), intent(in)  :: wopos
      type(r3p),    intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      integer(ip)               :: ielty,jelty,ielem,nelem,igaus,ndime
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty,pgaus,pnode
      character(6), optional    :: auxString
      real(rp)                  :: rz
      class(FiniteElement), pointer :: e => NULL()
      real(rp), allocatable     :: elvar(:,:),gpvar(:), gaussArray(:,:)
 
      !If not open, open the postfile
      call EndPrePos(a,itste,ttime)
 
      call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
      call Mesh%GetNdime(ndime)
      call Mesh%GetNelem(nelem)
 
      !Gauss point headers

      call Mesh%ElementAlloc(e,a%Memor,'DefaultRule','postpr_VTK')
 
      rz = 0.0_rp
      jelty = 1
      if (present(auxString)) then
      else    
         call a%Memor%alloc(ndime,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
         call a%Memor%alloc(ndime,e%mnode,elvar,'elvar','postpr_VTK') 
         call a%Memor%alloc(ndime,gpvar,'gpvar','postpr_VTK') 
         call VTK_BeginVectorGP(a%VTKInterface,len_trim(wopos),trim(wopos))
         do ielty=1,nelty
            do ielem = 1,nelem
               call Mesh%GetIelty(ielem,jelty)
               if(jelty==ielty) then
                  call Mesh%ElementLoad(ielem,e)
                  call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
                  gaussArray(:,1:pgaus) = bridge(ielem)%a(:,1,1:pgaus)
                  call e%GaussToNodes(gaussArray,elvar)
                  call e%interpc(ndime,elvar,gpvar)
                  if (ndime == 3) rz = gpvar(3)
                  call VTK_WriteVectorGP(a%VTKInterface,ielem,gpvar(1),gpvar(2),rz)
               end if
            enddo
         enddo
         call a%Memor%dealloc(ndime,Mesh%mgaus,gaussArray,'gaussArray','postpr_VTK') 
         call a%Memor%dealloc(ndime,e%mnode,elvar,'elvar','postpr_VTK') 
         call a%Memor%dealloc(ndime,gpvar,'gpvar','postpr_VTK') 
         call Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','postpr_VTK')
         call VTK_EndVectorGP(a%VTKInterface)
      end if   
      a%kfl_written = a%kfl_written + 1
   end subroutine pgpr3p

   subroutine pgpipsca(a,bridge,wopos,itste,ttime,Mesh)
      !-----------------------------------------------------------------------
      !
      ! Write a scalar at Gauss points in postprocess file
      !
      !-----------------------------------------------------------------------
      implicit none
      class(PostprFile_VTK) :: a
      class(FemMesh)    :: Mesh  
      character(*), intent(in)  :: wopos
      integer(ip),     intent(in)  :: bridge(:)
      integer(ip) , intent(in)  :: itste
      real(rp)    , intent(in)  :: ttime
      integer(ip)               :: ielty,jelty,ielem,nelem
      
      integer(ip), pointer      :: nnode(:) => NULL(), ngaus(:) => NULL(), ltopo(:) => NULL()
      integer(ip)               :: nelty
      
      !If not open, open the postfile
      call EndPrePos(a,itste,ttime)
      call Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
      call Mesh%GetNelem(nelem)
      
      call VTK_BeginScalarGP(a%VTKInterface,len_trim(wopos),trim(wopos))
      
      jelty = 1
      do ielty=1,nelty
         do ielem = 1,nelem
            call Mesh%GetIelty(ielem,jelty)
            if(jelty==ielty) then
               call VTK_WriteScalarGP(a%VTKInterface,ielem-1,real(bridge(ielem),rp))
            end if
         enddo
      enddo
      call VTK_EndScalarGP(a%VTKInterface)
      a%kfl_written = a%kfl_written + 1

   end subroutine pgpipsca

#endif  

end module Mod_postpr_VTK
