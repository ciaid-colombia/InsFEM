module Mod_ParallelSystemInterface
   use typre
   use Mod_Memor
   use Mod_ParallelOrderingInterface
   use Mod_Element
   implicit none
   private
   public ParallelSystemInterface
   
   type ParallelSystemInterface
   
      integer(ip) :: ndofn,npoinLocal,gnpoin,npoinGhost
      type(MemoryMan), pointer   :: Memor => NULL()
      integer(ip) :: kfl_flush = 0   !Default is do not flush
      integer(ip) :: MPIsize, MPIrank, MPIroot,MPIcomm
   
contains
      procedure :: SetMPI
      procedure :: SetFlush
      procedure :: Init
      procedure :: Tozero
      procedure :: ToZeroMatrix
      procedure :: TozeroRhs
      procedure :: Solve
      procedure :: Assembly
      procedure :: AssemblyRowBlock
      procedure :: AssemblyElmat
      procedure :: AssemblyPointRhs
      procedure :: AssemblyEigenMat
      procedure :: AssemblyElrhs
      procedure :: AssemblyForceGhosts
      procedure :: Deallocate
      procedure :: MatGetDiagonal
      procedure :: CopyToRHS
      procedure :: Matvec
      procedure :: RhsDot
      procedure :: Init2

   end type
   
contains

   subroutine SetMPI(a,comm,size,root,irank)
      implicit none
      class(ParallelSystemInterface) :: a
      integer(ip) :: irank,size,root,comm
      
      a%MPIcomm = comm
      a%MPIrank = irank
      a%MPIsize = size
      a%MPIroot = root
   end subroutine
   
   subroutine SetFlush(a,kfl_flush)
      use typre
      implicit none
      class(ParallelSystemInterface) :: a
      integer(ip) :: kfl_flush
      
      a%kfl_flush = kfl_flush
   end subroutine
      

   subroutine Init(a,npoinLocal,npoinGhost,gnpoin,ndofn,issym,ia,ja,LocalOrdering,Memor,optionsfile,exmod,lun_solve,opti)
      use typre
      use Mod_Memor
      use Mod_ParallelOrderingInterface
      implicit none
      class(ParallelSystemInterface) :: a
      class(ParallelOrderingInterface) :: LocalOrdering
      type(MemoryMan), target   :: Memor
      integer(ip) :: gnpoin,npoinLocal,npoinGhost,ndofn,issym,ia(*),ja(*)
      character(150) :: exmod,optionsfile
      integer(ip) :: lun_solve
      character(5), optional :: opti
      
      call runend('ParallelSystem, Init procedure not implemented')
      
   end subroutine  

   subroutine ToZero(a)
      implicit none
      class(ParallelSystemInterface) :: a
      
      call runend('ParallelSystem, ToZero procedure not implemented')
      
   end subroutine
   
   subroutine ToZeroMatrix(a)
      implicit none
      class(ParallelSystemInterface) :: a
      
      call runend('ParallelSystem, ToZero procedure not implemented')
      
   end subroutine
   
   subroutine ToZeroRhs(a)
      implicit none
      class(ParallelSystemInterface) :: a
      
      call runend('ParallelSystem, ToZeroRhs procedure not implemented')
      
   end subroutine
   
   subroutine MatGetDiagonal(a,ipoin,diag,iwa,rwa)
      use typre
      implicit none
      class(ParallelSystemInterface) :: a
      integer(ip) :: ipoin
      real(rp) :: diag(a%ndofn)
      integer(ip) :: iwa(*)
      real(rp)    :: rwa(*)

      call runend('ParallelSystem, MatGetDiagonal procedure not implemented')
   end subroutine   

   subroutine Solve(a,unkno)
      implicit none
      class(ParallelSystemInterface) :: a
      real(rp) :: unkno(a%ndofn,a%npoinLocal)

      call runend('ParallelSystem, Solve procedure not implemented')

   end subroutine  

   subroutine Assembly(a,e,elmat,elrhs)
      implicit none
      class(ParallelSystemInterface) :: a
      class(FiniteElement) :: e
      real(rp), target :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)

      call runend('ParallelSystem, Assembly procedure not implemented')

   end subroutine
   
   subroutine AssemblyRowBlock(a,e,irowBlock,elmat,elrhs)
      implicit none
      class(ParallelSystemInterface) :: a
      class(FiniteElement) :: e
      real(rp), target :: elmat(a%ndofn,1,a%ndofn,e%mnode), elrhs(a%ndofn,1)
      integer(ip) :: irowBlock

      call runend('ParallelSystem, AssemblyRowBlock procedure not implemented')

   end subroutine
   
   subroutine AssemblyElmat(a,e,elmat)
      implicit none
      class(ParallelSystemInterface) :: a
      class(FiniteElement) :: e
      real(rp), target :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode)

      call runend('ParallelSystem, AssemblyElmat procedure not implemented')

   end subroutine
   
   subroutine AssemblyPointRhs(a,rhSize,pointId,pointRhs)
      implicit none
      class(ParallelSystemInterface) :: a
      integer(ip) :: rhSize,pointId(rhSize)
      real(rp):: pointRhs(a%ndofn,rhSize)
      
      call runend('ParallelSystem, AssemblyPointRhs procedure not implemented')
      
   end subroutine

   subroutine AssemblyElrhs(a,e,elrhs)
      implicit none
      class(ParallelSystemInterface) :: a
      class(FiniteElement) :: e
      real(rp), target :: elrhs(a%ndofn,e%mnode)
      
      call runend('ParallelSystem, AssemblyElrhs procedure not implemented')
      
   end subroutine
   
   subroutine AssemblyForceGhosts(a,e,elmat,elrhs)
      use Mod_Element
      implicit none
      class(ParallelSystemInterface) :: a
      class(FiniteElement) :: e
      real(rp), target :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)
      
      call runend('ParallelSystem, Assembly procedure not implemented')
      
   end subroutine  
   
   subroutine CopyToRHS(a,rhsvec)
      class(ParallelSystemInterface) :: a
      real(rp) :: rhsvec(a%ndofn,a%npoinLocal)
      
      call runend('ParallelSystem, CopyToRHS procedure not implemented')
      
   end subroutine
   
   subroutine Deallocate(a)
      implicit none
      class(ParallelSystemInterface) :: a
      
      call runend('ParallelSystem, Deallocate procedure not implemented')
      
   end subroutine  
   
   subroutine Matvec(a,vec,vecout)
      use typre
      implicit none
      class(ParallelSystemInterface) :: a
      real(rp) :: vec(a%ndofn,a%npoinLocal), vecout(a%ndofn,a%npoinLocal)
      
      call runend('ParallelSystem, Matvec not implemented')
   end subroutine
   
   subroutine RhsDot(a,vec,vdot)
      use typre
      implicit none
      class(ParallelSystemInterface) :: a
      real(rp) :: vec(a%ndofn,a%npoinLocal), vdot
      
      call runend('ParallelSystem, RhsDot not implemented')
   end subroutine

   subroutine Init2(a,kfl_massMatrix,kfl_Precondition) !For projected LS
      use typre
      implicit none
      class(ParallelSystemInterface) :: a
      integer(ip)  :: kfl_Precondition    !Use a Precondition projection
      logical      :: kfl_massMatrix      !Mass matrix projection in basis
      
      call runend('ParallelSystem, Init2 procedure not implemented')
      
   end subroutine  

   subroutine AssemblyEigenMat(a,e,massmat,kmat)
      implicit none
      class(ParallelSystemInterface) :: a
      type(FiniteElement) :: e
      real(rp), target :: massmat(a%ndofn,e%mnode,a%ndofn,e%mnode)
      real(rp), target :: kmat(a%ndofn,e%mnode,a%ndofn,e%mnode)
    
      call runend('ParallelSystem, AssemblyEigenMat procedure not implemented')
    
   end subroutine

end module
