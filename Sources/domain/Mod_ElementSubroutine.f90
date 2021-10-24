module Mod_ElementLoop
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none


      class(FiniteElement), pointer :: e => NULL()
      
      class(FemMesh), pointer    :: Mesh => NULL()
      type(MemoryMan), pointer   :: Memor => NULL()

      
      integer(ip) :: ndofn = 1
      real(rp), allocatable :: elmat(:,:,:,:), elrhs(:,:)
   
      procedure(), pointer :: Initializations  => NULLsub
      procedure(), pointer :: ElmatsToZero     => NULLsub
      procedure(), pointer :: Gathers          => NULLsub
      procedure(), pointer :: PreGauss         => NULLsub
      procedure(), pointer :: Interpolates     => NULLsub
      procedure(), pointer :: InGauss          => NULLsub
      procedure(), pointer :: InGaussElmats    => NULLsub
      procedure(), pointer :: PostGauss        => NULLsub
      procedure(), pointer :: PostGaussElmats  => NULLsub
      procedure(), pointer :: PreDirichlet     => NULLsub
      procedure(), pointer :: Dirichlet        => NULLsub
      procedure(), pointer :: AssemblySystem   => NULLsub
      procedure(), pointer :: Finalizations    => NULLsub
      
      procedure(), pointer :: GaussHOOK   => NULLsub
      procedure(), pointer :: ElementHook    => NULLsub
             

contains

   subroutine SetMesh(eMesh)
      class(FemMesh), target :: eMesh
      
      Mesh => eMesh
   end subroutine
      
   subroutine SetMemor(eMemor)
      type(MemoryMan), target :: eMemor
      
      Memor => eMemor
   end subroutine
   
   subroutine SetNdofn(endofn)
      
      integer(ip) :: endofn
      
      ndofn = endofn
   end subroutine
   
   subroutine NULLSUB
      
   end subroutine
   
   subroutine SetPointersToNullSub
      Initializations  => NULLsub
      ElmatsToZero     => NULLsub
      Gathers          => NULLsub
      PreGauss         => NULLsub
      Interpolates     => NULLsub
      InGauss          => NULLsub
      InGaussElmats    => NULLsub
      PostGauss        => NULLsub
      PostGaussElmats  => NULLsub
      PreDirichlet     => NULLsub
      Dirichlet        => NULLsub
      AssemblySystem   => NULLsub
      Finalizations    => NULLsub
      
      ElementHook      => NULLSUB
      GAussHOok        => NULLSUB
      
   end subroutine
      
   subroutine ElementLoop
      integer(ip) :: ielem, nelem

      !Allocate  element
      call Mesh%ElementAlloc(e,Memor,'DefaultRule','nsm_EnditeElmope')
      !Allocate Matrices just in case it they are necessary
      call Memor%alloc(ndofn,e%mnode,ndofn,e%mnode,elmat,'elmat','AllocateMatrices')
      call Memor%alloc(ndofn,e%mnode,elrhs,'elrhs','AllocateMatrices')
      
      
      call Initializations
      
      call Mesh%GetNelem(nelem)
      elements : do ielem = 1,nelem
         !Load Element
         call Mesh%ElementLoad(ielem,e)  
         
         call ElementHook
      enddo
      
      call Finalizations
      
      !Allocate Matrices just in case it they are necessary
      call Memor%dealloc(ndofn,e%mnode,ndofn,e%mnode,elmat,'elmat','AllocateMatrices')
      call Memor%dealloc(ndofn,e%mnode,elrhs,'elrhs','AllocateMatrices')
      !Allocate  element
      call Mesh%ElementDeAlloc(e,Memor,'DefaultRule','nsm_EnditeElmope')
   end subroutine
   
   
   subroutine GaussLoop
      integer(ip) :: igaus
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         call GaussHook
      enddo gauss_points
   end subroutine
   
   end subroutine
         
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
         !Element length at center of gravity
         call e%elmlen
                  
         call ElmatsToZero
         
         call Gathers
         
         call PreGauss
   
         !Gauss Point Loop
         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            
            call e%elmder
            
            call Interpolates
            
            call InGauss
            
            call InGaussElmats
         enddo gauss_points
         
         call PostGauss
         
         call PostGaussElmats
         
         call PreDirichlet
         
         call Dirichlet
         
         call AssemblySystem
      enddo elements
      
      
   
   end subroutine
   
  


end module

