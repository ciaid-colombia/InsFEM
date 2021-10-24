module Mod_nsm_BoundarySubscales
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersBoundarySubscales

   type, extends(PointerSetter) :: SPBoundarySubscales
contains
      procedure :: SpecificSet => SpecificSetBoundarySubscales
   end type
   type(SPBoundarySubscales) :: SetPointersBoundarySubscales

   procedure(), pointer :: ProcPointer_OBSGS => NULL()
   
   real(rp), allocatable :: ngradn(:),gpreb(:),elreb(:,:),grveb(:,:)
   real(rp), allocatable :: wreprb(:,:)
   real(rp)              :: temom,gpprb
   
contains

   subroutine SpecificSetBoundarySubscales(d)
      implicit none
      class(SPBoundarySubscales) :: d
   
      !ResidualProjection
      if (a%kfl_bousg /= 0) then
         call ConcatenateProcedures(ProcHook_Initializations,AllocBoundarySGS)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocBoundarySGS)
         call ConcatenateProcedures(ProcPointer_PostGaussElmats,PostGaussElmatsBoundary)
         if (a%kfl_repro == 1) then
            call ConcatenateProcedures(ProcHook_Initializations,AllocRep)
            call ConcatenateProcedures(ProcHook_Finalizations,DeallocRep)
            
            call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroRep)
            
            !We need to compute the residual at the gauss point
            call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpRes)
            
            !Now we assembly everything for computing the projection
            call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,GpresToElres)
            
            call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyResidual)
            call ConcatenateProcedures(ProcHook_PostLoop,ProjectResidual)
         endif
      endif
   end subroutine

   subroutine AllocBoundarySGS
      implicit none
      
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      
      call a%Memor%alloc(ndime,ngradn,'ngradn','nsm_BoundarySubscales')
   end subroutine

   subroutine DeallocBoundarySGS
      implicit none
      
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      
      call a%Memor%dealloc(ndime,ngradn,'ngradn','nsm_BoundarySubscales')
   end subroutine

  !For Boundary SGS
   subroutine PostGaussElmatsBoundary
      implicit none
      real(rp)    :: dsurf,dsurf0,ngradn(e%pnodb)
      integer(ip) :: igaub,iface,inodb
     
      call ComputeTauBoundary(e,chale,acvis,temom) 
      !call ComputeNGradN(e,ngradN) 
      do iface = 1,e%pface
         e%iface = iface
         call a%Mesh%FaceLoad(iface,ielem,e)
         dsurf0 = 0.0_rp
         do igaub = 1,e%pgaub
            e%igaub = igaub
            call e%elmderb
            call e%elenor
            dsurf=e%weigb(e%igaub)*e%eucta
            dsurf0 = dsurf0 + dsurf

            !Compute n.gradN
            !do inodb=1,e%pnodb
            !   ngradN(inodb) = dot_product(e%baloc(:,e%ndime),e%cartb(:,inodb))
            !end do
      
            !Compute contributions to LHS:
            !call nsm_elmbuv_bou(e,dsurf,temom,acvis,ngradN,wrmat1)
            !call nsm_elmbpv_bou(e,dsurf,temom,acvis,ngradN,elmpv)
            !call nsm_elmbuq_bou(e,dsurf,temom,acvis,ngradN,elmuq)
            !call nsm_elmbpq_bou(e,dsurf,temom,elmpq)

            !Orthogonal SGS
            call ProcPointer_OBSGS
         end do
      end do

   end subroutine
   
  !For Orthogonal Boundary SGS
   subroutine AllocRep
      implicit none
      integer(ip) :: ndime,npoin
      
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(e%ndime,npoin,wreprb,'wreprb','nsm_EndElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elreb,'elreb','nsm_BoundarySGS')
      call a%Memor%alloc(e%ndime,gpreb,'gpreb','nsm_BoundarySGS')
      call a%Memor%alloc(e%ndime,e%ndime,grveb,'grveb','nsm_BoundarySGS')
   end subroutine
   
   subroutine DeAllocRep
      implicit none

      call a%Memor%dealloc(e%ndime,e%mnode,elreb,'elreb','nsm_BoundarySGS')
      call a%Memor%dealloc(e%ndime,gpreb,'gpreb','nsm_BoundarySGS')
      call a%Memor%dealloc(e%ndime,e%ndime,grveb,'grveb','nsm_BoundarySGS')
   end subroutine

   subroutine ComputeGpRes

   end subroutine ComputeGpRes
 
   subroutine ElmatsToZeroRep
      implicit none
      
      elreb = 0.0_rp
   end subroutine

   subroutine GpresToElres
      implicit none
      integer(ip) :: inode
   
   !   do inode = 1,e%pnode
   !      elreb(:,inode) = elreb(:,inode) + e%shape(inode,e%igaus)*gpreb*dvol
   !   enddo
   end subroutine

   subroutine AssemblyResidual
      implicit none
      
      call a%Mesh%AssemblyToArray(e,size(wreprb,1),elreb,wreprb) 
   end subroutine

   subroutine ProjectResidual
      implicit none
      integer(ip) :: ndime,npoin,ipoin
      
      if (a%kfl_SwitchOff==1) then
         call a%Project_disc(size(wreprb,1),wreprb,a%kfl_fixno(1,:)) 
      else
         call a%Project(size(wreprb,1),wreprb) 
      endif
      a%repro = wreprb
      call move_alloc(wreprb,a%repro)
      !Dealloc wreprb
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%deallocObj(0,'wreprb','nsm_EndElmope',e%ndime*npoin*rp)
   end subroutine

   !This suborutine computes the residual at a Gauss point
   subroutine ComputeBoundaryResidual
      implicit none
      integer(ip) :: idime,jdime
      
      gpreb = 0.0_rp
      gpreb(1:e%ndime) = gpreb(1:e%ndime) + e%baloc(1:e%ndime,e%ndime)*gpprb
      do idime = 1,e%ndime
         gpreb(1:e%ndime) = gpreb(1:e%ndime) + e%baloc(jdime,e%ndime)*acvis*grveb(:,jdime) 
      end do
   end subroutine ComputeBoundaryResidual
end module
