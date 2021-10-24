 module Mod_nsm_ComputeTauSmoothing
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersComputeTauSmoothing

   type, extends(PointerSetter) :: SPComputeTauSmoothing
contains
      procedure :: SpecificSet => SpecificSetComputeTauSmoothing
   end type
   type(SPComputeTauSmoothing) :: SetPointersComputeTauSmoothing
   
   real(rp), allocatable :: elTauSmoothing(:,:), elSmoothedTau2(:,:)
   real(rp), allocatable :: wTausmo(:,:)
   real(rp) :: gptau(2)
   
contains
   
   subroutine SpecificSetComputeTauSmoothing(d)
      implicit none
      class(SPComputeTauSmoothing) :: d

      call ConcatenateProcedures(ProcHook_Initializations,AllocTauSmoothing)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocTauSmoothing)
      call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroTauSmoothing)
      call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,GpTauToElTauSmoothing)
      call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyTausmo)
      call ConcatenateProcedures(ProcHook_PostLoop,SmoothTausmo)
   end subroutine   
   
   !---------------------------------------------------------------------------
    subroutine AllocTauSmoothing
      implicit none
      integer(ip) :: npoin
      call a%Memor%alloc(2,e%mnode,elTauSmoothing,'elTauSmoothing','nsm_EndElmope')
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(2,npoin,wTausmo,'wTausmo','nsm_EndElmope')
   end subroutine
   
   subroutine DeallocTauSmoothing
      implicit none
      integer(ip) :: npoin
      call a%Memor%dealloc(2,e%mnode,elTauSmoothing,'elTauSmoothing','nsm_EndElmope')
   end subroutine
   
   subroutine ElmatsToZeroTauSmoothing
      implicit none
      elTauSmoothing = 0.0_rp
   end subroutine
   
   subroutine GpTauToElTauSmoothing
      implicit none
      real(rp)    :: timomaux,tidivaux
      integer(ip) :: inode
      
      !We compute the Taus at the Gauss point
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timomaux)
      call ComputeTauDiv(e,a%staco(4),chale,timomaux,tidivaux)
      !And we assembly them
      do inode = 1,e%pnode
         elTauSmoothing(1,inode) = elTauSmoothing(1,inode) + e%shape(inode,e%igaus)*timomaux*dvol
         elTauSmoothing(2,inode) = elTauSmoothing(2,inode) + e%shape(inode,e%igaus)*tidivaux*dvol
      enddo
   end subroutine
   
   subroutine AssemblyTausmo
      implicit none
      !wTausmo(:,e%lnods(1:e%pnode)) = wTausmo(:,e%lnods(1:e%pnode)) + elTauSmoothing(:,1:e%pnode)
      call a%Mesh%AssemblyToArray(e,size(wTausmo,1),elTauSmoothing,wtausmo) 
   end subroutine
   
   subroutine SmoothTausmo
      implicit none
      integer(ip) :: npoin, inode,ismooth
      class(FiniteElement), pointer :: e => NULL()
      
      call a%Mesh%Smooth(2,wTausmo) 
      call move_alloc(wTausmo,a%Tausmo)
      !deallocate(wTausmo)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%deallocObj(0,'wTausmo','nsm_EndElmope',2*npoin*rp)
      
      if (a%kfl_Tausm > 1) then
         call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supm_EnditeElmope')
         call a%Memor%alloc(2,e%mnode,elSmoothedTau2,'elSmoothedTau2','nsm_EndElmope')
         call a%Memor%alloc(2,e%mnode,elTauSmoothing,'elTauSmoothing','nsm_EndElmope')
         call a%Mesh%GetNelem(nelem)
         do ismooth = 2,a%kfl_Tausm
            call a%Mesh%GetNpoin(npoin)
            call a%Memor%alloc(2,npoin,wTausmo,'wTausmo','nsm_EndElmope')
            do ielem = 1,nelem
               call a%Mesh%ElementLoad(ielem,e)
               call e%elmdcg
               call e%gather(2,elSmoothedTau2,a%Tausmo)
               elTauSmoothing = 0.0_rp
               do igaus = 1,e%pgaus
                  e%igaus = igaus
                  dvol = e%weigp(e%igaus)*e%detjm
                  call e%interpg(2,elSmoothedTau2,gptau)
                  !And we assembly them
                  do inode = 1,e%pnode
                     elTauSmoothing(:,inode) = elTauSmoothing(:,inode) + e%shape(inode,e%igaus)*gptau(1:2)*dvol
                  enddo
               end do
               call a%Mesh%AssemblyToArray(e,size(wTausmo,1),elTauSmoothing,wTausmo)
            end do
            a%Tausmo=0.0_rp
            call a%Mesh%Smooth(2,wTausmo)
            call move_alloc(wTausmo,a%Tausmo)
            !deallocate(wTausmo)
            call a%Mesh%GetNpoin(npoin)
            call a%Memor%deallocObj(0,'wTausmo','supm_EndElmope',2*npoin*rp)
         enddo
         call a%Memor%dealloc(2,e%mnode,elTauSmoothing,'elTauSmoothing','nsm_EndElmope')
         call a%Memor%dealloc(2,e%mnode,elSmoothedTau2,'elSmoothedTau2','nsm_EndElmope')
         call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supm_EnditeElmope')
      endif
   end subroutine

end module
