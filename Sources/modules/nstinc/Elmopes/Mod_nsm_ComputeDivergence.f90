  module Mod_nsm_ComputeDivergence
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   implicit none
   private
   public SetPointersComputeDivergence

   type, extends(PointerSetter) :: SPComputeDivergence
contains
      procedure :: SpecificSet => SpecificSetComputeDivergence
   end type
   type(SPComputeDivergence) :: SetPointersComputeDivergence
   
   real(rp), allocatable :: elvort(:,:), gpvort(:)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SpecificSetComputeDivergence(d)
      implicit none
      class(SPComputeDivergence) :: d
      logical :: aux_logic
      
      !Logics for deciding if we need to compute Divergence
      aux_logic = .false.
      if (a%npp_stepi(11) /= 0) then   
         if (mod(a%istep,a%npp_stepi(11))==0) then
            aux_logic = .true.
         endif
      endif
      
      !If Divergence needs to be computed
      !We prepend it so that it is computed before the element cuts
      if (aux_logic .eqv. .true.) then
         call PrependProcedure(ProcHook_PreGauss,ComputeDivergence)
      endif  
   end subroutine   
   
   !Divergence 
   subroutine ComputeDivergence
      implicit none
      integer(ip)  :: idime,igaus
      
      !We recompute everything so that it works also in the case of cut elements
      call e%elmdel
      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         call e%gradient(e%ndime,elvel(:,:,1),grvel)    !Vel. gradient
         !Velocity divergence
         a%Divergence(ielem)%a(e%igaus) = 0.0_rp
         do idime = 1,e%ndime
            a%Divergence(ielem)%a(e%igaus) = a%Divergence(ielem)%a(e%igaus) + grvel(idime,idime)
         enddo
      enddo
   end subroutine

end module


 
