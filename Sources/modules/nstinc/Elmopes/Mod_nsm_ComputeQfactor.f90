  module Mod_nsm_ComputeQfactor
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   implicit none
   private
   public SetPointersComputeQFactor

   type, extends(PointerSetter) :: SPComputeQFactor
contains
      procedure :: SpecificSet => SpecificSetComputeQFactor
   end type
   type(SPComputeQFactor) :: SetPointersComputeQFactor

   real(rp), allocatable :: elqfac(:,:)
   real(rp):: gpqfac(1)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SpecificSetComputeQFactor(d)
      implicit none
      class(SPComputeQFactor) :: d
      logical :: aux_logic
      
      !Logics for deciding if we need to compute Q-factor
      aux_logic = .false.
      if (a%npp_stepi(21) /= 0) then   
         if (mod(a%istep,a%npp_stepi(21))==0) then
            aux_logic = .true.
         endif
      endif
      
      !If Q-factor needs to be computed
      if (aux_logic .eqv. .true.) then
         call SetPointersInterpolateGradients%Set
         call ConcatenateProcedures(ProcHook_Initializations,AllocQfactor)
         call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,ComputeQfactor)
         call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyQfactor)
         call ConcatenateProcedures(ProcHook_Finalizations,DeallocQfactor)
      endif  
   end subroutine   
   
   subroutine AllocQfactor
      implicit none
      
      call a%Memor%alloc(1,e%mnode,elqfac,'elqfac','nsm_EndElmope')
      a%qfac(:) = 0_rp
   end subroutine
   
   subroutine ComputeQfactor
      implicit none
      integer(ip)  :: inode,idime,jdime
      real(rp) :: gpdef

      gpdef = 0.0_rp
      do idime=1,e%ndime
         do jdime=1,e%ndime
            if (idime/=jdime) then
               gpdef = gpdef + grvel(jdime,idime)*grvel(idime,jdime)
            endif
         enddo
      enddo 

      gpqfac(1) = 0.5_rp*(divvel**2_ip - gpdef)
           
      do inode=1,e%pnode
         elqfac(1,inode)=elqfac(1,inode) + e%shape(inode,e%igaus)*gpqfac(1)*dvol
      end do
   end subroutine

   subroutine AssemblyQfactor
      implicit none
      integer(ip)  :: inode
      
      call a%Mesh%AssemblyToArray(e,1_ip,elqfac,a%qfac) 
      elqfac(:,:)=0_rp
   end subroutine

   subroutine DeallocQfactor
      implicit none

      call a%Memor%dealloc(1,e%mnode,elqfac,'elqfac','nsm_EndElmope')
   end subroutine
end module 
