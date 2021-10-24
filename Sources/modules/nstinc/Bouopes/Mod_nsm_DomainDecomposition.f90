module Mod_nsm_DomainDecomposition
   use Mod_PointerSetter
   use Mod_nsm_BaseBouope
   implicit none
   private
   public SetPointersDomainD 

   type, extends(PointerSetter) :: SPDomainDecomposition
contains
      procedure :: SpecificSet => SpecificSetDomainDecomposition
   end type
   type(SPDomainDecomposition) :: SetPointersDomainD

   real(rp),    allocatable, dimension(:,:) :: eltrac
   real(rp),    allocatable, dimension(:)   :: gptrac

contains

   subroutine SpecificSetDomainDecomposition(d)
      implicit none
      class(SPDomainDecomposition) :: d
   
      call ConcatenateProcedures(ProcHook_Initializations,AllocateArraysDD)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocateArraysDD)
      call ConcatenateProcedures(ProcHook_Gathers,GathersDD)
      call ConcatenateProcedures(ProcHook_Interpolates,InterpolateDD)
      call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,AssemblyDD)
   end subroutine

   subroutine AllocateArraysDD
      implicit none
      call a%Memor%alloc(e%ndime,e%mnodb,eltrac,'eltrac','nsm_bouope')
      call a%Memor%alloc(e%ndime,gptrac,'gptrac','nsm_bouope')
   end subroutine

   subroutine DeallocateArraysDD
      implicit none
      call a%Memor%dealloc(e%ndime,e%mnodb,eltrac,'eltrac','nsm_bouope')
      call a%Memor%dealloc(e%ndime,gptrac,'gptrac','nsm_bouope')
   end subroutine

   subroutine GathersDD
      implicit none
      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         if (associated(a%etraction)) call e%gatherb(e%ndime,eltrac,a%etraction(:,:))
      end if
   end subroutine
   
   subroutine InterpolateDD
      implicit none
      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) call e%interpb(e%ndime,eltrac,gptrac)
   end subroutine
   
   subroutine AssemblyDD
      implicit none
      integer(ip) :: idime,inode,inodb
      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            elrhs(1:e%ndime,inode)=elrhs(1:e%ndime,inode) -&
               e%shapb(inodb,e%igaub)*dsurf*gptrac(1:e%ndime)
         end do
      end if
   end subroutine

end module
