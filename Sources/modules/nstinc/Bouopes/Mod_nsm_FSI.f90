module Mod_nsm_FSI
   use Mod_PointerSetter
   use Mod_nsm_BaseBouope
   implicit none
   private
   public SetPointersFSI

   type, extends(PointerSetter) :: SPFSI
contains
      procedure :: SpecificSet => SpecificSetFSI
   end type
   type(SPFSI) :: SetPointersFSI

   real(rp)                                 :: alfa_inv=0.0_rp
   real(rp),    allocatable, dimension(:,:) :: etraction

contains

   subroutine SpecificSetFSI(d)
      implicit none
      class(SPFSI) :: d
            
      if(associated(a%etraction) .and. a%doRobin .eqv. .true.) then
         call ConcatenateProcedures(ProcHook_Gathers,GatherFSI)
         call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,AssembleFSI)
      endif
   end subroutine

   subroutine GatherFSI
      alfa_inv=0.0_rp
      if(a%alfa_robin < zensi) then
         if(a%MPIrank == a%MPIroot) write(*,*) &
         'Warning alpha robin = 0, imposing standard Dirichlet B.C on Fluid'
         ! Reset flag so it does not add traction terms
         a%doRobin = .false.
      else
         alfa_inv = 1.0_rp/a%alfa_robin
         etraction = 0.0_rp
         call e%gatherb(e%ndime,etraction,a%etraction(:,:))
      endif
   end subroutine
   
   subroutine AssembleFSI
      implicit none
      integer(ip) :: idime,inode,inodb

      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         elrhs(1:e%ndime,inode)=elrhs(1:e%ndime,inode) -&
            e%shapb(inodb,e%igaub)*dsurf*alfa_inv*etraction(1:e%ndime,inodb)
      end do
   end subroutine

end module
