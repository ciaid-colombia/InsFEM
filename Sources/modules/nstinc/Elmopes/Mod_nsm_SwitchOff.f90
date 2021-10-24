  module Mod_nsm_SwitchOff
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersSwitchOff

   type, extends(PointerSetter) :: SPSwitchOff
contains
      procedure :: SpecificSet => SpecificSetSwitchOff
   end type
   type(SPSwitchOff) :: SetPointersSwitchOff
   
   integer(ip) :: ndofn,disc

contains
   
   subroutine SpecificSetSwitchOff(d)
      implicit none
      class(SPSwitchOff) :: d
      
      call ConcatenateProcedures(ProcHook_Initializations,ComputeSize)
      call ConcatenateProcedures(ProcHook_PreGauss,InitElem)
      call ConcatenateProcedures(ProcHook_InGauss,Switchoff)
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine ComputeSize
      implicit none
      ndofn = size(elmat(:,1,1,1))
   end subroutine

   subroutine InitElem
      implicit none
      disc = 0_ip
   end subroutine

   subroutine SwitchOff
      implicit none
      integer(ip) :: idofn,inode
      
      do inode = 1,e%pnode
         if (a%kfl_fixno(1,e%lnods(inode))==-3) then
            disc = 1_ip
            e%detjm = 0.0_rp
            !elmat(:,inode,:,:) = 0.0_rp
         endif
      enddo
      if (disc==1_ip) then
         do inode = 1,e%pnode
            do idofn = 1,ndofn   
               if (a%kfl_fixno(1,e%lnods(inode))==-3) then         
                  elmat(idofn,inode,idofn,inode) = 1.0_rp
               endif
            enddo 
         enddo
      endif 
   end subroutine
   
end module

