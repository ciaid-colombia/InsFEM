submodule (Mod_nsm_BaseBouope) Mod_sup_BaseBouope

   implicit none

contains

   !*********************!
   ! General subroutines !
   !*********************!

   module subroutine AllocateBaseBouopeArraysSUP
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call a%Memor%alloc(tn,e%mnodb,bosig,'bosig','sup_forces')
      call a%Memor%alloc(tn,sigau,'sigau','sup_forces')
      call a%Memor%alloc(e%mnode,eltem,'eltem','sup_forces')
   end subroutine

   module subroutine DeallocateBaseBouopeArraysSUP
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call a%Memor%dealloc(tn,e%mnodb,bosig,'bosig','sup_forces')
      call a%Memor%dealloc(tn,sigau,'sigau','sup_forces')
      call a%Memor%dealloc(e%mnode,eltem,'eltem','sup_forces')
  end subroutine

   module subroutine BoundaryInterpolatesSUP
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call e%interpb(1,bopre,prgau(:))
      call e%interpb(tn,bosig,sigau)
   end subroutine

end submodule Mod_sup_BaseBouope
