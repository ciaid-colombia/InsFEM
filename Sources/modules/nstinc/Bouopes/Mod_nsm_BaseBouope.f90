module Mod_nsm_BaseBouope
   use typre
   use Mod_nsm_BaseElmope
   use Mod_NavierStokesBoundary
   implicit none
   
   integer(ip) :: igaub,iboun,nboun,ibody,nbody
   real(rp)    :: weightfactor
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: gpbpr(1)
   real(rp)    :: prgau(1),lambda=0.0_rp
   real(rp)    :: vista=0.0_rp,beta=0.0_rp
      
   real(rp),    allocatable, dimension(:,:)     :: bovel,gpbve,grbve,bosig
   real(rp),    allocatable, dimension(:)       :: gpcod,bopre,tract,vmassBoundary
   real(rp)   , allocatable, dimension(:)       :: sigau,eltem

  !Iterfaces for sldsup base procedures
  interface

      module subroutine AllocateBaseBouopeArraysSUP
      end subroutine 

      module subroutine DeallocateBaseBouopeArraysSUP 
      end subroutine 

      module subroutine BoundaryInterpolatesSUP
      end subroutine 

  end interface

contains

   !*********************!
   ! General subroutines !
   !*********************!

   subroutine AllocateBaseBouopeArrays
      implicit none
      call a%Memor%alloc(e%ndime,1,gpbve,'gpbve','nsm_bouope')
      call a%Memor%alloc(e%ndime,e%mnode,1,elvel,'elvel','nsm_bouope')
      call a%Memor%alloc(e%ndime,e%mnodb,bovel,'bovel','nsm_bouope')
      call a%Memor%alloc(e%ndime,e%ndime,grbve,'grbve','nsm_bouope')
      call a%Memor%alloc(e%mnodb,bopre,'bopre','nsm_bouope')
      call a%Memor%alloc(e%ndime,tract,'tract','nsm_bouope')   
      call a%Memor%alloc(e%ndime,gpcod,'gpcod','nsm_bouope')
   end subroutine

   subroutine AllocateBaseBouopeMatrices
      implicit none
      call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsm_bouope')
      call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','nsm_bouope')
   end subroutine

   subroutine DeallocateBaseBouopeArrays
      implicit none
      call a%Memor%dealloc(e%ndime,1,gpbve,'gpbve','nsm_bouope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,elvel,'elvel','nsm_bouope')
      call a%Memor%dealloc(e%ndime,e%mnodb,bovel,'bovel','nsm_bouope')
      call a%Memor%dealloc(e%ndime,e%ndime,grbve,'grbve','nsm_bouope')
      call a%Memor%dealloc(e%mnodb,bopre,'bopre','nsm_bouope')
      call a%Memor%dealloc(e%ndime,tract,'tract','nsm_bouope')   
      call a%Memor%dealloc(e%ndime,gpcod,'gpcod','nsm_bouope')
   end subroutine

   subroutine DeallocateBaseBouopeMatrices
      implicit none
      call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsm_bouope')
      call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','nsm_bouope')
   end subroutine

   subroutine BouopeElMatsToZero
      implicit none
      elmat=0.0_rp
      elrhs=0.0_rp
   end subroutine

   subroutine BoundaryGathers
      implicit none
      call e%gather(e%ndime,elvel,a%veloc(:,:,1))
      call e%gatherb(e%ndime,bovel,a%veloc(:,:,1))
      call e%gatherb(1_ip   ,bopre,a%press(:,1))       
   end subroutine

   subroutine BoundaryInterpolates
      implicit none
      call e%interpb(e%ndime,e%bocod,gpcod)
      call e%interpb(e%ndime,bovel,gpbve(:,1))
      call e%interpb(1,bopre,gpbpr)
      call e%gradientb(e%ndime,elvel,grbve)
   end subroutine

end module

#include "../../sigmaup/Bouopes/Mod_sup_BaseBouope.f90"   
