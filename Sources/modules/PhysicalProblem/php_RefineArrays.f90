module Mod_phpRefineArrays
   use typre
   use Mod_PhysicalProblem
   implicit none
   private
   public php_RefineArrays, php_RefineArrays_b
   
   interface php_RefineArrays
      module procedure php_RefineArraysReal1, php_RefineArraysReal2, php_RefineArraysReal0
   end interface
   
   interface php_RefineArrays_b
      module procedure php_RefineArraysReal1b, php_RefineArraysReal2b
   end interface
   
contains

   subroutine php_RefineArraysReal0(a,itask,array,vanam)
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: itask
      real(rp), allocatable :: array(:)
      character(*) :: vanam
      
      integer(ip) :: icomp,newnpoin,oldnpoin
      
      real(rp), allocatable :: auxarray(:)
      
      call a%Mesh%GetNpoin(newnpoin)
      oldnpoin = size(array,1)


      
      !array
      call a%Memor%alloc(newnpoin,auxarray,vanam,'php_Refine')
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(1_ip,array,auxarray)
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(1_ip,array,auxarray)
         endif   
      call move_alloc(auxarray,array)
      call a%Memor%deallocObj(0,vanam,'php_Refine',rp*oldnpoin)
   end subroutine


   subroutine php_RefineArraysReal1(a,itask,array,vanam)
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: itask
      real(rp), allocatable :: array(:,:)
      character(*) :: vanam
      
      integer(ip) :: ncomp,icomp,newnpoin,oldnpoin
      
      real(rp), allocatable :: auxarray(:,:)
      
      ncomp = size(array,2)
      
      call a%Mesh%GetNpoin(newnpoin)
      oldnpoin = size(array,1)


      
      !array
      call a%Memor%alloc(newnpoin,ncomp,auxarray,vanam,'php_Refine')
      do icomp = 1,ncomp
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(1_ip,array(:,icomp),auxarray(:,icomp))
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(1_ip,array(:,icomp),auxarray(:,icomp))
         endif   
      enddo
      call move_alloc(auxarray,array)
      call a%Memor%deallocObj(0,vanam,'php_Refine',rp*oldnpoin*ncomp)
   end subroutine


   subroutine php_RefineArraysReal2(a,itask,array,vanam)
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: itask
      real(rp), allocatable :: array(:,:,:)
      character(*) :: vanam
      
      real(rp), allocatable :: auxarray(:,:,:)
      integer(ip) :: ncomp,icomp,oldnpoin,newnpoin,ndime
      
      
      
      call a%Mesh%GetNpoin(newnpoin)
      ndime = size(array,1)
      oldnpoin = size(array,2)
      ncomp = size(array,3)
      
      !array
      call a%Memor%alloc(ndime,newnpoin,ncomp,auxarray,vanam,'php_Refine')
      do icomp = 1,ncomp
         if (itask == 'Refine') then
            call a%Refiner%UpdateVariable(ndime,array(:,:,icomp),auxarray(:,:,icomp))
         elseif (itask == 'Rebala') then
            call a%Refiner%RebalanceVariable(ndime,array(:,:,icomp),auxarray(:,:,icomp))
         endif   
      enddo
      call move_alloc(auxarray,array)
      call a%Memor%deallocObj(0,vanam,'php_Refine',rp*ndime*oldnpoin*ncomp)
   end subroutine
   
   
   
   subroutine php_RefineArraysReal1b(a,itask,array,vanam)
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: itask
      real(rp), allocatable :: array(:,:)
      character(*) :: vanam
      
      integer(ip) :: ncomp,icomp,newnpoin,oldnpoin
      
      real(rp), allocatable :: auxarray(:,:)
      integer(ip) :: ndime
      
      ndime = size(array,1)
      
      call a%Mesh%GetNpoin(newnpoin)
      oldnpoin = size(array,2)


      
      !array
      call a%Memor%alloc(ndime,newnpoin,auxarray,vanam,'php_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(ndime,array,auxarray)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(ndime,array,auxarray)
      endif   
      call move_alloc(auxarray,array)
      call a%Memor%deallocObj(0,vanam,'php_Refine',rp*oldnpoin*ndime)
   end subroutine
   
   subroutine php_RefineArraysReal2b(a,itask,array,vanam)
      implicit none
      class(PhysicalProblem) :: a
      character(6) :: itask
      real(rp), allocatable, target :: array(:,:,:)
      character(*) :: vanam
      integer(ip) :: ncomp,icomp,newnpoin,oldnpoin
      real(rp), allocatable, target :: auxarray(:,:,:)
      
      real(rp), pointer :: pointerArray(:,:) => NULL(), pointerAuxArray(:,:) => NULL()
      integer(ip) :: ndime1,ndime2
      
      ndime1 = size(array,1)
      ndime2 = size(array,2)
      
      call a%Mesh%GetNpoin(newnpoin)
      oldnpoin = size(array,3)
      
      
      call a%Memor%alloc(ndime1,ndime2,newnpoin,auxarray,vanam,'php_Refine')
      
      pointerArray(1:ndime1*ndime2,1:oldnpoin) => array
      pointerAuxArray(1:ndime1*ndime2,1:newnpoin) => auxarray
      
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(ndime1*ndime2,Pointerarray,Pointerauxarray)
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(ndime1*ndime2,Pointerarray,Pointerauxarray)
      endif   
      call move_alloc(auxarray,array)
      call a%Memor%deallocObj(0,vanam,'php_Refine',rp*oldnpoin*ndime1*ndime2)
   end subroutine
   
end module

