subroutine sldup_refine(a,itask)
   use typre
   use Mod_Element
   use Mod_phpRefineArrays
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   character(6) :: itask

   integer(ip) :: oldnpoin,newnpoin,icomp,ndime
   real(rp), allocatable    :: auxstrain(:,:) ,auxpress(:,:)

   !This subroutine transforms the tempe arrays from one mesh to the other
   !Adaptive Mesh Refinement
   !Old Dimensions
   oldnpoin = size(a%disp,2)

   !We need to modify all the arrays
   !We assume that the mesh has already been updated
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNdime(ndime)

   !press
   call a%Memor%alloc(newnpoin,a%ncomp,auxpress,'press','sld_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%press(:,icomp),auxpress(:,icomp))
      endif
   enddo
   call move_alloc(auxpress,a%press)
   call a%Mesh%InterpolateHangingValues(1_ip,a%press(:,1))
   call a%Memor%deallocObj(0,'press','sld_Refine',rp*oldnpoin*a%ncomp)

   !strain
   !call a%Memor%alloc(sz,newnpoin,auxstrain,'strain','sld_Refine')
   !   if (itask == 'Refine') then
   !      call a%Refiner%UpdateVariable(sz,a%strain(:,:),auxstrain(:,:))
   !   elseif (itask == 'Rebala') then
   !      call a%Refiner%RebalanceVariable(sz,a%strain(:,:),auxstrain(:,:))
   !   endif
   !call move_alloc(auxstrain,a%strain)
   !call a%Mesh%InterpolateHangingValues(sz,a%strain(:,:))
   !call a%Memor%deallocObj(0,'strain','sld_Refine',rp*sz*oldnpoin)

end subroutine

subroutine sldup_refineGauss(a,itask)
   use typre
   use Mod_Element
   use Mod_Mesh
   use Mod_UPSolids
   use Mod_phpRefineArrays
   implicit none
   class(UPSolidsProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   character(6) :: itask
   integer(ip) :: ndime,ielem,pgaus,nelem,pnode

   nelem=size(a%btraction,1)

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
       deallocate(a%extForce(ielem)%a)
       deallocate(a%btraction(ielem)%a)
       a%btraction(ielem)%a => NULL()
       a%extForce(ielem)%a  => NULL()
   end do
   call a%Memor%dealloc(nelem,a%btraction,'btraction','sld_turnof')
   call a%Memor%dealloc(nelem,a%extForce,'extForce','sld_turnof')

   call a%Mesh%GetNelem(nelem)

   call a%Memor%alloc(nelem,a%btraction,'btraction','sld_memall')

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_memall')
   call a%Memor%alloc(nelem,a%extForce,'extForce','sld_memall')

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

       call a%Mesh%ElementLoad(ielem,e)  
       allocate(a%extForce(ielem)%a(ndime,e%mnode))
       allocate(a%btraction(ielem)%a(ndime,pgaus))

       a%btraction(ielem)%a = 0.0_rp
       a%extForce(ielem)%a  = 0.0_rp

   end do

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_memall')

end subroutine
