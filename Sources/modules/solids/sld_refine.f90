subroutine sld_Refine(a,itask)
   use typre
   use Mod_Element
   use Mod_phpRefineArrays
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   character(6) :: itask
   integer(ip)  :: oldnpoin,newnpoin,icomp,ndime,sz
   integer(ip)  :: iboun,newnboun,oldnboun
   real(rp), allocatable    :: auxdisp(:,:,:),auxveloc(:,:,:),auxaccel(:,:,:)
   real(rp), allocatable    :: aux_cp(:,:),auxbtrac(:,:)
   integer(ip), allocatable :: auxkfl_fixrs(:),auxkfl_bours(:)
   integer(ip), pointer     :: BoundaryNewToOld(:) => NULL(),iauxBoundaryMatch1(:) => NULL()
   type(i1p), pointer       :: iauxBoundaryMatch2(:) => NULL()

   !This subroutine transforms the tempe arrays from one mesh to the other
   !Adaptive Mesh Refinement
   !Old Dimensions
   oldnpoin = size(a%disp,2)
   oldnboun = size(a%kfl_bours)

   !We need to modify all the arrays
   !We assume that the mesh has already been updated
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNboun(newnboun)
   call a%Mesh%GetNdime(ndime)
   sz = (ndime*(ndime+1))/2

   !disp
   call a%Memor%alloc(ndime,newnpoin,a%ncomp,auxdisp,'disp','sld_Refine')
   do icomp = 1,a%ncomp
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(ndime,a%disp(:,:,icomp),auxdisp(:,:,icomp))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(ndime,a%disp(:,:,icomp),auxdisp(:,:,icomp))
      endif
   enddo
   call move_alloc(auxdisp,a%disp)
   call a%Mesh%InterpolateHangingValues(ndime,a%disp(:,:,1))
   call a%Memor%deallocObj(0,'disp','sld_Refine',rp*ndime*oldnpoin*a%ncomp)

   if (a%kfl_timei==1 ) then
       !veloc
       call a%Memor%alloc(ndime,newnpoin,2,auxveloc,'veloc','sld_Refine')
       do icomp = 1,2
           if (itask == 'Refine') then
               call a%Refiner%UpdateVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
           elseif (itask == 'Rebala') then
               call a%Refiner%RebalanceVariable(ndime,a%veloc(:,:,icomp),auxveloc(:,:,icomp))
           endif
       !call a%Mesh%InterpolateHangingValues(ndime,a%veloc(:,:,icomp))
       enddo
       call move_alloc(auxveloc,a%veloc)
       call a%Memor%deallocObj(0,'veloc','sld_Refine',rp*ndime*oldnpoin*2)

       !accel
       call a%Memor%alloc(ndime,newnpoin,2,auxaccel,'accel','sld_Refine')
       do icomp = 1,2
           if (itask == 'Refine') then
               call a%Refiner%UpdateVariable(ndime,a%accel(:,:,icomp),auxaccel(:,:,icomp))
           elseif (itask == 'Rebala') then
               call a%Refiner%RebalanceVariable(ndime,a%accel(:,:,icomp),auxaccel(:,:,icomp))
           endif
       !call a%Mesh%InterpolateHangingValues(ndime,a%accel(:,:,icomp))
       enddo
       call move_alloc(auxaccel,a%accel)
       call a%Memor%deallocObj(0,'accel','sld_Refine',rp*ndime*oldnpoin*2)
   endif

   if (a%kfl_docoupconv) then

       !disp_cp
       call a%Memor%alloc(ndime,newnpoin,aux_cp,'disp_cp','sld_Refine')
       if (itask == 'Refine') then
           call a%Refiner%UpdateVariable(ndime,a%disp_cp(:,:),aux_cp(:,:))
       elseif (itask == 'Rebala') then
           call a%Refiner%RebalanceVariable(ndime,a%disp_cp(:,:),aux_cp(:,:))
       endif
       call move_alloc(aux_cp,a%disp_cp)
       call a%Memor%deallocObj(0,'disp_cp','sld_Refine',rp*ndime*oldnpoin)

   endif

   !btraction
   call a%Memor%alloc(ndime,newnpoin,auxbtrac,'btrac','sld_Refine')
   if (itask == 'Refine') then
       call a%Refiner%UpdateVariable(ndime,a%btraction_nodal(:,:),auxbtrac(:,:))
   elseif (itask == 'Rebala') then
       call a%Refiner%RebalanceVariable(ndime,a%btraction_nodal(:,:),auxbtrac(:,:))
   endif
   call move_alloc(auxbtrac,a%btraction_nodal)
   call a%Memor%deallocObj(0,'btrac','sld_Refine',rp*ndime*oldnpoin)

   !kfl_fixrs
  ! call a%Memor%alloc(newnpoin,auxkfl_fixrs,'kfl_fixrs','sld_Refine')
  ! if (itask == 'Refine') then
  !    call a%Refiner%UpdateVariable(1_ip,a%kfl_fixrs,auxkfl_fixrs,'minab')
  ! elseif (itask == 'Rebala') then
  !    call a%Refiner%RebalanceVariable(1_ip,a%kfl_fixrs,auxkfl_fixrs)
  !endif
  !call move_alloc(auxkfl_fixrs,a%kfl_fixrs)
  !call a%Memor%deallocObj(0,'kfl_fixrs','sld_Refine',ip*oldnpoin)

  !if (itask == 'Refine') then
  !    call a%Mesh%GetBoundaryNewToOld(BoundaryNewToOld)
  !    !kfl_bours
  !    call a%Memor%alloc(newnboun,auxkfl_bours,'kfl_bours','sld_Refine')
  !    do iboun = 1,newnboun
  !        if (BoundaryNewToOld(iboun) > 0) then
  !            auxkfl_bours(iboun) = a%kfl_bours(BoundaryNewToOld(iboun))
  !        else
  !            auxkfl_bours(iboun) = 0
  !        endif
  !    enddo
  !    call move_alloc(auxkfl_bours,a%kfl_bours)
  !    call a%Memor%deallocObj(0,'kfl_bours','sld_Refine',ip*oldnboun)
  !elseif (itask == 'Rebala') then

  !    call a%Mesh%GetRebalanceInverseBoundaryMatches(iauxBoundaryMatch1,iauxBoundaryMatch2)
  !    call a%Memor%alloc(newnboun,auxkfl_bours,'kfl_bours','sld_Refine')

  !    call a%Refiner%RebalanceVariableBoundary(oldnboun,iauxBoundaryMatch1,iauxBoundaryMatch2,1,a%kfl_bours,auxkfl_bours)

  !    call move_alloc(auxkfl_bours,a%kfl_bours)
  !    call a%Memor%deallocObj(0,'kfl_bours','sld_Refine',oldnboun*ip)
  !endif

  !Setting the displacements to the mesh
  if(a%sld_type.eq.'NONLI') then
      call a%Mesh%SetALE(1_ip)
      call a%Mesh%SetDisplacements(a%disp)
  endif

  call a%SolidSpecificRefine(itask)

  call a%SolidSpecificRefineGauss(itask)


end subroutine

subroutine sld_refineGauss(a,itask)
   use typre
   use Mod_Element
   use Mod_Mesh
   use Mod_Solids
   use Mod_phpRefineArrays
   implicit none
   class(SolidsProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   character(6) :: itask
   integer(ip)  :: ndime,ielem,pgaus,sz,nelem,pnode

   call a%Mesh%GetNdime(ndime)
   sz = (ndime*(ndime+1))/2

   nelem=size(a%strain_g,1)

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
       deallocate(a%strain_g(ielem)%a)
       deallocate(a%extForce(ielem)%a)
       deallocate(a%pstrain(ielem)%a)
       deallocate(a%stress_g(ielem)%a)
       deallocate(a%btraction(ielem)%a)
       if(a%kfl_printPrincipalStresses) then
           deallocate(a%printPstrain(ielem)%a)
       end if
       a%strain_g(ielem)%a    => NULL()
       a%extForce(ielem)%a  => NULL()
       a%pstrain(ielem)%a   => NULL()
       a%stress_g(ielem)%a    => NULL()
       a%btraction(ielem)%a => NULL()
       if(a%kfl_printPrincipalStresses) then
           a%printPstrain(ielem)%a => NULL()
       end if
   end do

   if(a%kfl_printPrincipalStresses) then
       call a%Memor%dealloc(nelem,a%printPStrain,'printPStrain','sld_turnof')
   end if
   call a%Memor%dealloc(nelem,a%extForce,'extForce','sld_turnof')
   call a%Memor%dealloc(nelem,a%strain_g,'strain_g','sld_turnof')
   call a%Memor%dealloc(nelem,a%pStrain,'pStrain','sld_turnof')
   call a%Memor%dealloc(nelem,a%stress_g,'stress_g','sld_turnof')
   call a%Memor%dealloc(nelem,a%btraction,'btraction','sld_turnof')

   call a%Mesh%GetNelem(nelem)

   call a%Memor%alloc(nelem,a%extForce,'extForce','sld_memall')
   call a%Memor%alloc(nelem,a%stress_g,'stress_g','sld_memall')
   call a%Memor%alloc(nelem,a%strain_g,'strain_g','sld_memall')
   call a%Memor%alloc(nelem,a%pStrain,'pStrain','sld_memall')
   call a%Memor%alloc(nelem,a%btraction,'btraction','sld_memall')

   if(a%kfl_printPrincipalStresses) then
       call a%Memor%alloc(nelem,a%printPStrain,'printPStrain','sld_memall')
   end if

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_memall')

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

       call a%Mesh%ElementLoad(ielem,e)  
       allocate(a%extForce(ielem)%a(ndime,e%mnode))
       allocate(a%stress_g(ielem)%a(sz,pgaus))
       allocate(a%strain_g(ielem)%a(sz,pgaus))
       allocate(a%pstrain(ielem)%a(ndime))
       allocate(a%btraction(ielem)%a(ndime,pgaus))

       if(a%kfl_printPrincipalStresses) then
           allocate(a%printPstrain(ielem)%a(sz,pgaus))
           a%printPstrain(ielem)%a = 0.0_rp
       end if

       a%extForce(ielem)%a  = 0.0_rp
       a%strain_g(ielem)%a    = 0.0_rp
       a%pstrain(ielem)%a   = 0.0_rp
       a%stress_g(ielem)%a    = 0.0_rp
       a%btraction(ielem)%a = 0.0_rp

   end do

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_memall')

end subroutine
