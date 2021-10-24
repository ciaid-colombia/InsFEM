subroutine php_Refine(a,itask)
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   character(6) :: itask
   
   integer(ip) :: newnpoinLocal,newnpoin,newgnpoin,newnelem,lnodsSize,newnboun
   
   integer(ip) :: oldnpoin,oldnboun, nptra
   
   integer(ip), allocatable :: auxkfl_fixno(:,:),auxkfl_funno(:)
   real(rp), allocatable    :: auxunkno(:,:),auxbvess(:,:,:)
   integer(ip), allocatable :: auxkfl_fixbo(:),auxkfl_funbo(:)
   type(r1p), allocatable :: auxbvnat(:)
   
   integer(ip) :: ibvess, nbvess,bvnatsize,iboun,AllocatedBytes,oldboun
   integer(ip), pointer :: BoundaryNewToOld(:) => NULL()
   
   integer(ip), pointer :: iauxBoundaryMatch1(:) => NULL()
   type(i1p), pointer   :: iauxBoundaryMatch2(:) => NULL()
   integer(ip), allocatable :: aux_kfl_fixbo(:)
   
   integer(ip) :: ibopo, isHanging, idofn,ipoin
   
   call a%Timer%Refine%Tic
   
   !Adaptive Mesh Refinement
   !Old Dimensions
   oldnpoin = size(a%unkno,2)
   oldnboun = size(a%kfl_fixbo)
   
   !We need to modify all the arrays
   !We assume that the mesh has already been updated
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNboun(newnboun)

   !kfl_Fixno
   call a%Memor%alloc(a%ndofbc,newnpoin,auxkfl_fixno,'kfl_fixno','php_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(a%ndofbc,a%kfl_fixno,auxkfl_fixno,'minim')
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(a%ndofbc,a%kfl_fixno,auxkfl_fixno)
   endif
   call move_alloc(auxkfl_fixno,a%kfl_fixno)
   call a%Memor%deallocObj(0,'kfl_fixno','php_Refine',ip*a%ndofbc*oldnpoin)
   
   !Additional checks in order to avoid fixno in the interior of the domain
   if (a%kfl_SwitchOff==0) then
      do ipoin = 1,newnpoin 
         call a%Mesh%GetIbopo(ipoin,ibopo)
         call a%Mesh%IsHangingNode(ipoin,isHanging)
         !Not a hanging node, and not in the boundary
         if ((ibopo == 0) .and. (isHanging == 0)) then
            do idofn = 1,a%ndofbc
               if (a%kfl_fixno(idofn,ipoin) > 0) a%kfl_fixno(idofn,ipoin) = 0
            enddo
         endif
      enddo
   endif
         
   !Bvess
   nbvess = size(a%bvess,3)
   call a%Memor%alloc(a%ndofbc,newnpoin,nbvess,auxbvess,'bvess','php_Refine')
   do ibvess = 1,nbvess
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(a%ndofbc,a%bvess(:,:,ibvess),auxbvess(:,:,ibvess))
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(a%ndofbc,a%bvess(:,:,ibvess),auxbvess(:,:,ibvess))
      endif
   enddo
   call move_alloc(auxbvess,a%bvess)
   call a%Memor%deallocObj(0,'bvess','php_Refine',rp*a%ndofbc*oldnpoin*nbvess)
   
   !Fixbo
   if (itask == 'Refine') then
      call a%Mesh%GetBoundaryNewToOld(BoundaryNewToOld)
      !kfl_fixbo
      call a%Memor%alloc(newnboun,auxkfl_fixbo,'kfl_fixbo','php_Refine')
      do iboun = 1,newnboun
         if (BoundaryNewToOld(iboun) > 0) then
            auxkfl_fixbo(iboun) = a%kfl_fixbo(BoundaryNewToOld(iboun))
         else
            auxkfl_fixbo(iboun) = 0
         endif
      enddo
      call move_alloc(auxkfl_fixbo,a%kfl_fixbo)
      call a%Memor%deallocObj(0,'kfl_fixbo','php_Refine',ip*oldnboun)
      
      !bvnat
      call a%Memor%alloc(newnboun,auxbvnat,'bvnat','php_Refine')
      bvnatsize = 0
      do iboun = 1,newnboun
         if (a%kfl_fixbo(iboun) /= 0 ) then
            oldboun = BoundaryNewToOld(iboun)
            if (oldboun /= 0) then
               if (associated(a%bvnat(oldboun)%a)) then
                  bvnatsize = bvnatsize+size(a%bvnat(BoundaryNewToOld(iboun))%a,1)
                  allocate(auxbvnat(iboun)%a(size(a%bvnat(BoundaryNewToOld(iboun))%a,1)))
                  auxbvnat(iboun)%a = a%bvnat(BoundaryNewToOld(iboun))%a
               endif
            endif
         endif
      enddo
      call a%Memor%allocObj(0,'bvnat%a','php_Refine',rp*bvnatsize)
      !dealloc old one
      bvnatsize = 0
      do iboun = 1,oldnboun         
         if (associated(a%bvnat(iboun)%a)) then
            bvnatsize = bvnatsize+size(a%bvnat(iboun)%a)
            deallocate(a%bvnat(iboun)%a)
         endif
      enddo
      call a%Memor%deallocObj(0,'bvnat%a','php_Refine',rp*bvnatsize)
      !Copy
      call move_alloc(auxbvnat,a%bvnat)
      call a%Memor%deallocObj(0,'bvnat','php_Refine',rp*oldnboun)
   elseif (itask == 'Rebala') then
      
      call a%Mesh%GetRebalanceInverseBoundaryMatches(iauxBoundaryMatch1,iauxBoundaryMatch2)
      call a%Memor%alloc(newnboun,aux_kfl_fixbo,'kfl_fixbo','php_Refine')
      call a%Refiner%RebalanceVariableBoundary(oldnboun,iauxBoundaryMatch1,iauxBoundaryMatch2,1,a%kfl_fixbo,aux_kfl_fixbo)
      call move_alloc(aux_kfl_fixbo,a%kfl_fixbo)
      call a%Memor%deallocObj(0,'kfl_fixbo','php_Refine',ip*oldnboun)

      !Now for bvnat
      call a%Memor%alloc(newnboun,auxbvnat,'bvnat','php_Refine')
      call a%Refiner%RebalanceVariableBoundary(oldnboun,iauxBoundaryMatch1,iauxBoundaryMatch2,a%bvnat,auxbvnat,AllocatedBytes)

!       a%kfl_fixbo = 0
      call a%Memor%allocObj(0,'bvnat%a','php_Refine',AllocatedBytes)
      bvnatsize = 0
      do iboun = 1,oldnboun         
         if (associated(a%bvnat(iboun)%a)) then
            bvnatsize = bvnatsize+size(a%bvnat(iboun)%a)
            deallocate(a%bvnat(iboun)%a)
         endif
      enddo
      call a%Memor%deallocObj(0,'bvnat%a','php_Refine',rp*bvnatsize)
      !Copy
      call move_alloc(auxbvnat,a%bvnat)
      call a%Memor%deallocObj(0,'bvnat','php_Refine',rp*oldnboun)
   endif

   
   !Non constant boundary conditions
   if (a%kfl_conbc /= 1) then
      !Funno 
      call a%Memor%alloc(newnpoin,auxkfl_funno,'kfl_funno','php_Refine')
      if (itask == 'Refine') then
         call a%Refiner%UpdateVariable(1_ip,a%kfl_funno,auxkfl_funno,'minim')
      elseif (itask == 'Rebala') then
         call a%Refiner%RebalanceVariable(1_ip,a%kfl_funno,auxkfl_funno)
      endif
      call move_alloc(auxkfl_funno,a%kfl_funno)
      call a%Memor%deallocObj(0,'kfl_funno','php_Refine',ip*oldnpoin)
      
      !Funbo
      !This is wrong since it needs to be coded!!!!!
      if (itask == 'Refine') then
         call a%Mesh%GetBoundaryNewToOld(BoundaryNewToOld)
         !kfl_fixbo
         call a%Memor%alloc(newnboun,auxkfl_funbo,'kfl_funbo','php_Refine')
         do iboun = 1,newnboun
            if (BoundaryNewToOld(iboun) > 0) then
               auxkfl_funbo(iboun) = a%kfl_funbo(BoundaryNewToOld(iboun))
            else
               auxkfl_funbo(iboun) = 0
            endif
         enddo
         call move_alloc(auxkfl_funbo,a%kfl_funbo)
         call a%Memor%deallocObj(0,'kfl_funbo','php_Refine',ip*oldnboun)
      elseif (itask == 'Rebala') then
         call a%Mesh%GetRebalanceInverseBoundaryMatches(iauxBoundaryMatch1,iauxBoundaryMatch2)
         call a%Memor%alloc(newnboun,aux_kfl_fixbo,'kfl_funbo','php_Refine')
         
         call a%Refiner%RebalanceVariableBoundary(oldnboun,iauxBoundaryMatch1,iauxBoundaryMatch2,1,a%kfl_funbo,aux_kfl_fixbo)
         
         call move_alloc(aux_kfl_fixbo,a%kfl_funbo)
         call a%Memor%deallocObj(0,'kfl_funbo','php_Refine',oldnboun*ip)
      endif
   endif
  
   !Reallocate the linear system
   if (a%kfl_SkipLinearSystemRefinement == 0) then
      call a%LinearSystemTurnof
      call a%LinearSystemMemall
   endif
   
   !Specific actions for each module
   call a%SpecificRefine(itask)
   
   !Tracking of points
   if(a%nptra>0) then 
      nptra = 0
      if (a%MPIrank == a%MPIroot) nptra = a%nptra
      call a%TrackingInterpolator%Finalize
      call a%TrackingInterpolator%Initialize_Interp(a%Mesh,a%cptra(:,1:nptra))
   endif
   
   call a%Timer%Refine%Toc
   
end subroutine
