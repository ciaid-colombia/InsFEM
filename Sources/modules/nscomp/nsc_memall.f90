subroutine nsc_memall(a)
  !-----------------------------------------------------------------------
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    compressible fluid flow equations
  !-----------------------------------------------------------------------
  use typre
  use Mod_Memor
  use Mod_Listen
  use Mod_Mesh
  use Mod_PhysicalProblem
  use Mod_NSCompressible
  implicit none
  
  class(NSCompressibleProblem) :: a
  
  integer(ip) ::  ndime,npoin,nelem,pgaus,ielem,pnode,nbody
  integer(ip) ::  coresgp_coun,moresgp_coun,enresgp_coun
  integer(ip) ::  cosgs_coun,mosgs_coun,ensgs_coun
  integer(ip) ::  scdiff_coun,sgdiff_coun
  
  call a%Mesh%GetNdime(ndime)
  call a%Mesh%GetNpoin(npoin)
  call a%Mesh%GetNelem(nelem)
  call a%Mesh%GetNbody(nbody)

  !Residual Projection
  if (a%kfl_repro >=1) then
     call a%Memor%alloc(ndime+2,npoin,a%repro,'repro','nsc_memall')
  endif

  !Gradient Orthogonal projection
  if(a%kfl_shock==2 .or. a%ErrorEstimatorTypeOfSubscales == 1) then
     call a%Memor%alloc(ndime+1,ndime,npoin,a%grprj,'grprj','nsc_memall')
  end if
  
  !Body forces and moments 
  if(a%kfl_outfm==1) then   
     if(a%npp_stepi(11) /= 0)then   
       call a%Memor%alloc(ndime,npoin,      a%forfl,'forfl','nsc_memall')
       call a%Memor%alloc(npoin,            a%hfxfl,'hfxfl','nsc_memall')
       if(ndime==2)then
          call a%Memor%alloc(1,npoin,a%momfl,'momfl','nsc_memall')
       elseif(ndime==3) then  
          call a%Memor%alloc(3,npoin,a%momfl,'momfl','nsc_memall')
       endif
     end if
     call a%Memor%alloc(nbody,a%lun_force,'lun_force','nsc_memall')
     call a%Memor%alloc(ndime+1,nbody,a%force,'force','nsc_memall')
     call a%Memor%alloc(3,nbody,a%formo,'formo','nsc_memall')
  end if

  !Postprocess the added shock capturing diffusiviy at the gauss points
  if (a%npp_stepi(9) /= 0) then
     scdiff_coun = 0
     call a%Memor%alloc(2_ip,nelem,a%scdiffarray,'scdiffarray','nsc_memall')
     do ielem=1,nelem
        call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
        allocate(a%scdiffarray(1,ielem)%a(pgaus))
        allocate(a%scdiffarray(2,ielem)%a(pgaus))
        a%scdiffarray(1,ielem)%a = 0.0_rp
        a%scdiffarray(2,ielem)%a = 0.0_rp
        scdiff_coun = scdiff_coun + 2*pgaus
     end do
     call a%Memor%allocObj(0,'scdiff%a','nsc_memall',scdiff_coun*rp)
  
  endif

  !Postprocess the added subgrid diffusiviy at the gauss points
  if (a%npp_stepi(10) /= 0) then
     sgdiff_coun = 0
     call a%Memor%alloc(2_ip,nelem,a%sgdiffarray,'sgdiffarray','nsc_memall')
     do ielem=1,nelem
        call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
        allocate(a%sgdiffarray(1,ielem)%a(pgaus))
        allocate(a%sgdiffarray(2,ielem)%a(pgaus))
        a%sgdiffarray(1,ielem)%a = 0.0_rp
        a%sgdiffarray(2,ielem)%a = 0.0_rp
        sgdiff_coun = sgdiff_coun + 2*pgaus
     end do
     call a%Memor%allocObj(0,'sgdiff%a','nsc_memall',sgdiff_coun*rp)
  
  endif

   !Subgrid scales
   if(a%kfl_trasg/=0) then
      cosgs_coun = 0
      mosgs_coun = 0
      ensgs_coun = 0
      !For explicit time integration scheme:
      !position 1: stage increment
      !position 2: stage contribution
      !position 3: previous time step 
      
      call a%Memor%alloc(nelem,a%cosgs,'cosgs','nsc_memall')
      call a%Memor%alloc(nelem,a%mosgs,'mosgs','nsc_memall')
      call a%Memor%alloc(nelem,a%ensgs,'ensgs','nsc_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%cosgs(ielem)%a(ncsgs,pgaus))
         allocate(a%mosgs(ielem)%a(ndime,ncsgs,pgaus))
         allocate(a%ensgs(ielem)%a(ncsgs,pgaus))
         a%cosgs(ielem)%a = 0.0_rp
         a%mosgs(ielem)%a = 0.0_rp
         a%ensgs(ielem)%a = 0.0_rp
         cosgs_coun = cosgs_coun + ncsgs*pgaus
         mosgs_coun = mosgs_coun + ndime*ncsgs*pgaus
         ensgs_coun = ensgs_coun + ncsgs*pgaus
      end do
      call a%Memor%allocObj(0,'cosgs%a','nsc_memall',cosgs_coun*rp)
      call a%Memor%allocObj(0,'mosgs%a','nsc_memall',mosgs_coun*rp)
      call a%Memor%allocObj(0,'ensgs%a','nsc_memall',ensgs_coun*rp)
   end if

   !Residual at the gauss points
   if ((a%kfl_shock == 1) .or. (a%npp_stepi(13) /= 0)) then
      coresgp_coun = 0
      moresgp_coun = 0
      enresgp_coun = 0
      call a%Memor%alloc(nelem,a%coResGP,'coResGP','nsc_memall')
      call a%Memor%alloc(nelem,a%moResGP,'moResGP','nsc_memall')
      call a%Memor%alloc(nelem,a%enResGP,'enResGP','nsc_memall')
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%coResGP(ielem)%a(pgaus))
         allocate(a%moResGP(ielem)%a(ndime,pgaus))
         allocate(a%enResGP(ielem)%a(pgaus))
         a%coResGP(ielem)%a = 0.0_rp
         a%moResGP(ielem)%a = 0.0_rp
         a%enResGP(ielem)%a = 0.0_rp
         coresgp_coun = coresgp_coun + pgaus
         moresgp_coun = moresgp_coun + ndime*pgaus
         enresgp_coun = enresgp_coun + pgaus
      end do
      call a%Memor%allocObj(0,'coresgp%a','nsc_memall',coresgp_coun*rp)
      call a%Memor%allocObj(0,'moresgp%a','nsc_memall',moresgp_coun*rp)
      call a%Memor%allocObj(0,'enresgp%a','nsc_memall',enresgp_coun*rp)
   
   endif

  !Statistics of fields
  if (a%npp_stepi(15) /= 0) then
     call a%Memor%alloc(npoin,    2,a%timad,'timad','nsc_memall')
     call a%Memor%alloc(ndime,npoin,2,a%timav,'timav','nsc_memall')
     call a%Memor%alloc(npoin,     2,a%timap,'timap','nsc_memall')
     call a%Memor%alloc(npoin,     2,a%timat,'timat','nsc_memall')
     call a%Memor%alloc(npoin,a%rmsd,'rmsd','nsc_memall')
     call a%Memor%alloc(ndime,npoin,a%rmsv,'rmsv','nsc_memall')
     call a%Memor%alloc(npoin,a%rmsp,'rmsp','nsc_memall')
     call a%Memor%alloc(npoin,a%rmst,'rmst','nsc_memall')
     call a%Memor%alloc(npoin,a%turbi,'turbi','nsc_memall')
  endif

  !FSI coupling
  call a%Memor%alloc(ndime,npoin,a%btraction,'btraction','nsc_memall')

  call a%SpecificNSCompMemall 

end subroutine nsc_memall 
