subroutine nsc_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the compressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   implicit none
   class(NSCompressibleProblem) :: a
   
   integer(ip) :: ndime,npoin,nelem,pgaus,ielem,pnode,nbody
   integer(ip) :: coresgp_coun,moresgp_coun,enresgp_coun
   integer(ip) :: cosgs_coun,mosgs_coun,ensgs_coun
   integer(ip) :: scdiff_coun,sgdiff_coun,ibody
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNbody(nbody)
   
   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outph,'(//,a,/)') '     * * * END OF COMPRESSIBLE NAVIER-STOKES RUN * * *'
      write(a%lun_solve,'(//,a,/)') '     * * * END OF COMPRESSIBLE NAVIER-STOKES RUN * * *'
   endif
   
   !Residual Projection
   if (a%kfl_repro >=1) then
      call a%Memor%dealloc(ndime+2,npoin,a%repro,'repro','nsc_turnof')
   endif

   !Gradient Orthogonal projection
   if(a%kfl_shock==2 .or. a%ErrorEstimatorTypeOfSubscales == 1) then
      call a%Memor%dealloc(ndime,ndime+1,npoin,a%grprj,'grprj','nsc_turnof')
   end if

   !Body forces and moments 
   if(a%kfl_outfm==1)then   
      if(a%npp_stepi(11) /= 0)then   
        call a%Memor%dealloc(ndime,npoin,      a%forfl,'forfl','nsc_turnof')
        call a%Memor%dealloc(npoin,            a%hfxfl,'hfxfl','nsc_turnof')
        if(ndime==2)then
           call a%Memor%dealloc(1,npoin,a%momfl,'momfl','nsc_turnof')
        elseif(ndime==3) then 
           call a%Memor%dealloc(3,npoin,a%momfl,'momfl','nsc_turnof')
        endif
      end if
      call a%Memor%dealloc(ndime+1,nbody,a%force,'force','nsc_turnof')
      call a%Memor%dealloc(3,nbody,a%formo,'formo','nsc_turnof')
      do ibody=1,nbody
         close(a%lun_force(ibody))
      end do
      call a%Memor%dealloc(nbody,a%lun_force,'lun_force','lmn_turnof')
   end if

   !Postprocess the added shock capturing diffusiviy at the gauss points
   if (a%npp_stepi(9) /= 0) then
      scdiff_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%scdiffarray(1,ielem)%a)
         deallocate(a%scdiffarray(2,ielem)%a)
         scdiff_coun = scdiff_coun + 2*pgaus
      end do
      call a%Memor%deallocObj(0,'scdiff%a','nsc_turnof',scdiff_coun*rp)
      call a%Memor%dealloc(2_ip,nelem,a%scdiffarray,'scdiffarray','nsc_turnof')
   
   endif

   !Postprocess the added subgrid diffusiviy at the gauss points
   if (a%npp_stepi(10) /= 0) then
      sgdiff_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%sgdiffarray(1,ielem)%a)
         deallocate(a%sgdiffarray(2,ielem)%a)
         sgdiff_coun = sgdiff_coun + 2*pgaus
      end do
      call a%Memor%deallocObj(0,'sgdiff%a','nsc_turnof',sgdiff_coun*rp)
      call a%Memor%dealloc(2_ip,nelem,a%sgdiffarray,'sgdiffarray','nsc_turnof')
   
   endif

   !Subgrid scales
   if(a%kfl_trasg/=0) then
      cosgs_coun = 0
      mosgs_coun = 0
      ensgs_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%cosgs(ielem)%a)
         deallocate(a%mosgs(ielem)%a)
         deallocate(a%ensgs(ielem)%a)
         cosgs_coun = cosgs_coun + ncsgs*pgaus
         mosgs_coun = mosgs_coun + ndime*ncsgs*pgaus
         ensgs_coun = ensgs_coun + ncsgs*pgaus
      end do
      call a%Memor%deallocObj(0,'cosgs%a','nsc_turnof',cosgs_coun*rp)
      call a%Memor%deallocObj(0,'mosgs%a','nsc_turnof',mosgs_coun*rp)
      call a%Memor%deallocObj(0,'ensgs%a','nsc_turnof',ensgs_coun*rp)
      call a%Memor%dealloc(nelem,a%cosgs,'cosgs','nsc_turnof')
      call a%Memor%dealloc(nelem,a%mosgs,'mosgs','nsc_turnof')
      call a%Memor%dealloc(nelem,a%ensgs,'ensgs','nsc_turnof')
   end if
   
   !Residual at the gauss points
   if ((a%kfl_shock == 1) .or. (a%npp_stepi(13) /= 0)) then
      coresgp_coun = 0
      moresgp_coun = 0
      enresgp_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%coResGP(ielem)%a)
         deallocate(a%moResGP(ielem)%a)
         deallocate(a%enResGP(ielem)%a)
         coresgp_coun = coresgp_coun + pgaus
         moresgp_coun = moresgp_coun + ndime*pgaus
         enresgp_coun = enresgp_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'coresgp%a','nsc_turnof',coresgp_coun*rp)
      call a%Memor%deallocObj(0,'moresgp%a','nsc_turnof',moresgp_coun*rp)
      call a%Memor%deallocObj(0,'enresgp%a','nsc_turnof',enresgp_coun*rp)
      call a%Memor%dealloc(nelem,a%coResGP,'coResGP','nsc_turnof')
      call a%Memor%dealloc(nelem,a%moResGP,'moResGP','nsc_turnof')
      call a%Memor%dealloc(nelem,a%enResGP,'enResGP','nsc_turnof')
   endif

  !Statistics of fields
   if (a%npp_stepi(15) /= 0) then
      call a%Memor%dealloc(npoin,    2,a%timad,'timad','nsc_turnof')
      call a%Memor%dealloc(ndime,npoin,2,a%timav,'timav','nsc_turnof')
      call a%Memor%dealloc(npoin,     2,a%timap,'timap','nsc_turnof')
      call a%Memor%dealloc(npoin,     2,a%timat,'timat','nsc_turnof')
      call a%Memor%dealloc(npoin,a%rmsd,'rmsd','nsc_turnof')
      call a%Memor%dealloc(ndime,npoin,a%rmsv,'rmsv','nsc_turnof')
      call a%Memor%dealloc(npoin,a%rmsp,'rmsp','nsc_turnof')
      call a%Memor%dealloc(npoin,a%rmst,'rmst','nsc_turnof')
      call a%Memor%dealloc(npoin,a%turbi,'turbi','nsc_turnof')
   endif

   !FSI coupling
   call a%Memor%dealloc(ndime,npoin,a%btraction,'btraction','nsc_turnof')

   !NRBCS
   if (a%kfl_nscnbc == 1) then
      call a%Memor%dealloc(npoin,a%bvmass,'BVMASS','nsc_turnof')
   endif

   call a%SpecificNSCompTurnof

end subroutine nsc_turnof

