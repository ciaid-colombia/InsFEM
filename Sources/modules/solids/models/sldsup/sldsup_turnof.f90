subroutine sldsup_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_SUPSolids
   implicit none
   class(SUPSolidsProblem) :: a
   integer(ip)             :: ndime,npoin,ielem,nelem,ncomp
   integer(ip)             :: sz,pgaus,pnode,resid_coun
   integer(ip)             :: usgs_coun,ssgs_coun,psgs_coun,ncsgs
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   ncomp = a%ncomp
   sz = (ndime*(ndime+1))/2

   call a%Memor%dealloc(sz   ,npoin,a%ncomp,a%sigma,'sigma','sldsup_turnof')

   call a%Memor%dealloc(sz   ,npoin  ,a%devstrain,'devstrain','sldsup_turnof')

   if(a%kfl_docoupconv) then 
       call a%Memor%dealloc(sz   ,npoin,a%sigma_cp,'sigma_cp','sldsup_turnof')
   endif

   if(a%kfl_printJ2Stresses) then 
       call a%Memor%dealloc(npoin  ,a%j2,'J2','sldsup_turnof')
   endif

   call a%SUPSolidSpecificTurnof

end subroutine sldsup_turnof

subroutine sldsup_Specificturnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_SUPSolids
   implicit none
   class(SUPSolidsProblem) :: a
   integer(ip)             :: ndime,npoin,ielem,nelem,ncomp
   integer(ip)             :: sz,pgaus,pnode,resid_coun
   integer(ip)             :: usgs_coun,ssgs_coun,psgs_coun,ncsgs
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   ncomp = a%ncomp
   sz = (ndime*(ndime+1))/2

   !Residual Projection
   if (a%kfl_repro >=1) then
      call a%Memor%dealloc(ndime+sz+1,npoin,a%repro,'repro','lmn_turnof')
   endif

   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      usgs_coun = 0
      ssgs_coun = 0
      psgs_coun = 0
      ncsgs     = 2               !Number of components in time for accuracy
      if(a%kfl_tacsg>0) ncsgs=a%ncomp-1
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%u_sgs(ielem)%a)
         deallocate(a%s_sgs(ielem)%a)
         deallocate(a%p_sgs(ielem)%a)
         usgs_coun = usgs_coun + ndime*ncsgs*pgaus
         ssgs_coun = ssgs_coun + sz*pgaus
         psgs_coun = psgs_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'u_sgs%a','sldsup_turnof',usgs_coun*rp)
      call a%Memor%deallocObj(0,'s_sgs%a','sldsup_turnof',ssgs_coun*rp)
      call a%Memor%deallocObj(0,'p_sgs%a','sldsup_turnof',psgs_coun*rp)
      call a%Memor%dealloc(nelem,a%u_sgs,'u_sgs','sldsup_turnof')
      call a%Memor%dealloc(nelem,a%s_sgs,'s_sgs','sldsup_turnof')
      call a%Memor%dealloc(nelem,a%p_sgs,'p_sgs','sldsup_turnof')
   end if

   !Postprocess the residual at the gauss points
   if (a%npp_stepi(16) /= 0) then
      resid_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%residualU(ielem)%a)
         deallocate(a%residualP(ielem)%a)
         deallocate(a%residualS(ielem)%a)
         resid_coun = resid_coun + (ndime+sz+1)*pgaus
      end do
      call a%Memor%deallocObj(0,'residual%a','lmn_turnof',resid_coun*rp)
      call a%Memor%dealloc(nelem,a%residualU,'residual','lmn_turnof')
      call a%Memor%dealloc(nelem,a%residualP,'residual','lmn_turnof')
      call a%Memor%dealloc(nelem,a%residualS,'residual','lmn_turnof')
   endif

end subroutine sldsup_Specificturnof
