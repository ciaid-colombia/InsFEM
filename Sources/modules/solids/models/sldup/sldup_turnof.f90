subroutine sldup_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem) :: a
   integer(ip)             :: ndime,npoin,ielem,nelem,ncomp
   integer(ip)             :: sz,pgaus,pnode,resid_coun
   integer(ip)             :: usgs_coun,ssgs_coun,psgs_coun,ncsgs
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   ncomp = a%ncomp

   call a%Memor%dealloc(npoin,      a%ncomp,a%press,'press','sldup_memall')

   if(a%kfl_docoupconv) then 
       call a%Memor%dealloc(npoin      ,a%press_cp,'press_cp','sldup_memall')
   endif

   call a%UPSolidSpecificTurnof
   call a%TurnofExtend

end subroutine sldup_turnof

subroutine sldup_Specificturnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem) :: a
   integer(ip)             :: ndime,npoin,ielem,nelem,ncomp
   integer(ip)             :: sz,pgaus,pnode,resid_coun
   integer(ip)             :: usgs_coun,ssgs_coun,psgs_coun,ncsgs
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   sz = (ndime*(ndime+1))/2

   ncomp = a%ncomp

   if(a%kfl_printNodalSigma) call a%Memor%dealloc(sz   ,npoin  ,1,a%sigma,'sigma','sldup_turnof')

   !Residual Projection
   if (a%kfl_repro >=1) then
      call a%Memor%dealloc(ndime+1,npoin,a%repro,'repro','sldup_turnof')
   endif

   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      usgs_coun = 0
      psgs_coun = 0
      ncsgs     = 2               !Number of components in time for accuracy
      if(a%kfl_tacsg>0) ncsgs=a%ncomp-1
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%u_sgs(ielem)%a)
         deallocate(a%p_sgs(ielem)%a)
         usgs_coun = usgs_coun + ndime*ncsgs*pgaus
         psgs_coun = psgs_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'u_sgs%a','sldup_turnof',usgs_coun*rp)
      call a%Memor%deallocObj(0,'p_sgs%a','sldup_turnof',psgs_coun*rp)
      call a%Memor%dealloc(nelem,a%u_sgs,'u_sgs','sldup_turnof')
      call a%Memor%dealloc(nelem,a%p_sgs,'p_sgs','sldup_turnof')
   end if

   !Postprocess the residual at the gauss points
   if (a%npp_stepi(16) /= 0) then
      resid_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%residualU(ielem)%a)
         deallocate(a%residualP(ielem)%a)
         resid_coun = resid_coun + (ndime+1)*pgaus
      end do
      call a%Memor%deallocObj(0,'residual%a','sldup_turnof',resid_coun*rp)
      call a%Memor%dealloc(nelem,a%residualU,'residual','sldup_turnof')
      call a%Memor%dealloc(nelem,a%residualP,'residual','sldup_turnof')
   endif

   if(a%kfl_printJ2Stresses) then 

       do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

       deallocate(a%j2_g(ielem)%a)
       end do

       call a%Memor%dealloc(nelem,a%j2_g,'J2','sldup_memall')
   endif

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
       if (a%kfl_printGaussSigma) then
           deallocate(a%sigma_g(ielem)%a)
       endif

   end do
   call a%Memor%dealloc(nelem,a%sigma_g,'sigma_g','sldup_turnof')

end subroutine sldup_Specificturnof
