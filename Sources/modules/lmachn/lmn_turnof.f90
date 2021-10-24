subroutine lmn_turnof(a)
   ! DESCRIPTION
   !    This routine closes the run for the incompressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   
   integer(ip) :: ndime,npoin,nelem,ncomp,ncsgs,pgaus,iaux,ielem,pnode,nboun,nbody,vesgs_coun,tesgs_coun,prsgs_coun,resid_coun,ibody
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNbody(nbody)
   ncomp = a%ncomp
   
   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outph,'(//,a,/)') '     * * * END OF LOW MACH RUN * * *'
      write(a%lun_solve,'(//,a,/)') '     * * * END OF LOW MACH RUN * * *'
   endif
   call a%Memor%dealloc(ndime,npoin,a%ncomp,a%veloc,'veloc','lmn_turnof')
   call a%Memor%dealloc(npoin,a%ncomp      ,a%press,'press','lmn_turnof')
   call a%Memor%dealloc(      a%ncomp      ,a%pther,'pther','lmn_turnof')
   call a%Memor%dealloc(npoin,a%ncomp      ,a%tempe,'tempe','lmn_turnof')
   call a%Memor%dealloc(npoin              ,a%itemp,'itemp','lmn_turnof')
   if(a%npp_stepi(4)>0) then     
      call a%Memor%dealloc(npoin              ,a%densf,'densf','lmn_turnof')
   end if   

   if (a%kfl_sourc == 2) then
      call a%Memor%dealloc(npoin,a%PointwiseSource,'PointwiseSource','lmn_reabcs')
   endif      
   
   !Transient subgrid scales
   if(a%kfl_trasg/=0) then
      ncsgs=1
      vesgs_coun = 0
      tesgs_coun = 0
      prsgs_coun = 0
      if(a%kfl_tacsg>0) ncsgs=a%ncomp-1
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%vesgs(ielem)%a)
         deallocate(a%tesgs(ielem)%a)
         deallocate(a%prsgs(ielem)%a)
         vesgs_coun = vesgs_coun + ndime*ncsgs*pgaus
         tesgs_coun = tesgs_coun + ncsgs*pgaus
         prsgs_coun = prsgs_coun + pgaus
      end do
      call a%Memor%deallocObj(0,'vesgs%a','lmn_turnof',vesgs_coun*rp)
      call a%Memor%deallocObj(0,'tesgs%a','lmn_turnof',tesgs_coun*rp)
      call a%Memor%deallocObj(0,'prsgs%a','lmn_turnof',prsgs_coun*rp)
      call a%Memor%dealloc(nelem,a%vesgs,'vesgs','lmn_turnof')
      call a%Memor%dealloc(nelem,a%tesgs,'tesgs','lmn_turnof')
      call a%Memor%dealloc(nelem,a%prsgs,'prsgs','lmn_turnof')
   end if
   
   !Residual Projection
   if (a%kfl_repro >=1) then
      iaux = max(ndime+2,a%kfl_repro*(ndime+2))
      call a%Memor%dealloc(iaux,npoin,a%repro,'repro','lmn_turnof')
   endif
   
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      call a%Memor%dealloc(ndime+1,ndime,npoin,a%grprj,'grprj','lmn_turnof')
   end if

   !Postprocess the residual at the gauss points
   if (a%npp_stepi(11) /= 0) then
      resid_coun = 0
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%residualU(ielem)%a)
         deallocate(a%residualP(ielem)%a)
         deallocate(a%residualT(ielem)%a)
         resid_coun = resid_coun + (ndime+1)*pgaus
      end do
      call a%Memor%deallocObj(0,'residual%a','lmn_turnof',resid_coun*rp)
      call a%Memor%dealloc(nelem,a%residualU,'residual','lmn_turnof')
      call a%Memor%dealloc(nelem,a%residualP,'residual','lmn_turnof')
      call a%Memor%dealloc(nelem,a%residualT,'residual','lmn_turnof')
   endif
   
   if(a%kfl_outfm==1)then   
      call a%Memor%dealloc(ndime,nbody,a%force,'force','lmn_turnof')
      call a%Memor%dealloc(3,nbody,a%momen,'momen','lmn_turnof')
      call a%Memor%dealloc(nbody,a%heatf,'heatf','lmn_turnof')
      do ibody=1,nbody
         close(a%lun_force(ibody))
      end do
      call a%Memor%dealloc(nbody,a%lun_force,'lun_force','lmn_turnof')
   end if

   if (a%MPIrank == a%MPIroot) then
      call iofile(two,a%lun_trap,a%fil_trap,'LMN TRACKING OF POINTS')
   endif

end subroutine lmn_turnof
