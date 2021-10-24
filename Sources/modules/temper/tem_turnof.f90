subroutine tem_turnof(a)
   use typre
   use Mod_Temperature
   use Mod_r1pElementAllocation
   use def_parame
   implicit none
   class(TemperatureProblem) :: a
   
   integer(ip) :: ndime,npoin,nelem,tesgs_coun,ncsgs,pnode,pgaus,ielem,nboun,ncomp,ibody,nbody, sig_coun

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNbody(nbody)
   ncomp = a%ncomp
   
   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outph,'(//,a,/)') '     * * * END OF TEMPERATURE RUN * * *'
      write(a%lun_solve,'(//,a,/)') '     * * * END OF TEMPERATURE RUN * * *'
   endif
   
   call a%Memor%dealloc(npoin,a%ncomp,a%tempe,'tempe','tem_turnof')
   
   if (a%kfl_sourc == 2) then
      call a%Memor%dealloc(npoin,a%PointwiseSource,'PointwiseSource','tem_reabcs')
   endif      
   
   if(a%kfl_trasg/=0) then
      tesgs_coun = 0
      ncsgs=2
      
      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         deallocate(a%tesgs(ielem)%a)
         tesgs_coun = tesgs_coun + ncsgs*pgaus
      end do
      call a%Memor%deallocObj(0,'tesgs%a','tem_memall',tesgs_coun*rp)
      call a%Memor%dealloc(nelem,a%tesgs,'tesgs','tem_memall')
   end if

   if(allocated(a%repro)) then
      call a%Memor%dealloc(npoin,a%repro,'repro','tem_memall')
   endif
  
   !Gradient Orthogonal projection
   if(a%kfl_adapsgs == 1) then
      call a%Memor%alloc(ndime,npoin,a%grprj,'grprj','tem_turnof')
   end if
   
   if (allocated(a%gradient)) then
      call a%Memor%dealloc(ndime,npoin,a%gradient,'gradient','tem_memall')
   endif
   
   if (allocated(a%GradientGaussPoints)) then
      call DeAllocateR2pElement(a%Mesh,ndime,a%GradientGaussPoints,a%Memor,'gradient')
   endif
   
   if (a%npp_stepi(8) /= 0) then
      call DeAllocateR1pElement(a%Mesh,a%ShockCapturingViscosity,a%Memor)
   endif
   
   if(a%kfl_dispa /= 0) then
      call a%Memor%dealloc(npoin,a%dissipation,'dissipation','tem_memall')
   endif
   
   if(a%kfl_outfm==1)then   
      call a%Memor%dealloc(nbody,a%heatf,'heatf','tem_turnof')
      do ibody=1,nbody
         close(a%lun_force(ibody))
      end do
      call a%Memor%dealloc(nbody,a%lun_force,'lun_force','tem_turnof')
   end if

 if (a%kfl_CouplingThreeField==1) then
      call a%Memor%dealloc(ndime,ndime,npoin,a%SmoothedVelocityGradient,'SmoothedVelocityGradient','tem_turnof')  
      if (a%npp_stepi(9) /= 0) then  
         sig_coun = 0
         do ielem=1,nelem
            call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)        
            deallocate(a%sigmatermarray(ielem)%a)
            sig_coun = sig_coun + pgaus        
         end do
         call a%Memor%deallocObj(0,'sigmatermarray%a','tem_turnof',sig_coun*rp)
         call a%Memor%dealloc(nelem,a%sigmatermarray,'sigmatermarray','tem_turnof')
         call a%Memor%dealloc(npoin,a%ViscousDissipation,'ViscousDissipation','tem_turnof')
      end if
  end if
   
   !Materials
   deallocate(a%Materials)
   call a%Memor%deallocObj(0,'Materials','tempe_memall',1*a%NumberOfMaterials)

   if (a%NumberOfMaterials > 1) then
      call a%Memor%dealloc(nelem,a%ElementMaterials,'ElementMaterials','tem_memall')
   endif
   
   
end subroutine
