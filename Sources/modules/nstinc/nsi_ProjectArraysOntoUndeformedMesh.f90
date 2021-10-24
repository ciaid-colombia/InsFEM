subroutine nsi_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
   use typre
   use Mod_MeshInterpolator
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   type(Interpolator) :: Interp
   integer(ip) :: itask
   
   integer(ip) :: ndime,icomp
   
   if (itask == 1) then
   
      call a%Mesh%GetNdime(ndime)
      do icomp = 1,size(a%veloc,3)
         call Interp%Interpolate(ndime,a%veloc(:,:,icomp),a%veloc(:,:,icomp))
         call Interp%Interpolate(1_ip,a%press(:,icomp),a%press(:,icomp))
      enddo

      if (a%kfl_StabilizeFreeSurface == 1) then
      
         call Interp%Interpolate(ndime*ndime,a%StabilizeFreeSurface_VelocityGradient,a%StabilizeFreeSurface_VelocityGradient)
         call Interp%Interpolate(ndime,a%StabilizeFreeSurface_PressureGradientAndForce,a%StabilizeFreeSurface_PressureGradientAndForce)
         
      endif
   
   
      !Transient subgrid scales
      if(a%kfl_trasg/=0) then
         call runend('nsi fmale not ready for trasgs')
      end if
      
      !Residual Projection
      if (a%kfl_repro >=1) then
         call Interp%Interpolate(size(a%repro,1),a%repro,a%repro)
      endif
      
      !Tau Smoothing
      if (a%kfl_Tausm >= 1) then
         do icomp = 1,2
            call Interp%Interpolate(1_ip,a%Tausmo(:,icomp),a%Tausmo(:,icomp))
         enddo
      endif
      
      !Postprocess the residual at the gauss points
      if (a%npp_stepi(18) /= 0) then
               
      endif
      
      !Gauss point dissipation
      if (a%npp_stepi(20) /= 0) then
         call runend('nsi gauss point dissipation not ready for fmale')
      endif
      
      
   elseif (itask == 2) then   
      !Nothing to do in the second pass
      
      
      
      
   endif
   
   


end subroutine



subroutine nsi_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
   use typre
   use Mod_Advector
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem), target :: a
   type(Advector) :: Advect
   integer(ip) :: itask
   
   integer(ip) :: ndime,icomp
   real(rp), pointer :: auxpointer(:,:) => NULL()
   
   if (itask == 1) then
!       call a%FilePostpr%postpr(a%veloc(:,:,1),'PreProjVeloc1',a%istep,a%ctime,a%Mesh)
!       call a%FilePostpr%postpr(a%veloc(:,:,3),'PreProjVeloc3',a%istep,a%ctime,a%Mesh)

   
      call a%Mesh%GetNdime(ndime)
      do icomp = 1,size(a%veloc,3)
         call Advect%Advect(ndime,a%veloc(:,:,icomp),a%veloc(:,:,icomp))
         
         auxpointer(1:1,1:size(a%press,1)) => a%press(:,icomp)
         call Advect%Advect(1_ip,auxpointer,auxpointer)
      enddo
!       call a%FilePostpr%postpr(a%veloc(:,:,1),'ProjectedVeloc1',a%istep,a%ctime,a%Mesh)
!       call a%FilePostpr%postpr(a%veloc(:,:,3),'ProjectedVeloc3',a%istep,a%ctime,a%Mesh)

      if (a%kfl_StabilizeFreeSurface == 1) then
!          call a%FilePostpr%postpr(a%StabilizeFreeSurface_PressureGradientAndForce(:,:),'grapress_PREinterp',a%istep,a%ctime,a%Mesh)
         
         auxpointer(1:ndime*ndime,1:size(a%StabilizeFreeSurface_VelocityGradient,3)) => a%StabilizeFreeSurface_VelocityGradient
         call Advect%Advect(ndime*ndime,auxpointer,auxpointer)
         call Advect%Advect(ndime,a%StabilizeFreeSurface_PressureGradientAndForce,a%StabilizeFreeSurface_PressureGradientAndForce)
         
!          call a%FilePostpr%postpr(a%StabilizeFreeSurface_PressureGradientAndForce(:,:),'grapress_interp',a%istep,a%ctime,a%Mesh)

      endif
   
   
      !Transient subgrid scales
      if(a%kfl_trasg/=0) then
         call runend('nsi fmale not ready for trasgs')
      end if
      
      !Residual Projection
      if (a%kfl_repro >=1) then
         call Advect%Advect(size(a%repro,1),a%repro,a%repro)
      endif
      
      !Tau Smoothing
      if (a%kfl_Tausm >= 1) then
         do icomp = 1,2
            auxpointer(1:1,1:size(a%Tausmo,1)) => a%Tausmo(:,icomp)
            call Advect%Advect(1_ip,auxpointer,auxpointer)
         enddo
      endif
      
      !Postprocess the residual at the gauss points
      if (a%npp_stepi(18) /= 0) then
               
      endif

      
      !Gauss point dissipation
      if (a%npp_stepi(20) /= 0) then
         call runend('nsi gauss point dissipation not ready for fmale')
      endif
      
      
   elseif (itask == 2) then   
      !Nothing to do in the second pass
      
   endif
   
end subroutine
