subroutine plcd_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* lmachn/plcd_endite
   ! NAME
   !    plcd_endite
   ! DESCRIPTION
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_PLCD
   use Mod_int2str
   use Mod_SmoothedFieldGradient
   use Mod_Debugging
   implicit none
   class(PLCDProblem) :: a
   integer(ip) :: itask
   real(rp), pointer :: TimeStep => NULL()

   interface
      subroutine plcd_ComputeForcesProjection(a,ForcesProjection)
         use typre
         use Mod_PLCD
         implicit none
         class(PLCDProblem) :: a
         real(rp) :: ForcesProjection(:,:)
      end subroutine
   end interface

   integer(ip) :: ndime,npoin
   real(rp), allocatable :: ForcesProjection(:,:)
   character(50) :: names

   interface
      subroutine plcd_EnditeElmope(b)
         import
         implicit none
         class(PLCDProblem), target :: b
      end subroutine

      subroutine plcd_PostprocessMaterialData(b,iiter_char)
         import
         implicit none
         class(PLCDProblem), target :: b
         character(5) :: iiter_char
      end subroutine

      subroutine plcd_ComputeForcesVectors(a)
         use typre
         use Mod_PLCD
         implicit none
         class(PLCDProblem) :: a
      end subroutine

   end interface

   if (itask == 10) then

      call a%Mesh%GetNdime(ndime)
      a%Displacement(:,:,1) = a%unkno(1:ndime,:)
         
      !UP Formulation
      !Update pressure and compute pressure gradient projection
      if (a%UseUPFormulation) then
         a%Pressure(:,1) = a%unkno(ndime+1,:)

         call ComputeSmoothedFieldGradient(a%Mesh,a%Memor,1_ip,a%Pressure(:,1),a%UPResidualProjection)

         if (a%kfl_exacs > 0) then
            !Forces Projection
            call a%Mesh%GetNdime(ndime)
            call a%Mesh%Getnpoin(npoin)
            call a%Memor%alloc(a%ndofn,npoin,ForcesProjection,'ForcesProjection','plcd_endite')
            call plcd_ComputeForcesProjection(a,ForcesProjection)

            a%UPResidualProjection(1:ndime,:) = a%UPResidualProjection(1:ndime,:) - ForcesProjection(1:ndime,:)

            !call a%FilePostpr%postpr(ForcesProjection,'ForcesProjection',a%istep,a%ctime,a%Mesh)
            !call a%FilePostpr%postpr(a%UPResidualProjection,'ResidualProjection',a%istep,a%ctime,a%Mesh)
            !call a%FilePostpr%postpr(a%PForcesProjection,'PForcesProjection',a%istep,a%ctime,a%Mesh)

            call a%Memor%dealloc(a%ndofn,npoin,ForcesProjection,'ForcesProjection','plcd_endite')
         endif
      endif


      !Compute the smoothed displacement gradient if necessary
      if (a%UseSmoothedDisplacementGradient) then

         call a%Mesh%GetNdime(ndime)
         call ComputeSmoothedFieldGradient(a%Mesh,a%Memor,ndime,a%Displacement(:,:,1),a%SmoothedDisplacementGradient)
         !call a%FilePostpr%postpr(a%SmoothedDisplacementGradient,'SmoothGrad'//int2str(a%itera),a%istep,a%ctime,a%Mesh)
      endif
      
      if (a%kfl_TransientProblem == 1) then
         TimeStep => a%css%TimeStep
      
         a%Acceleration(:,:,1) = (1.0_rp/(a%Beta*TimeStep*TimeStep))*(a%Displacement(:,:,1)-a%Displacement(:,:,3))
         a%Velocity(:,:,1) = a%Velocity(:,:,3) + a%Gamma*TimeStep*a%Acceleration(:,:,1)
      
      endif

      !compute the Forces Vectors
      call plcd_ComputeForcesVectors(a)

      if (a%kfl_PostprocessMatDataAtEachIteration == 1) then
         call plcd_PostprocessMaterialData(a,int2str(a%itera))
      endif

      if (a%kfl_PostprocessDisplacementAtEachIteration == 1) then
         names = 'displacement'//trim(adjustl(int2str(a%itera)))
         call a%FilePostpr%postpr(a%Displacement(:,:,1),names,a%istep,a%ctime,a%Mesh)
      endif


   elseif (itask == 1) then


   endif


end subroutine plcd_endite



