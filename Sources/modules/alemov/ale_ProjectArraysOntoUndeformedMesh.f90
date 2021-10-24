subroutine ale_ProjectArraysOntoUndeformedMesh(a,Interp,itask)
   use typre
   use Mod_Alemov
   use Mod_MeshInterpolator
   implicit none
   class(AlemovProblem) :: a
   type(Interpolator) :: Interp
   integer(ip) :: itask
   
   integer(ip) :: ndime,icomp,ipoin,npoin,idime
   real(rp), allocatable :: auxdisplacement(:,:,:), auxvelocity(:,:,:)
   
   if (itask == 1) then
   
      call a%Memor%alloc(size(a%Displacement,1),size(a%Displacement,2), size(a%Displacement,3),auxdisplacement,'auxdisplacement','ale_ProjectArraysOntoUndeformedMesh')
      call a%Memor%alloc(size(a%velocity,1),size(a%velocity,2), size(a%velocity,3),auxvelocity,'auxvelocity','ale_ProjectArraysOntoUndeformedMesh')
   
      auxdisplacement = a%Displacement
      auxvelocity = a%Velocity
   
      do icomp = 1,size(a%Displacement,3)
         call Interp%Interpolate(size(a%Displacement,1),a%Displacement(:,:,icomp),a%Displacement(:,:,icomp))
      enddo
      
      do icomp = 1,size(a%Velocity,3)
         call Interp%Interpolate(size(a%Velocity,1),a%Velocity(:,:,icomp),a%Velocity(:,:,icomp))
      enddo
   
      !If the boundary condition is prescribed, then do not let the projection override the values
      !This avoids problems in the case of folded elements (FMALE)
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         do idime = 1,size(a%kfl_fixno,1)
            if (a%kfl_fixno(idime,ipoin) == 1) then
               a%Displacement(idime,ipoin,:) = auxdisplacement(idime,ipoin,:)
               a%Velocity(idime,ipoin,:) = auxvelocity(idime,ipoin,:)
            endif
         enddo
      enddo
   
      call a%Memor%dealloc(size(a%Displacement,1),size(a%Displacement,2), size(a%Displacement,3),auxdisplacement,'auxdisplacement','ale_ProjectArraysOntoUndeformedMesh')
      call a%Memor%dealloc(size(a%velocity,1),size(a%velocity,2), size(a%velocity,3),auxvelocity,'auxvelocity','ale_ProjectArraysOntoUndeformedMesh')

   elseif (itask == 2) then
   
      !The mesh velocity does not change, but the displacement is zero
      !The previous displacements need to be substracted the current displacement
      do icomp = 2,size(a%Displacement,3)
         a%Displacement(:,:,icomp) = a%Displacement(:,:,icomp)- a%Displacement(:,:,1_ip)
      enddo
      a%Displacement(:,:,1) = 0.0_rp

      !we can now recompute vmass and Extnor
      !Recompute vmass and extnorlpoty, since we go back to undeformed mesh!
      !Dealloc and recompute ExtnorLpoty and Vmass
      call a%Mesh%DeallocExnorLpoty
      call a%Mesh%DeallocVmass
      
      !Recompute
      call a%Mesh%ComputeVmass
      call a%Mesh%ExtnorLpoty
   
   endif
   
end subroutine




subroutine ale_AdvectArraysOntoUndeformedMesh(a,Advect,itask)
   use typre
   use Mod_Mesh
   use Mod_Alemov
   use Mod_Advector
   implicit none
   class(AlemovProblem) :: a
   type(Advector) :: Advect
   integer(ip) :: itask
   
   integer(ip) :: ndime,icomp,ipoin,npoin,idime
   real(rp), allocatable :: auxdisplacement(:,:,:), auxvelocity(:,:,:)
   
   if (itask == 1) then
   
       call a%FilePostpr%postpr(a%Displacement(:,:,1),'aledisppreproj',a%istep,a%ctime,a%Mesh)
       call a%FilePostpr%postpr(a%Velocity(:,:,1),'aleveopreproj',a%istep,a%ctime,a%Mesh)
   
   
      call a%Memor%alloc(size(a%Displacement,1),size(a%Displacement,2), size(a%Displacement,3),auxdisplacement,'auxdisplacement','ale_ProjectArraysOntoUndeformedMesh')
      call a%Memor%alloc(size(a%velocity,1),size(a%velocity,2), size(a%velocity,3),auxvelocity,'auxvelocity','ale_ProjectArraysOntoUndeformedMesh')
   
      auxdisplacement = a%Displacement
      auxvelocity = a%Velocity
   
      do icomp = 1,size(a%Displacement,3)
         call Advect%Advect(size(a%Displacement,1),a%Displacement(:,:,icomp),a%Displacement(:,:,icomp))
      enddo
       
      do icomp = 1,size(a%Velocity,3)
         call Advect%Advect(size(a%Velocity,1),a%Velocity(:,:,icomp),a%Velocity(:,:,icomp))
      enddo
      !a%Velocity = auxvelocity
      
   
      !If the boundary condition is prescribed, then do not let the projection override the values
      !This avoids problems in the case of folded elements (FMALE)
      call a%Mesh%GetNpoin(npoin)
      do ipoin = 1,npoin
         do idime = 1,size(a%kfl_fixno,1)
            if (a%kfl_fixno(idime,ipoin) == 1) then
               a%Displacement(idime,ipoin,:) = auxdisplacement(idime,ipoin,:)
               a%Velocity(idime,ipoin,:) = auxvelocity(idime,ipoin,:)
            endif
         enddo
      enddo
   
      call a%Memor%dealloc(size(a%Displacement,1),size(a%Displacement,2), size(a%Displacement,3),auxdisplacement,'auxdisplacement','ale_ProjectArraysOntoUndeformedMesh')
      call a%Memor%dealloc(size(a%velocity,1),size(a%velocity,2), size(a%velocity,3),auxvelocity,'auxvelocity','ale_ProjectArraysOntoUndeformedMesh')
   
   
       call a%FilePostpr%postpr(a%Displacement(:,:,1),'aledisppostproj',a%istep,a%ctime,a%Mesh)
       call a%FilePostpr%postpr(a%Velocity(:,:,1),'aleveopostproj',a%istep,a%ctime,a%Mesh)
    
    
   elseif (itask == 2) then
   
!       !The mesh velocity does not change, but the displacement is zero
!       !The previous displacements need to be substracted the current displacement
!        do icomp = 2,size(a%Displacement,3)
!           a%Displacement(:,:,icomp) = a%Displacement(:,:,icomp)- a%Displacement(:,:,1_ip)
!        enddo
!        a%Displacement(:,:,1) = 0.0_rp
      
!       call a%FilePostpr%postpr(a%Displacement(:,:,1),'aledisppostproj2',a%istep,a%ctime,a%Mesh)
!       call a%FilePostpr%postpr(a%Velocity(:,:,1),'aleveopostproj2',a%istep,a%ctime,a%Mesh)
!    
!       
       !we can now recompute vmass and Extnor
       !Recompute vmass and extnorlpoty, since we go back to undeformed mesh!
       !Dealloc and recompute ExtnorLpoty and Vmass
       call a%Mesh%DeallocExnorLpoty
       call a%Mesh%DeallocVmass
       
       !Recompute
       call a%Mesh%ComputeVmass
       call a%Mesh%ExtnorLpoty
       
 !    elseif (itask == 3) then
 !       !Update boundary conditions for ALE (necessary if Fixed-Mesh ALE)
 !       call a%Updbcs
     
    endif
    
   
end subroutine
