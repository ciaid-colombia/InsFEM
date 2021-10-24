subroutine plcd_ComputeInitialAcceleration(b)
   use Mod_plcd_BaseElmope
   use Mod_PLCD
   use Mod_plcd_BMatrixFactory
   use Mod_plcd_Stages
   use Mod_plcd_elmdir
   use Mod_php_AssemblyVectorToSystem
   use Mod_Debugging
   
   implicit none
   class(PLCDProblem), target :: b
   integer(ip) :: idofn, inode,jnode,ipoin,idime,jdime, ielem2
   class(PLCDMaterial), pointer :: Material
   real(rp), pointer :: density => NULL()
   real(rp), pointer :: TimeStep => NULL()
   real(rp) :: aux
   
   interface
      subroutine plcd_ComputeForcesVectors(a)
         use typre
         use Mod_PLCD
         implicit none
         class(PLCDProblem) :: a
      end subroutine
   end interface
   

   !deb_PostprocessMatrix = 1

   a => b

   call plcd_ComputeForcesVectors(a)
   call php_AssemblyVectorToSystem(a,a%ResidualForcesVector)
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_ComputeInitialAcceleration')
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','plcd_ComputeInitialAcceleration')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','plcd_ComputeInitialAcceleration')

   
   call a%Mesh%GetNelem(nelem)
   do ielem = 1,nelem
      ielem2 = ielem

      !Load Element
      call a%Mesh%ElementLoad(ielem,e)

      elmat = 0.0_rp
      elrhs = 0.0_rp

      !Compute linear derivatives
      call e%elmdel

      ElementMatData => a%ElementMaterialsData(ielem)%p

      !Gausspoint loop
      GaussPointLoop : do igaus = 1,e%pgaus
         e%igaus = igaus

         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         
         TimeStep => a%css%TimeStep
         call ElementMatData%GetMaterialPointer(Material)
         density => Material%density
      
         aux = 1.0_rp/(a%Beta*TimeStep*TimeStep)
      
         forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime)
            elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + aux*density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol
         end forall
         
      enddo GaussPointLoop

      !Dirichlet Boundary Conditions
      call plcd_elmdir(a,e,elmat,elrhs)

      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   enddo
   
   !We solve for initial accelerations
   call a%LinearSystem%Solve(a%Acceleration(:,:,1))
   
   call a%Mesh%ArrayCommunicator%GhostCommunicate(size(a%Acceleration,1),a%Acceleration(:,:,1))

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','plcd_ComputeInitialAcceleration')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','plcd_ComputeInitialAcceleration')
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_ComputeInitialAcceleration')

end subroutine 