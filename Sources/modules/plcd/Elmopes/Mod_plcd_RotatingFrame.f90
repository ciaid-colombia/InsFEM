module Mod_plcd_RotatingFrame
   use typre  
   use Mod_plcd_BaseElmope
   use Mod_PointerSetter
   implicit none
   private
   public SetPointersRotatingFrame, SetPointersRotatingFrameEndite 
   
   type, extends(PointerSetter) :: SPRotatingFrame
contains
      procedure :: SpecificSet => SpecificSetPointersRotatingFrame
   end type 
   type(SPRotatingFrame) :: SetPointersRotatingFrame
   
   type, extends(PointerSetter) :: SPRotatingFrameEndite
contains
      procedure :: SpecificSet => SpecificSetPointersRotatingFrameEndite
   end type 
   type(SPRotatingFrameEndite) :: SetPointersRotatingFrameEndite
   
   real(rp) :: vorticity(3,3)
   
   real(rp), allocatable :: elcoriolis(:,:,:,:),elmass(:,:,:,:),elcent(:,:,:,:),elcentrifugal(:,:,:,:)
   real(rp), allocatable :: GaussElCoriolisForces(:,:),GaussElCentrifugalForces(:,:)
   
   real(rp), allocatable :: elveloc(:,:)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SpecificSetPointersRotatingFrame(d)
      implicit none
      class(SPRotatingFrame) :: d
            !-------------------------------------------------------
            !Rotating Frame of Reference
            if (a%kfl_RotatingFrame == 1) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyCentrifugalTermtoJacobian)
               if (a%kfl_TransientProblem /= 0) then
                  call ConcatenateProcedures(ProcHook%Initializations,ComputeVorticity)
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyCoriolisTermtoJacobian)
               endif
            endif
   end subroutine
   
   subroutine ComputeAndAssemblyCoriolisTermtoJacobian
      implicit none
      integer(ip) :: idime, jnode, inode, jdime, ldime
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      real(rp), pointer :: TimeStep => NULL()
      real(rp) :: aux
      
      TimeStep => a%css%TimeStep
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      aux = a%Gamma/(a%Beta*TimeStep)
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime, jdime=1:e%ndime)
         elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) - aux*2.0_rp*density*e%shape(inode,e%igaus)*vorticity(idime,jdime)*e%shape(jnode,e%igaus)*dvol
      end forall

   end subroutine ComputeAndAssemblyCoriolisTermtoJacobian
   
   subroutine ComputeAndAssemblyCentrifugalTermtoJacobian
      implicit none
      integer(ip) :: idime, jnode, inode, jdime, ldime
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      real(rp), pointer :: TimeStep => NULL()
      real(rp) :: aux
      
      TimeStep => a%css%TimeStep
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime)
         elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + a%angvelocitynorm2*density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol 
      end forall
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime,jdime=1:e%ndime)
         elmat(idime,inode,jdime,jnode) = elmat(idime,inode,jdime,jnode) - a%angvelocity(idime)*a%angvelocity(jdime)*density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol 
      endforall

   end subroutine ComputeAndAssemblyCentrifugalTermtoJacobian
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ENDITE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine SpecificSetPointersRotatingFrameEndite(d)
      implicit none
      class(SPRotatingFrameEndite) :: d
            !-------------------------------------------------------
            !Rotating Frame of Reference
            if (a%kfl_RotatingFrame /= 0) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocateCentrifugalTerms)
               call ConcatenateProcedures(ProcHook%Finalizations,DeAllocateCentrifugalTerms)
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyCentrifugalForces)
               
               if (a%kfl_TransientProblem /= 0) then
                  call ConcatenateProcedures(ProcHook%Initializations,AllocateCoriolisTerms)
                  call ConcatenateProcedures(ProcHook%Finalizations,DeAllocateCoriolisTerms)
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyCoriolisForces)
                  call ConcatenateProcedures(ProcHook%PreGauss,GathersElementalVelocity)
               endif
         
            endif
   end subroutine   
   
   subroutine AllocateCentrifugalTerms
      call a%Memor%alloc(a%ndofn,e%mnode,GaussElCentrifugalForces,'GaussElCentrifugalForces','plcd_RotatingFrame')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmass,'elmass','plcd_RotatingFrame')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elcentrifugal,'elcentrifugal','plcd_RotatingFrame')

   end subroutine AllocateCentrifugalTerms
   
   subroutine AllocateCoriolisTerms
      call a%Memor%alloc(a%ndofn,e%mnode,elveloc,'elveloc','plcd_RotatingFrame')
      call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elcoriolis,'elcoriolis','plcd_RotatingFrame')
      call a%Memor%alloc(a%ndofn,e%mnode,GaussElCoriolisForces,'GaussElCoriolisForces','plcd_RotatingFrame')
      
      call ComputeVorticity

   end subroutine AllocateCoriolisTerms
   
   subroutine DeAllocateCentrifugalTerms
      call a%Memor%dealloc(a%ndofn,e%mnode,GaussElCentrifugalForces,'GaussElCentrifugalForces','plcd_RotatingFrame')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmass,'elmass','plcd_RotatingFrame')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elcentrifugal,'elcentrifugal','plcd_RotatingFrame')
   
   end subroutine DeAllocateCentrifugalTerms
   
   subroutine DeAllocateCoriolisTerms
      call a%Memor%dealloc(e%ndime,e%mnode,elveloc,'elveloc','plcd_RotatingFrame')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elcoriolis,'elcoriolis','plcd_RotatingFrame')
      call a%Memor%dealloc(a%ndofn,e%mnode,GaussElCoriolisForces,'GaussElCoriolisForces','plcd_RotatingFrame')
   
   end subroutine DeAllocateCoriolisTerms
   
   subroutine GathersElementalVelocity
      call e%gather(e%ndime,elveloc(:,:),a%Velocity(:,:,1))
   
   end subroutine GathersElementalVelocity
   
   subroutine ComputeAndAssemblyCoriolisForces
      implicit none
      integer(ip) :: idime, jnode, inode, jdime, ldime
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      
      GaussElCoriolisForces = 0.0_rp
      
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime, jdime=1:e%ndime)
         elcoriolis(idime,inode,jdime,jnode) = 2.0_rp*density*e%shape(inode,e%igaus)*vorticity(idime,jdime)*e%shape(jnode,e%igaus)*dvol
      end forall
      
      call matvecmult(e%ndime,e%pnode,elcoriolis,elveloc,GaussElCoriolisForces)
      
      GaussElInternalForces = GaussElInternalForces - GaussElCoriolisForces
      
   end subroutine ComputeAndAssemblyCoriolisForces
   
   subroutine ComputeAndAssemblyCentrifugalForces
      implicit none
      integer(ip) :: idime, jnode, inode, jdime
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      real(rp) :: aux(e%ndime,e%pnode)
      
      GaussElCentrifugalForces = 0.0_rp
      
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime, jdime=1:e%ndime)
         elcentrifugal(idime,inode,jdime,jnode) = a%angvelocity(idime)*a%angvelocity(jdime)*density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol
      end forall
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime)
         elmass(idime,inode,idime,jnode) = density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol
      end forall
      
      elcentrifugal = elcentrifugal-(a%angvelocitynorm2)*elmass 
      if (a%kfl_LargeStrains == 2) then
         aux(:,:) = e%elcod(:,:) + eldisp(:,:,1)
      else
         aux(:,:) = e%elcod(:,:) 
      endif
      
      call matvecmult(e%ndime,e%pnode,elcentrifugal,aux,GaussElCentrifugalForces)
      GaussElInternalForces = GaussElInternalForces - GaussElCentrifugalForces

   end subroutine ComputeAndAssemblyCentrifugalForces
   
   subroutine ComputeVorticity
   implicit none
   
   vorticity = 0.0_rp
   
   vorticity(1,2) = -a%angvelocity(3)
   vorticity(1,3) = a%angvelocity(2)
   vorticity(2,1) = -vorticity(1,2)
   vorticity(2,3) = -a%angvelocity(1)
   vorticity(3,1) = -vorticity(1,3)
   vorticity(3,2) = -vorticity(2,3)
   
   end subroutine ComputeVorticity

end module
