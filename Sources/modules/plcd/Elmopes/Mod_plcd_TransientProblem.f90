module Mod_plcd_TransientProblem
   use typre  
   use Mod_plcd_BaseElmope
   use Mod_PointerSetter
   implicit none
   private
   public SetPointersTransientProblem, SetPointersTransientProblemEndite 
   
   type, extends(PointerSetter) :: SPTransientProblem
contains
      procedure :: SpecificSet => SpecificSetPointersTransientProblem
   end type 
   type(SPTransientProblem) :: SetPointersTransientProblem
   
   type, extends(PointerSetter) :: SPTransientProblemEndite
contains
      procedure :: SpecificSet => SpecificSetPointersTransientProblemEndite
   end type 
   type(SPTransientProblemEndite) :: SetPointersTransientProblemEndite
   
   real(rp), allocatable :: elaccel(:,:), GaussElTransientForces(:,:)
   real(rp), allocatable :: elmass(:,:,:,:)
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SpecificSetPointersTransientProblem(d)
      implicit none
      class(SPTransientProblem) :: d
            !-------------------------------------------------------
            !TransientProblem
            if (a%kfl_TransientProblem == 1) then
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyMassMatrix)
            endif
   end subroutine
   
   subroutine ComputeAndAssemblyMassMatrix
      implicit none
      integer(ip) :: idime, jnode, inode
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      real(rp), pointer :: TimeStep => NULL()
      real(rp) :: aux
      
      TimeStep => a%css%TimeStep
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      aux = 1.0_rp/(a%Beta*TimeStep*TimeStep)
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime)
         elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + aux*density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol
      end forall

   end subroutine ComputeAndAssemblyMassMatrix
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !ENDITE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine SpecificSetPointersTransientProblemEndite(d)
      implicit none
      class(SPTransientProblemEndite) :: d
            !-------------------------------------------------------
            !Temporal Derivative
            if (a%kfl_TransientProblem == 1) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocateElementalAcceleration)
               call ConcatenateProcedures(ProcHook%Finalizations,DeAllocateElementalAcceleration)
               call ConcatenateProcedures(ProcHook%PreGauss,GathersElementalAcceleration)
               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeAndAssemblyTransientForces)
            endif
   end subroutine   
   
   subroutine AllocateElementalAcceleration
      call a%Memor%alloc(e%ndime,e%mnode,elaccel,'elaccel','plcd_TransientProblem')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%pnode,elmass,'elmass','plcd_TransientProblem')
      call a%Memor%alloc(a%ndofn,e%mnode,GaussElTransientForces,'GaussElTransientForces','plcd_RotatingFrame')

   end subroutine AllocateElementalAcceleration
   
   subroutine DeAllocateElementalAcceleration
      call a%Memor%dealloc(e%ndime,e%mnode,elaccel,'elaccel','plcd_TransientProblem')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmass,'elmass','plcd_TransientProblem')
      call a%Memor%dealloc(a%ndofn,e%mnode,GaussElTransientForces,'GaussElTransientForces','plcd_RotatingFrame')
   
   end subroutine DeAllocateElementalAcceleration
   
   subroutine GathersElementalAcceleration
      call e%gather(e%ndime,elaccel(:,:),a%Acceleration(:,:,1))
   
   end subroutine GathersElementalAcceleration
   
   subroutine ComputeAndAssemblyTransientForces
      implicit none
      integer(ip) :: idime, jnode, inode
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      
      GaussElTransientForces = 0.0_rp
      
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime)
         elmass(idime,inode,idime,jnode) = density*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol
      end forall
      
      call matvecmult(e%ndime,e%pnode,elmass,elaccel,GaussElTransientForces)
     
      GaussElInternalForces = GaussElInternalForces + GaussElTransientForces


   end subroutine ComputeAndAssemblyTransientForces

end module
