module Mod_plcd_ExternalForces
   use typre
   use Mod_plcd_BaseElmope
   use Mod_plcdExacso
   implicit none
   private
   public SetPointersExternalForces, GetExternalForces
   
   type, extends(PointerSetter) :: SPExternalForces
contains
      procedure :: SpecificSet => SpecificSetExternalForces
   end type
   
   type(SPExternalForces) :: SetPointersExternalForces
   
   integer(ip), allocatable :: kfl_IsSet
   
   type(plcdExacso) :: Exacso
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetExternalForces(d)
      class(SPExternalForces) :: d
         
      !Element size required by EMDS
      if (a%kfl_exacs > 0) then
         call ConcatenateProcedures(ProcHook%InGaussElmats,ExacsoExternalForces)
      endif
      
      !Gravity Force
      if (a%kfl_GravityForce == 1) then
         call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeGravityForces)
      endif
      
   end subroutine   
   
   
   !For actual computations
   !------------------------------------------------------
   !Hanging nodes
   subroutine ExacsoExternalForces
      implicit none
      
      integer(ip) ::  inode
      real(rp) :: elext(4), gpcod(3)
      
      
      elext = 0.0_rp   
      call e%interpg(e%ndime,e%elcod,gpcod)
      call Exacso%ComputeSolution(e%ndime,gpcod,a)
      call Exacso%GetForce(e%ndime,elext,a)
      do inode = 1,e%pnode
         GaussElExternalForces(1:a%ndofn,inode) = GaussElExternalForces(1:a%ndofn,inode) + elext(1:a%ndofn)*e%shape(inode,e%igaus)*dvol
      enddo
      
      
   end subroutine
   
   subroutine GetExternalForces(e2,a2,elext2)
      implicit none
      class(FiniteElement) :: e2
      class(PLCDProblem) :: a2
      real(rp) :: elext2(*)
      
      
      real(rp) :: gpcod(3)
      
      elext2(1:e2%ndime+1) = 0.0_rp   
      if (a2%kfl_exacs > 0) then
          call e2%interpg(e2%ndime,e2%elcod,gpcod)
          call Exacso%ComputeSolution(e2%ndime,gpcod,a2)
          call Exacso%GetForce(e2%ndime,elext2(1:e2%ndime+1),a2)
      endif
   end subroutine
   
   subroutine ComputeGravityForces
      implicit none
      integer(ip) :: inode
      class(PLCDMaterial), pointer :: Material
      real(rp), pointer :: density => NULL()
      
      call ElementMatData%GetMaterialPointer(Material)
      density => Material%density
      
      do inode = 1,e%pnode
         GaussElExternalForces(1:a%ndofn,inode) = GaussElExternalForces(1:a%ndofn,inode) + e%shape(inode,e%igaus)*density*a%gravity(1:a%ndofn)*dvol
      enddo
      
   end subroutine ComputeGravityForces
   
end module
