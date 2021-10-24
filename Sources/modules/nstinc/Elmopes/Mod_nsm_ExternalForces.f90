 module Mod_nsm_ExternalForces
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersExternalForces

   type, extends(PointerSetter) :: SPExternalForces
contains
      procedure :: SpecificSet => SpecificSetExternalForces
   end type
   type(SPExternalForces) :: SetPointersExternalForces

   real(rp) :: acceleration_module
   logical :: aux_logic

contains

   subroutine SpecificSetExternalForces(d)
      implicit none
      class(SPExternalForces) :: d
      real(rp), external :: funcre

      ProcPointer_ExternalForces  => nsiForces

      !ExactSol
      if (a%kfl_exacs/=0) then
         ProcPointer_ExternalForces => ExactSolutionForces
      endif

      if (a%kfl_TurbulentBodyForces == 1) then
         call ConcatenateProcedures(ProcPointer_ExternalForces,nsm_TurbulentBodyForces)
      endif

      if (a%kfl_ExitBodyForces == 1) then
         call ConcatenateProcedures(ProcPointer_ExternalForces,nsm_ExitBodyForces)
      endif
      
      if (a%kfl_Plasma == 1) then
         call ConcatenateProcedures(ProcPointer_ExternalForces,nsm_PlasmaForces)
      endif

      if (a%kfl_FORAcceleration == 1) then
         !We compute the acceleration only once
         acceleration_module = funcre(a%FORParam,size(a%FORParam),a%kfl_FORFunty,a%ctime)
         call ConcatenateProcedures(ProcPointer_ExternalForces,AccelerationForces)
      endif
      
      !Coriolis
      if (a%kfl_CoriolisForce == 1) then
         call ConcatenateProcedures(ProcPointer_ExternalForces,CoriolisForces)
      endif
      
      !Logics for deciding if we need postprocess forces
      aux_logic = .false.
      if (a%npp_stepi(25) /= 0) then   
         if (mod(a%istep,a%npp_stepi(25))==0) then
            aux_logic = .true.
         endif
      endif
      
      if (aux_logic) then
         call ConcatenateProcedures(ProcHook_Initializations,PostprocessExternalForcesToZero)
         call ConcatenateProcedures(ProcPointer_ExternalForces,PostprocessExternalForces)
         call ConcatenateProcedures(ProcHook_Finalizations,PostprocessExternalForcesFinalizations)
      endif

   end subroutine

   !---------------------------------------------------------------------------
   !Non Exact Solution Case
   subroutine  nsiForces
      implicit none
      !Compute vector of external forces
      call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext)
   end subroutine

   subroutine  ExactSolutionForces
      implicit none
      real(rp)    :: gpcod(e%ndime)
      !Interpolate
      call e%interpg(e%ndime,e%elcod,gpcod)
      !Compute vector of external forces
      call exacso%nsi_ComputeSolution(e%ndime,gpcod,a%ctime,a)
      call exacso%nsi_GetForce(e%ndime,elext,a)
   end subroutine

   subroutine AccelerationForces
      call nsi_ComputeExternalForces(e,acden,acceleration_module,a%FORAcceleration,elext)
   end subroutine

   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine nsm_TurbulentBodyForces
       implicit none
       real(rp) ::  gpcod(3), force(3)

       call e%interpg(e%ndime,e%elcod,gpcod)
       force(1:e%ndime) = 0.0_rp
       call a%TBF%GetForce(a%ctime,gpcod,force)
       elext(1:e%ndime) = elext(1:e%ndime) + force(1:e%ndime)
   end subroutine

   subroutine nsm_ExitBodyForces
      implicit none
      real(rp) :: gpcod(3), force(3)

      call e%interpg(e%ndime,e%elcod,gpcod)
      force(1:e%ndime) = 0.0_rp
      call a%EBF%GetForce(a%ctime,gpcod,gpvel,force)
      elext(1:e%ndime) = elext(1:e%ndime) + force(1:e%ndime)
   end subroutine
   
   subroutine nsm_PlasmaForces
      implicit none
      real(rp) :: gpcod(3), Pforce(3)

      call e%interpg(e%ndime,e%elcod,gpcod)
      Pforce(1:e%ndime) = 0.0_rp
      call a%PAC%PlasmaForce(gpcod,Pforce)
      elext(1:e%ndime) = elext(1:e%ndime) + Pforce(1:e%ndime) 
   end subroutine
   
   subroutine PostprocessExternalForcesToZero
      implicit none
      
      a%ExternalForcesArray = 0.0_rp
   end subroutine
   
   subroutine PostprocessExternalForces
      implicit none
      real(rp) :: elextArray(e%ndime,e%pnode)
      integer(ip) :: inode
      
      do inode = 1,e%pnode
         elextArray(:,inode) = e%shape(inode,e%igaus)*elext(1:e%ndime)*dvol
      enddo   
      call a%Mesh%AssemblyToArray(e,e%ndime,elextArray,a%ExternalForcesArray) 
   end subroutine
   
   subroutine PostprocessExternalForcesFinalizations
      implicit none
      integer(ip) :: ndime
      
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%Smooth(ndime,a%ExternalForcesArray)
   end subroutine
   
   subroutine CoriolisForces
      implicit none
      real(rp) :: gpcod(3), force(3), w(3),aux(3),auxgpvel(3)
      
      w = a%CoriolisW
      !Centrifugal force
      gpcod = 0.0_rp
      call e%interpg(e%ndime,e%elcod,gpcod)
      force = 0.0_rp
      call vecpro(w,gpcod,aux,3)
      call vecpro(-w,aux,force,3)
      
      !Coriolis force
      !We implement it explicitly as a first approximation
      auxgpvel = 0.0_rp
      auxgpvel(1:e%ndime) = gpvel(1:e%ndime,1)
      call vecpro(w,auxgpvel,aux,3)
      force = force - 2*aux
      
      !Multiply by the density
      force = force*acden
      elext(1:e%ndime) = elext(1:e%ndime) + force(1:e%ndime)
   end subroutine
      
end module
