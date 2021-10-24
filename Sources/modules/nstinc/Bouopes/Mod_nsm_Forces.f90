module Mod_nsm_Forces
   use Mod_PointerSetter
   use Mod_nsm_HangingNodes
   use Mod_nsm_BaseBouope
   implicit none
   private
   public SetPointersForces

   type, extends(PointerSetter) :: SPForces
contains
      procedure :: SpecificSet => SpecificSetForces
   end type
   type(SPForces) :: SetPointersForces

   real(rp), allocatable, dimension(:,:) :: traction
   real(rp), allocatable, dimension(:)   :: elvmass

contains

   subroutine SpecificSetForces(d)
      implicit none
      class(SPForces) :: d

      call ConcatenateProcedures(ProcHook_Initializations,AllocForces)
      call ConcatenateProcedures(ProcHook_Initializations,ToZero)
      call ConcatenateProcedures(ProcHook_ElmatsToZero,BouToZero)
      call ConcatenateProcedures(ProcHook_InGauss,CalculateForce)
      call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,TractionGauss)
      call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyTraction)
      call ConcatenateProcedures(ProcHook_Finalizations,ExtraOp)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocForces)

   end subroutine

   subroutine AllocForces
      implicit none
      integer(ip) :: npoin

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)     
      call a%Memor%alloc(e%mnode,elvmass,'elvmass','nsi_forces')
      call a%Memor%alloc(npoin,vmassBoundary,'vmassBoundary','nsi_forces')
      call a%Memor%alloc(ndime,e%mnode,traction,'traction','nsi_forces')   
   end subroutine

   subroutine DeallocForces
      implicit none
      integer(ip) :: npoin

      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoin(npoin)     
      call a%Memor%dealloc(npoin,vmassBoundary,'vmassBoundary','nsi_forces')
      call a%Memor%dealloc(e%mnode,elvmass,'elvmass','nsi_forces')
      call a%Memor%dealloc(ndime,e%mnode,traction,'traction','nsi_forces')   
   end subroutine

   subroutine ToZero
      implicit none

      a%btraction = 0.0_rp
      vmassBoundary = 0.0_rp
   end subroutine

   subroutine BouToZero
      implicit none

      elvmass = 0.0_rp
      traction = 0.0_rp
   end subroutine

   subroutine CalculateForce
      implicit none
      integer(ip) :: jdime

      !  t=n*S  S=-p*I+2mu*gradS(u)                  
      tract = 0.0_rp
      do idime=1,ndime
          tract(idime)= tract(idime) - e%baloc(idime,ndime)*gpbpr(1)
          do jdime=1,ndime
              tract(idime) = tract(idime) + acvis*e%baloc(jdime,ndime)*(grbve(idime,jdime) + grbve(jdime,idime))
          end do
      end do
   end subroutine

   subroutine TractionGauss
      implicit none
      integer(ip) :: inodb,inode
      do inodb = 1,e%pnodb
         inode= e%lboel(inodb)
         do idime=1,e%ndime
            traction(idime,inode) = traction(idime,inode) + tract(idime)*e%shapb(inodb,e%igaub)*dsurf
         end do
         elvmass(inode) = elvmass(inode) + e%shapb(inodb,e%igaub)*dsurf
      end do
   end subroutine

   subroutine AssemblyTraction
      implicit none
      call a%Mesh%AssemblyToArray(e,1_ip,elvmass,vmassBoundary)
      call a%Mesh%AssemblyToArray(e,e%ndime,traction,a%btraction)
   end subroutine
   
   subroutine ExtraOp
      implicit none
      integer(ip) :: kfl_HangingNodes, kfl_perio
      integer(ip) :: ipoin,npoin,ndime

      call a%Mesh%GetNpoin(npoin)     
      call a%Mesh%GetHanging(kfl_HangingNodes)
      call a%Mesh%GetPerio(kfl_perio)
      call a%Mesh%GetNdime(ndime)
      do ipoin = 1,npoin
         if (vmassBoundary(ipoin) > 0.0_rp) a%btraction(:,ipoin) = a%btraction(:,ipoin)/vmassBoundary(ipoin)
      enddo
      call a%Mesh%ArrayCommunicator%GhostCommunicate(ndime,a%btraction)
      !Hanging nodes
      if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(ndime,a%btraction)
      if (kfl_perio == 1) call a%Mesh%MasterToSlave(ndime,a%btraction)
   end subroutine

end module
