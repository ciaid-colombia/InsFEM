module Mod_tem_SigmaupCoupling
   use Mod_tem_BaseElmope
   use Mod_tem_ExternalForces
   use typre
   implicit none
   private
   public :: SetPointersSigmaupCoupling
   real(rp), allocatable :: sigNS(:,:)
   real(rp), allocatable :: gpsigNS(:)
   integer(ip) :: ndime, npoin
   integer(ip), allocatable :: kfl_IsSet 
   real(rp), allocatable, target :: elSmoothVelGradient(:,:,:),auxGradDisp(:,:)
   real(rp),allocatable :: elViscoDiss(:)
   real(rp)             :: gpViscoDiss

contains

   subroutine SetPointersSigmaupCoupling(itask,task)
      implicit none
      integer(ip) :: itask
      character(6)::task
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            !Pointers are already set
            kfl_IsSet = 1
            if (a%kfl_CouplingThreeField == 1) then
               call ConcatenateProcedures(ProcHook%Initializations,AllocSup)
               call ConcatenateProcedures(ProcHook%Gathers,GathersSigmaTerm)
               call ConcatenateProcedures(ProcHook%Interpolates,InterpolateSigmaNS)
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocSup)
               call ConcatenateProcedures(ProcPointer%ExternalForces,SigmaTerm) 
            endif
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-------------------------------------------------------------------
   !Computation Subroutines
   subroutine AllocSup
      implicit none
      call a%Memor%alloc(e%mnode,elViscoDiss,'elViscoDiss','tem_elmope')
      call a%Memor%alloc((e%ndime-1)*(e%ndime-1)+2,e%mnode,sigNS,'sigNS','tem_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elSmoothVelGradient,'elSmoothVelGradient','ComputeDisplacementGradients')
      call a%Memor%alloc((e%ndime-1)*(e%ndime-1)+2,gpsigNS,'gpsigNS','tem_elmope')
      call a%Memor%Alloc(e%ndime,e%ndime,auxGradDisp,'auxGradVeloc','tem_elmope')
   end subroutine
   
   subroutine DeallocSup
      implicit none
      call a%Memor%dealloc(e%mnode,elViscoDiss,'elViscoDiss','tem_elmope')
      call a%Memor%dealloc((e%ndime-1)*(e%ndime-1)+2,e%mnode,sigNS,'sigNS','tem_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elSmoothVelGradient,'elSmoothVelGradient','ComputeDisplacementGradients')
      call a%Memor%dealloc(e%ndime,e%ndime,auxGradDisp,'auxGradVeloc','tem_elmope')
      call a%Memor%dealloc((e%ndime-1)*(e%ndime-1)+2,gpsigNS,'gpsigNS','tem_elmope')
   end subroutine
   

   subroutine GathersSigmaTerm
      call e%gather((e%ndime-1)*(e%ndime-1)+2,sigNS(:,:),a%sigmaNS(:,:))
      call e%gather(e%ndime*e%ndime,elSmoothVelGradient,a%SmoothedVelocityGradient)
   end subroutine
   
   subroutine InterpolateSigmaNS
      implicit none
      call e%interpg((e%ndime-1)*(e%ndime-1)+2,sigNS,gpsigNS)
      call e%interpg(e%ndime*e%ndime,elSmoothVelGradient,auxGradDisp)
   end subroutine  
   
   subroutine SigmaTerm
      implicit none
        call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
        call tem_SigmaTerm(e,gpsigNS,auxGradDisp,acvis,acden,acsph,elext,gpViscoDiss) !!! ojo modificar
        if (a%npp_stepi(9) /= 0) a%sigmatermarray(ielem)%a(e%igaus)=gpViscoDiss
   end subroutine  
   
   subroutine UpdateSpecificHeatAndConductivity
      implicit none
      real(rp) :: c_p0, c_ps, ks, k0
      !Only for ºC
      c_p0  = 1.2122
      c_ps  =-0.00112
      ks    = 0.00118
      k0    = 0.7753
      
!       actco= a%tcond*(k0 + ks*gptem(1)) 
!       acsph= a%sphea*(c_p0 + c_ps*gptem(1))
   end subroutine
   !------------------------------------------------------
   
   subroutine GaussPointAssembly
      implicit none
      integer:: inode
      do inode = 1,e%pnode
         elViscoDiss(inode) = elViscoDiss(inode) + e%shape(inode,e%igaus)*gpViscoDiss*dvol
      enddo
   end subroutine  
   
   
   subroutine AssemblyViscousDissipation
      implicit none
       !call a%Mesh%AssemblyToArray(e,size(a%ViscousDissipation,1),elViscoDiss, a%ViscousDissipation)
      a%ViscousDissipation(e%lnods(1:e%pnode)) = a%ViscousDissipation(e%lnods(1:e%pnode)) + elViscoDiss(1:e%pnode)    
   end subroutine
   
   
   subroutine tem_SigmaTerm(e,sig,grvel,visco,dens,sph,elext,term)
      !Computes the source term from the three field problem sigma:grad_sym(u)
      use Mod_Element
      use typre
      implicit none
      class(FiniteElement) :: e
      real(rp), intent(in) :: sig((e%ndime-1)*(e%ndime-1)+2), grvel(e%ndime,e%ndime)         
      real(rp), intent(in) ::  dens, sph, visco
      real(rp), intent(out):: term
      real(rp), intent(inout) :: elext
      integer(ip) :: idime, jdime
      real(rp)    :: sigcomp(e%ndime,e%ndime), grvelSym(e%ndime,e%ndime)
   
      term=0.0_rp

      call sup_SigmaMatrix(e%ndime,(e%ndime-1)*(e%ndime-1)+2,sig,sigcomp)
      grvelSym=0.5*(grvel+transpose(grvel))
      
      do idime=1,e%ndime
         do jdime=1,e%ndime
            term= term+(sigcomp(idime,jdime))*grvelSym(idime,jdime)
         end do
      end do
      !2*(0.0667)*visco*grvelSym(idime,jdime)
   
      term=term/(sph*acden) !ojo dens fuera de la división!!***
      
      elext=term+elext
      
   end subroutine

      
end module
