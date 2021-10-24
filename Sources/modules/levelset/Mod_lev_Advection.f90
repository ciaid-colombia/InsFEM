module Mod_lev_Advection
   use typre
   use Mod_lev_BaseElmope
   use Mod_LevelSet
   implicit none
   private
   public SetPointersAdvection
   
   integer(ip), allocatable :: kfl_IsSet
   
   !ALE
   logical :: isALE
   real(rp), allocatable :: elmve(:,:),gpmve(:)
   real(rp), pointer :: meshve(:,:,:) => NULL()
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersAdvection(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            !--------------------------------------------------------------------
            !Advection (nonsense without convective term)
            if(a%kfl_advec /=0) then
               if(a%kfl_advec == 1)then
                  call ConcatenateProcedures(ProcHook%Gathers,GatherExternalVeloc)
                  call ConcatenateProcedures(ProcHook%Interpolates,VelocInterpolates) 
               else
                  call ConcatenateProcedures(ProcHook%Interpolates,GaussPointAnalyticalVeloc)
               end if
               
               !call ConcatenateProcedures(ProcHook%PreAssembly,EnsureEntryBoundaryConditions)
               call ConcatenateProcedures(ProcHook%InGaussElmats,EnsureEntryBoundaryConditions)
               
            end if
            
            !ALE
            call a%Mesh%GetALE(isALE)
            if (isALE .and. a%kfl_ForceEulerianAdvection == 0) then
               call a%Mesh%GetMeshVeloc(meshve)
               call ConcatenateProcedures(ProcHook%Initializations,AllocMeshVeloc)
               call ConcatenateProcedures(ProcHook%Gathers,GatherMeshVeloc)
               call ConcatenateProcedures(ProcHook%Interpolates,InterpolateMeshVeloc)
               call ConcatenateProcedures(ProcHook%Interpolates,ALEAdvectionVelocity)
               call ConcatenateProcedures(ProcHook%Finalizations,DeallocMeshVeloc)
            end if
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   



   !--------------------------------------------------------------------
   !FOR ADVECTION
   subroutine GatherExternalVeloc
      real(rp), pointer :: coord(:) => NULL()
      real(rp) :: dist(2)
      
      !Gathers
      call e%gather(e%ndime,elvel(:,:),a%veloc(:,:))
   end subroutine
   
   subroutine VelocInterpolates
      !Interpolate
      call e%interpg(e%ndime,elvel(:,:),gpvel(:))
      gpadv(1:e%ndime) = gpvel(1:e%ndime)
   end subroutine  
   
   !to evaluate the problem without Navier-Stokes Coupling
   subroutine GaussPointAnalyticalVeloc
      implicit none
      
      call tem_elmvel(e,a%kfl_advec,gpvel)
      gpadv = gpvel
      
   end subroutine
   
   
   !ALE
   subroutine AllocMeshVeloc
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elmve,'elmve','nsm_AdvectionVelocity')
      call a%Memor%alloc(e%ndime,gpmve,'gpmve','nsm_AdvectionVelocity')
   end subroutine
   
   subroutine GatherMeshVeloc
      implicit none
      
      call e%gather(e%ndime,elmve,meshve(:,:,1))
   end subroutine
   
   subroutine InterpolateMeshVeloc
      implicit none
      
      call e%interpg(e%ndime,elmve,gpmve)
   end subroutine
   
   subroutine ALEAdvectionVelocity
      implicit none
      
      gpadv = gpadv-gpmve
   end subroutine
   
   subroutine DeallocMeshVeloc
      implicit none
      
      call a%Memor%dealloc(e%ndime,e%mnode,elmve,'elmve','nsm_AdvectionVelocity')
      call a%Memor%dealloc(e%ndime,gpmve,'gpmve','nsm_AdvectionVelocity')
   end subroutine
   
   subroutine EnsureEntryBoundaryConditions
      implicit none
      
      real(rp), pointer :: exnor(:,:) => NULL()
      real(rp) :: adiag
      integer(ip) :: inode,ipoin,ibopo,jnode,idime
      real(rp) :: prod,raux,vnor,gradient(e%ndime), gnor
      
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if (ibopo /= 0) then
            prod = dot_product(exnor(:,1),elvel(:,inode))
            
            !Nodal velocity norm
            call vecnor(elvel(:,inode),e%ndime,vnor,2)
            vnor = vnor*10.0_rp !We ensure that diffusivity dominates
            !If we have inflow advection velocity, fix the level set value
            do idime = 1,e%ndime
               gradient(idime) = dot_product(e%cartd(idime,1:e%pnode),ellev(1:e%pnode,1))
            enddo
            call vecnor(gradient,e%ndime,gnor,2)
            !We want the new gradient to have the same direction and norm 1
            if (prod < -1e-6 .and. gnor > 1e-6) then
               do jnode = 1,e%pnode
                  raux = dot_product(e%cartd(1:e%ndime,inode),e%cartd(1:e%ndime,jnode))*dvol*vnor
                  elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + raux
                  elrhs(1,inode) = elrhs(1,inode) + raux*ellev(jnode,1)/gnor
               enddo     
               
               !Dirichlet boundary conditions
               !adiag=elmat(1,inode,1,inode)
               !elrhs(1,inode)=adiag*ellev(inode,1)
               !elmat(1,inode,:,1:e%pnode)=0.0_rp
               !elmat(1,inode,1,inode)=adiag
            end if
         endif
      end do
   end subroutine
         
            
   

end module
