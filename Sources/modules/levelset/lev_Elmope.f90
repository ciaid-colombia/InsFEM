

module Mod_LevElmope
   use typre
   use Mod_lev_BaseElmope
   use Mod_lev_HangingNodes
   use Mod_lev_Advection

contains

   !SetPointers
   subroutine SetPointers
      use typre 
      implicit none
      
      integer(ip) :: kfl_nonlinear,nelty    
      
      !External Procedures
      procedure() :: NULLSUB  


      call ResetProcedureComposition
      
      !--------------------------------------------------------------------
      !Defaults
      
      !Pointers
      ProcPointer%tem_elmbuv => tem_elmbuv
      ProcPointer%tem_elmrhu => tem_elmrhu
      ProcPointer%ComputeTau => ComputeStaticTau
      ProcPointer%ExternalForces => ExternalForces
      
      !Hooks
      ProcHook%Initializations => NULLSUB
      ProcHook%Gathers         => NULLSUB
      ProcHook%OnIeltyChange   => NULLSUB
      ProcHook%PreGauss        => NULLSUB
      ProcHook%Interpolates    => NULLSUB
      ProcHook%Elext           => NULLSUB
      ProcHook%InGauss         => NULLSUB
      ProcHook%InGaussElmats   => NULLSUB
      ProcHook%PostGaussElmats => NULLSUB
      ProcHook%Testf           => NULLSUB
      ProcHook%PreDirichlet    => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      ProcHook%PhysicalProp    => NULLSUB
      ProcHook%PreAssembly     => NULLSUB
      
      call SetPointersHangingNodes(0)
      call SetPointersAdvection(0)
      
      !--------------------------------------------------------------------
      !Advection (nonsense without convective term)
      call SetPointersAdvection(1)
      
      !-------------------------------------------------------------------
      !Non-Linear Elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call ConcatenateProcedures(ProcHook%InGauss,InGaussNonLinear)
      endif 
      !------------------------------------------------------
      !Hanging Nodes
      call SetPointersHangingNodes(1)
      
      
      
      
      
      
      call SetPointersHangingNodes(100)
      call SetPointersAdvection(100)
      
   end subroutine

   !For computing Tau
   subroutine ComputeStaticTau
      ! acvis=0.0, acrcp=0.0 and acden=1.0
      call ComputeTauCDR(e,acden,acvis,acrcp,gpvno,a%staco,chale,tilev)
      
   end subroutine  
   
  
   !---------------------------------------------------------------
   !For external Forces
   subroutine ExternalForces
      implicit none

   end subroutine
   
   subroutine InterpolateCoord
      implicit none
      
      call e%interpg(e%ndime,e%elcod,gpcod)
   end subroutine   
   
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
   

end module


subroutine lev_elmope(LevelProblem)
   use typre
   use Mod_LevelSet
   use Mod_lev_BaseElmope
   use Mod_LevElmope
   use Mod_TemperatureElement
   implicit none
   class(LevelSetProblem), target :: LevelProblem
   logical :: isALE
   a=>LevelProblem
   
   !If we need to force NoALE, do it
   if (a%kfl_ForceEulerianAdvection == 1) then
      call a%Mesh%GetALE(isALE)
      call a%Mesh%SetALE(0)
   endif
   
   !call a%FilePostpr%postpr(a%veloc,'lev_veloc',a%istep,a%ctime,a%Mesh)
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   ReferenceDtinv = a%dtinv
   if (a%kfl_tsche_1st_current == 'CN   ') ReferenceDtinv = a%dtinv*2    
      
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_elmope')

   !Matrices Alloc
   call a%Memor%alloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','lev_elmope')
   call a%Memor%alloc(1_ip,e%mnode,elrhs,'elrhs','lev_elmope')
   
   !Other Arrays Alloc
   call a%Memor%alloc(e%mnode,a%ncomp-1,ellev,'ellev','lev_elmope')
   call a%Memor%alloc(a%ncomp-1,gplev,'gplev','lev_elmope')
   call a%Memor%alloc(e%mnode,testf,'testf','lev_elmope')
   !Advection
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','lev_elmope')
   call a%Memor%alloc(e%ndime,gpvel,'gpvel','lev_elmope')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','lev_elmope')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','lev_elmope')
   gpvno = 0.0_rp
   
   !Hook
   call ProcHook%Initializations
   
   !do itest = 1,100
   call a%Mesh%GetNelem(nelem) 
   
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)      
      !Hook
      call ProcHook%OnIeltyChange
      
      !Hook
      call ProcHook%PreGauss
      
      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp
      
      !Gathers
      call e%gather(1,ellev(:,1),a%level(:,1))
      call e%gather(1,ellev(:,2),a%level(:,3))
      do itime = 3,nsteps ! Time bdf2 and others
         call e%gather(1,ellev(:,itime),a%level(:,itime+1)) 
      enddo
      
      !Hook
      call ProcHook%Gathers
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      
      !Element length at center of gravity
      call e%elmlen

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elvel,chale)
      
      !DvolsToZero
      dvolt0=0.0_rp
      dvolt1=0.0_rp
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
         e%igaus = igaus
         
         !Hook
         call ProcHook%InGauss

         dvol = e%weigp(e%igaus)*e%detjm
         
         !Interpolate
         call e%interpg(1,ellev(:,1),gplev(1))
         call e%interpg(1,ellev(:,2),gplev(2))
         do itime = 3,nsteps ! Time bdf2 and others
            call e%interpg(1,ellev(:,itime),gplev(itime))
         enddo
         !Hook
         call ProcHook%Interpolates
         
         !Advection velocity norm
         call vecnor(gpadv,e%ndime,gpvno,2)
      
         !Compute a·grad(V)
         call ComputeAGradV(e,gpadv,AGradV)
         
         !PhysicalProperties for use the CDR options
         acvis = 0.0_rp
         acrcp = 0.0_rp
         acden = 1.0_rp
         !Hook
         call ProcHook%PhysicalProp
         
         !Compute the stability parameters
         call ProcPointer%ComputeTau
         
         !Adjoint Test Function
         !Stabilization terms : -tau L*v
         call tem_ComputeTestf(e,acden,tilev,AGradV,acvis,testf)
         !Hook
         call ProcHook%Testf
         
         !Volumes Times Taus
         dvolt0=dvol       + dvolt0                ! w(gp)*detjm
         dvolt1=dvol*tilev + dvolt1                ! w(gp)*detjm*tau1
         
         !Compute Elext, Temporal Derivatives
         elext=0.0_rp
         !Time integration
         call tem_TimeIntegrationToElext(e,Integrator,acden,a%dtinv,gplev,elext)
         !Hook
         call ProcHook%Elext
         
         !InGaussElmats
         !Compute contributions to RHS : Block U
         call ProcPointer%tem_elmrhu(e,dvol,testf,elext,elrhs)
               
         !Compute contributions to elemental matrix : Block U,V
         call ProcPointer%tem_elmbuv(e,dvol,acden,acrcp,LHSdtinv,AGradV,testf,elmat(1,:,1,:))
         
         !Hook
         call ProcHook%InGaussElmats
   
      enddo gauss_points      

      !PreDirichlet
      call ProcHook%PreDirichlet
      
      !Dirichlet Boundary Conditions
      call php_elmdir(a,e,1_ip,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
      
      !PreAssembly
      call ProcHook%PreAssembly
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo elements
   
   !Finalizations
   !Hook
   call ProcHook%Finalizations
   
   !Matrices Alloc
   call a%Memor%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','lev_elmope')
   call a%Memor%dealloc(1_ip,e%mnode,elrhs,'elrhs','lev_elmope')
   
   !Other Arrays Alloc
   call a%Memor%dealloc(e%mnode,a%ncomp-1,ellev,'ellev','lev_elmope')
   call a%Memor%dealloc(a%ncomp-1,gplev,'gplev','lev_elmope')
   call a%Memor%dealloc(e%mnode,testf,'testf','lev_elmope')
   !Advection
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','lev_elmope')
   call a%Memor%dealloc(e%ndime,gpvel,'gpvel','lev_elmope')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','lev_elmope')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','lev_elmope')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_elmope')
   
   !unforce NoALE for the rest of the operations and modules
   if (a%kfl_ForceEulerianAdvection == 1) then
      if (isALE .eqv. .true.) then
         call a%Mesh%SetALE(1_ip)
      else
         call a%Mesh%SetALE(0_ip)
      endif
   endif
   
end subroutine
