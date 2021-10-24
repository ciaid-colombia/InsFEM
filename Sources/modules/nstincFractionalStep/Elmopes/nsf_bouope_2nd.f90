module Mod_nsf_bouope_2nd
   use Mod_nsm_baseElmope
   use Mod_nsm_baseBouope
   use Mod_nsm_HangingNodes
   implicit none 
   
contains

   subroutine SetPointers
      implicit none
      
      integer(ip) :: nelty
      
       !External Procedures
      procedure() :: NULLSUB
      
      call ResetProcedureComposition
      
      ProcHook_Initializations => NULLSUB
      ProcHook_OnIeltyChange => NULLSUB
      ProcHook_ElmatsToZero => NULLSUB
      ProcHook_PreGauss => NULLSUB
      ProcHook_Interpolates => NULLSUB
      ProcHook_InGaussElmats=> NULLSUB
      ProcHook_PreDirichlet => NULLSUB
      ProcHook_Finalizations => NULLSUB
      
      call SetPointersHangingNodes%Initialize
      
      call SetPointersHangingNodes%Set
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook_OnIeltyChange => OnIeltyChange
      
      call SetPointersHangingNodes%Finalize

   end subroutine
   
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
end module

subroutine nsf_bouope_2nd(b)
   use typre
   use Mod_Element
   use Mod_NSFractionalStep
   use Mod_php_SetTimeIntegrator 
   use Mod_TimeIntegrator   
   use Mod_NSF_Element
   use Mod_nsf_bouope_2nd
   implicit none
   class(NSFractionalStepProblem), target :: b

   real(rp), allocatable :: bopre1(:,:)
   real(rp), allocatable :: gppre1(:)
   integer(ip) :: inodb,inode,istep,jnode,jnodb,npoin
   real(rp) :: raux,ValRHS(1)
   
   a => b
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)

   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmat,'elmat','nsm_bouop0')
   call a%Memor%alloc(1,e%mnode,elrhs,'elrhs','nsm_bouop0')
   
   call a%Memor%alloc(e%mnode,1,elpre,'elpre','nsm_bouop0')
   call a%Memor%alloc(ndime,e%mnodb,bovel,'bovel','nsm_bouop0')
   call a%Memor%alloc(e%mnodb,nsteps,bopre1,'bopre1','nsm_bouop0')
   call a%Memor%alloc(1,e%ndime,grpre,'grpre','nsm_bouop0')
   call a%Memor%alloc(nsteps,gppre1,'gppre1','nsm_bouop0')
   
   call a%Memor%alloc(e%ndime,1,gpbve,'gpbve','nsf_bouope_2nd')
   
   !SetPointers
   call SetPointers
   
   
   !Hook
   call ProcHook_Initializations
   
   
   !Physical Parameters
   call a%GetPhysicalParameters(imat,acden,acvis)
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !Hook
      call ProcHook_OnIeltyChange
      
      !Initialize
      elmat=0.0_rp
      elrhs=0.0_rp

      call ProcHook_ElmatsToZero
      
      !Gathers
      call e%gather(1,elpre,a%press(:,3)) ! p_n
      call e%gatherb(e%ndime,bovel,a%veloc(:,:,1)) ! int_u_n+1

      call e%gatherb(1,bopre1(:,1),a%press(:,1))
      do istep=2,nsteps
         call e%gatherb(1,bopre1(:,istep),a%press(:,istep+1))
      enddo

      call e%elmdel
      
      !Hook
      call ProcHook_PreGauss
   
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         !Calculate exterior Normal
         call e%bounor
         
         !Derivatives at the boundary
         call e%elmderb
         
         !Interpolate
         call e%interpb(e%ndime,bovel,gpbve)

         !Hook
         call ProcHook_Interpolates
         
         dsurf=e%weigb(e%igaub)*e%eucta
         
         
         if (a%kfl_fixbo(iboun) == 6) then
            raux = 1.0_rp/(acden*LHSdtinv)
            
            !Assembly open boundary to LHS (v,n*gradp)
            do inodb = 1,e%pnodb
               inode = e%lboel(inodb)
               do jnodb = 1,e%pnodb
                  jnode = e%lboel(jnodb)
                  elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) - raux*e%shapb(inodb,e%igaub)*dot_product(e%cartb(1:e%ndime,jnode),e%baloc(1:e%ndime,e%ndime))*dsurf
               enddo
            enddo            
         endif
         
         if (a%kfl_fixbo(iboun)==6) then

            !We are in second order fractional step, assembly also contribution from pressure at the previous step
            call e%gradientb(1_ip,elpre,grpre)
            raux = 1.0_rp/(acden*LHSdtinv)           
            do inodb = 1,e%pnodb
               inode = e%lboel(inodb)
               elrhs(1,inode) = elrhs(1,inode) + raux*e%shapb(inodb,e%igaub)*dot_product(grpre(1,1:e%ndime),e%baloc(1:e%ndime,e%ndime))*dsurf
            enddo
         endif
         
         
         !Hooks
         call ProcHook_InGaussElmats
         
      enddo
      
      !Hook
      call ProcHook_PreDirichlet
   
      !Dirichlet Boundary conditions
      call nsf_elmdir_press(b,e,elmat,elrhs)

      ! Assembly
      call b%LinearSystemP%Assembly(e,elmat,elrhs)
      
   end do boundaries
   
   !Hook
   call ProcHook_Finalizations
   
   call a%Memor%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','nsm_bouop0')
   call a%Memor%dealloc(1_ip,e%mnode,elrhs,'elrhs','nsm_bouop0')
  
   call a%Memor%dealloc(e%mnode,1,elpre,'elpre','nsm_bouop0')
   call a%Memor%dealloc(ndime,e%mnodb,bovel,'bovel','nsm_bouop0')
   call a%Memor%dealloc(e%mnodb,nsteps,bopre1,'bopre1','nsm_bouop0')
   call a%Memor%dealloc(1,e%ndime,grpre,'grpre','nsm_bouop0')
   call a%Memor%dealloc(nsteps,gppre1,'gppre1','nsm_bouop0')
  
  
   call a%Memor%dealloc(e%ndime,1,gpbve,'gpbve','nsf_bouope_2nd')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   
   
end subroutine
