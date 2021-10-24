module Mod_nsf_bouope_1st
   use typre
   use Mod_nsm_BaseElmope
   use Mod_nsm_BaseBouope
   use Mod_nsm_NeumannBC
   use Mod_nsm_ViscosityBouope
  implicit none

contains
   subroutine SetBCPointers
      implicit none
      call ResetProcedureComposition
      
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetPointersAndHooksToNULLSUB
      call SetPointersNeumannBC%Initialize
      call SetPointersViscosity%Initialize
      call SetPointersNeumannBC%Set
      call SetPointersViscosity%Set
      call SetPointersNeumannBC%Finalize
      call SetPointersViscosity%Finalize
   end subroutine
end module

subroutine nsf_bouope_1st(b)
   use typre
   use Mod_Element
   use Mod_NavierStokes
   use Mod_NSFractionalStep
   use Mod_nsm_elmdir
   use Mod_php_elmdir
   use Mod_nsm_BaseElmope
   use Mod_nsm_HangingNodes
   use Mod_nsf_bouope_1st
   implicit none
   class(NSFractionalStepProblem), target :: b
   
   real(rp), allocatable :: welmat(:,:,:,:)
   real(rp), allocatable :: welrhs(:,:)
   
   integer(ip) :: inode,jnodb,inodb,ipoin,jnode,jpoin
   integer(ip) :: kfl_HangingNodes
   
   interface
       subroutine nsf_ExtrapolatePressure(a)
         use Mod_NSFractionalStep
         implicit none
         class(NSFractionalStepProblem), target   :: a
       end subroutine
    end interface

    a=> b

    !Set Pointers
    call SetBCPointers
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,welmat,'welmat','nsf_bouope_1st')
   call a%Memor%alloc(e%ndime,e%mnode,welrhs,'welrhs','nsf_bouope_1st')
   call AllocateBaseBouopeMatrices
   call AllocateBaseBouopeArrays
   
   call a%Mesh%GetHanging(kfl_HangingNodes)
   
   call ProcHook_Initializations
   
   !Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      call BouopeElMatsToZero 
      call ProcHook_ElmatsToZero

      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)
      call ProcHook_PhysicalProp

      call BoundaryGathers
      call ProcHook_Gathers
      
      !Cartesian derivatives
      call e%elmdel
      !Element length
      call e%elmlen
      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,bovel,chale,a%kfl_hdifumin)
      
      dsurf0 = 0.0_rp
      
      call ProcHook_PreGauss
      
      !Gauss-Point Loop
      GaussPointLoop: do igaub=1,e%pgaub
         e%igaub = igaub
      
         !Derivatives at the boundary
         call e%elmderb
      
         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         wmatr = 0.0_rp
         tract = 0.0_rp
         
         call BoundaryInterpolates 
         call ProcHook_Interpolates
         call ProcHook_InGauss

         !Velocity norm
         call vecnor(gpbve(:,1),e%ndime,gpvno,2)
         
         call ProcHook_InGaussElmats

         call ProcHook_InGaussElmatsAssembly

      end do GaussPointLoop

      call ProcHook_PostLoop
  
      !At this point we do nothing with the pressure terms (either equation and unknown)
      welmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
      welrhs(1:e%ndime,1:e%pnode) = elrhs(1:e%ndime,1:e%pnode) 

      if (kfl_HangingNodes == 1) call ModifyMatricesHanging
      
      !Boundary conditions
      call nsm_rotdir(a,e,e%ndime,welmat,welrhs)
      call php_elmdir(a,e,e%ndime,a%ndofbc,a%ndofbcstart,1_ip,welmat,welrhs)

      !Assembly
      call a%LinearSystem%Assembly(e,welmat,welrhs)
      
   end do boundaries
   
   call ProcHook_Finalizations

   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,welmat,'welmat','nsf_bouope_1st')
   call a%Memor%dealloc(e%ndime,e%mnode,welrhs,'welrhs','nsf_bouope_1st')
   call DeallocateBaseBouopeArrays
   call DeallocateBaseBouopeMatrices
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')

end subroutine nsf_bouope_1st
