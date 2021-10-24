module Mod_ale_Elmope
   use typre
   use Mod_ale_BaseElmope
   implicit none
   
contains
   subroutine SetPointers
      use Mod_ale_FixedMesh
      implicit none
      
      
      procedure() :: NULLSUB
   
      call ResetProcedureComposition
   
      ProcHook%Initializations => NULLSUB
      ProcHook%PostGaussElmats => NULLSUB
      ProcHook%PreDirichlet    => NULLSUB
      ProcHook%Finalizations   => NULLSUB
      
      !Pointers Initialization
      call SetPointersFixedMesh%Initialize
      
      
      !Set the Pointers
      call SetPointersFixedMesh%Set
      
      
      !Finalize the pointer setter
      call SetPointersFixedMesh%Finalize
   end subroutine



end module

subroutine ale_elmope(b,externalcurrentbvess)
   !-----------------------------------------------------------------------
   !> This routine performs an element loop and computes the elemental matrix 
   !! and RHS for the mesh displacements, applying boundary conditions. The free
   !! surface condition (kfl_fixno=3) is first changed to 1 so it can enter the 
   !! php_elmdir, where the Dirichlet boundary conditions are applied. Later on, 
   !! they are reset to 3.
   !-----------------------------------------------------------------------
   use typre
   use Mod_ConvectiveElement  
   use Mod_php_elmdir
   use Mod_Alemov
   use Mod_ale_BaseElmope
   use Mod_ale_Elmope
   
   implicit none
   class(AlemovProblem), target :: b
   integer(ip) :: externalcurrentbvess, poinStatus
   
   a=> b
   currentbvess = externalcurrentbvess
   
   !SetPointers
   call SetPointers
   
   !Call the Mesh to allocate element
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','ale_elmope')
   !Matrices Alloc
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','ale_elmope')
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','ale_elmope')
   

   !Initializations
   call ProcHook%Initializations
   
   
   
   !For a purely Lagrangian Free Surface
   call a%Mesh%GetNpoin(npoin)
   do ipoin=1,npoin                              ! 4 is used in FSI where disp are prescribed
      if (a%kfl_fixno(currentbvess,ipoin)==3 .or. a%kfl_fixno(currentbvess,ipoin)==4) then 
         a%kfl_fixno(currentbvess,ipoin)=1
      end if
   end do

   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e) 

      !ElmatsToZero
      elmat=0.0_rp
      elrhs=0.0_rp

      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg

      dvol = sum(e%weigp(1:e%pgaus))*e%detjm

      !Viscosity terms : we only consider mu*(grad v, grad u)
      call elmvis(e,dvol,1.0_rp,elmat)
      
      !PostGaussElmats
      call ProcHook%PostGaussElmats
      
      !PreDirichlet
      call ProcHook%PreDirichlet
      
      !Dirichlet Boundary Conditions
      call php_elmdir(a,e,a%ndofn,1_ip,a%ndofbcstart,currentbvess,elmat,elrhs)
      
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo elements

   a%kfl_fixno(currentbvess,1:npoin)=a%kfl_fixno0(currentbvess,1:npoin)

   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','ale_elmope')
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','ale_elmope')
   
   !Finalizations
   call ProcHook%Finalizations
   
   
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','ale_elmope')
end subroutine
