subroutine lmn_bouop0(LMProblem)
!-----------------------------------------------------------------------
!****f* lmachn/lmn_bouope
! NAME 
!    lmn_bouop0
! DESCRIPTION
!    Low Mach number boundary elemental operations:
!    1. Compute elemental matrix and RHS 
!    2. Impose Dirichlet boundary conditions
!    3. Assemble them
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element 
   use Mod_LowMach
   use Mod_lmn_elmdir
   use Mod_lmn_BaseElmope
   use Mod_lmn_HangingNodes
   implicit none
   class(LowMachProblem), target :: LMProblem
   integer(ip) :: iboun,nboun,ndime
   integer(ip) :: kfl_HangingNodes
   integer                  :: ierr,ibopo,ipoin,inodb
   real(rp) :: vnor
   real(rp), pointer :: exnor(:,:)=> NULL()
   logical  :: cycleflag
   
   interface 
      subroutine lmn_bouope(a,iboun,e,elmat,elrhs)
         !This subroutine computes the boundary contribution
         !-----------------------------------------------------------------------
            use typre
            use Mod_LowMach
            use Mod_Element
            implicit none
            class(LowMachProblem)  :: a
            integer(ip)            :: iboun
            class(FiniteElement)   :: e
            real(rp),intent(inout) :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode)              ! Element matrices
            real(rp),intent(inout) :: elrhs(a%ndofn,e%mnode)
      end subroutine
   end interface
   
   !-----------------------------------------------------------

   a=>LMProblem

   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   
   call a%Mesh%GetHanging(kfl_HangingNodes)
   
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','lmn_bouop0')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','lmn_bouop0')
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      cycleflag = .false.
      do inodb = 1,e%pnodb
         ipoin = e%lnodb(inodb)
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if (ibopo == 0) then
            cycleflag =  .true.
         else
            call vecnor(exnor(:,1),e%ndime,vnor,2)
            if (vnor == 0.0_rp) cycleflag =  .true. 
         end if
      end do

      if (cycleflag) cycle

      
      !Initialize
      elmat=0.0_rp
      elrhs=0.0_rp
      
      !Compute boundary terms
      call lmn_bouope(a,iboun,e,elmat,elrhs)
      
      if (kfl_HangingNodes == 1) call ModifyMatricesHanging

      !Boundary conditions
      call lmn_elmdir(a,e,elmat,elrhs)

      ! Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   end do boundaries

  
   
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','lmn_bouop0')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','lmn_bouop0')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
end subroutine lmn_bouop0
