subroutine plcd_EndBouope(b)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!  Calculates tractions in the boundary to check convergence for the 
!  FSI loops
!
!-----------------------------------------------------------------------
   use typre 
   use Mod_Mesh
   use Mod_Element
   use Mod_memor
   use Mod_iofile
   use Mod_int2str
   use Mod_PLCD  
   use Mod_plcd_BaseElmope
   use Mod_plcd_BMatrix
   use Mod_plcd_BMatrixFactory
   use MPI
   
   implicit none
   class(PLCDProblem), target :: b   
   
   integer(ip) :: idime,iboun,igaub,nboun,ndime,inodb,npoin,inode,ipoin
   real(rp)              :: dsurf
   real(rp), allocatable :: tract(:),StressTensor(:,:)
   class(PLCDMaterial), pointer :: Material

!    if (a%kfl_FSI == 1) then
!    
!    a => b
! 
!    !Initializations
!    call a%Mesh%GetNelem(nelem)
!    call a%Mesh%GetNdime(ndime)
!    call a%Mesh%GetNboun(nboun)       
!    call a%Mesh%GetNpoin(npoin)     
!    
!    call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_EndBouope')      
!    
!    !Allocation
!    call a%Memor%alloc(e%ndime,e%mnode,3,eldisp,'eldisp','plcd_EndBouope')
!    call a%Memor%alloc(e%ndime,e%ndime,StressTensor,'StressTensor','plcd_EndBouope')
!    call a%Memor%alloc(ndime,tract,'tract','plcd_EndBouope')      
! 
!    !Inializations
!    a%btraction=0.0_rp
!    
!    call CreateBMatrix(a,a%Memor,BMat)
!    call BMat%Alloc(e,a%Memor)
!    
!    !GradDispHistory might be different than gradDisp, for instance if smoothed gradients are used
!    gradDispHistory => gradDisp
!    
!    ! Loop over boundaries
!    boundaries: do iboun=1,nboun
!       !Load Element
!       call a%Mesh%BoundaryLoad(iboun,e)     
! 
!       ElementMatData => a%ElementMaterialsData(ielem)%p
!       call ElementMatData%GetMaterialPointer(Material)
! 
!       !Gathers in boundary elements 
!       call e%gather(e%ndime,eldisp(:,:,1),a%Displacement(:,:,1))
!       
!       !Compute linear derivatives
!       call e%elmdel
! 
!       !Gauss-Point Loop
!       GaussPointLoop : do igaub=1,e%pgaub
!           e%igaub = igaub
! 
!           !Calculate exterior Normal
!           call e%bounor
! 
!           dsurf=e%weigb(e%igaub)*e%eucta
! 
!           !Derivatives at the boundary
!           call e%elmderb         
!           
!           !ComputeDisplacementGradients
!           call e%gradient(e%ndime,eldisp(:,:,1),gradDisp)
! 
!           !Stress calculation
!           call BMat%Setup(e)
!           call ElementMatData%ComputeHistoryAndConstitutiveTensor(e%igaus,gradDispHistory)
!           call ElementMatData%ComputeStress(e%igaus,gradDisp)
!           call ElementMatData%GetStressTensorPointer(e%igaus,stress)
!           call Material%CT%GetStressTensor(stress,StressTensor)
! 
!           !----------------------------------------------------------------
!           !Traction value
!           !  t=n*S  S=depends on constitutive model
!           do idime =1,e%ndime
!               tract(idime) = dot_product(e%baloc(:,ndime),StressTensor(:,idime))
!           end do
! 
!           !Now we fill traction vector
!           do inodb = 1,e%pnodb
!              inode = e%lboel(inodb)
!              a%btraction(:,inodb) = a%btraction(:,inodb) + tract(:)*e%shapb(inodb,e%igaub)*dsurf
!           end do
! 
!       end do GaussPointLoop 
! 
!   end do boundaries 
! 
!   call a%Mesh%ArrayCommunicator%GhostCommunicate(e%ndime,a%btraction)    
! 
!   !Dellocation
!   call a%Memor%dealloc(ndime,tract,'tract','plcd_EndBouope')     
!   call a%Memor%dealloc(e%ndime,e%ndime,StressTensor,'StressTensor','plcd_EndBouope')
! 
!   !Element Deallocation
!   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','plcd_EndBouope')
! 
!    endif
end subroutine
