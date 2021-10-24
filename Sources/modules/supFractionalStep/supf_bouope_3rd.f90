subroutine supf_bouope_3rd(a)
   use typre
   use Mod_Memor
   use Mod_SUPFractionalStep
   use Mod_php_SetTimeIntegrator 
   use Mod_TimeIntegrator 
   use Mod_SUPF_Element   
   use Mod_Element
   use Mod_nsm_elmdir
   use Mod_php_elmdir
   implicit none
   class(SUPFractionalStepProblem), pointer :: a
   
   class(FiniteElement), pointer :: e => NULL()
   type(TimeIntegratorDt1) :: Integrator   
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv
   real(rp)    :: acden,acvis,raux   
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp), allocatable :: welmat(:,:,:,:)
   real(rp), allocatable :: welrhs(:,:)
   
   integer(ip) :: iboun,nboun,nelem,ndime,idime,inode,bcstart,auxtens,auxnbcs
   integer(ip) :: currentbvess,auxbcstart,aux1bc,aux2bc
   
   !todo multy materials
   integer(ip) :: imat=1
   
   interface    

       subroutine supf_bouope_3rdnormal(a,iboun,e,elrhs,elmat)
          !This subroutine computes the boundary contribution, incompressible Navier-Stokes equations in Three Field Case for Stress
          !-----------------------------------------------------------------------
             use typre
!             use Mod_NavierStokes          
             use Mod_Element
             use Mod_SUPFractionalStep                
             implicit none
             class(SUPFractionalStepProblem) :: a
             integer(ip)                :: iboun
             class(FiniteElement)        :: e
             real(rp),intent(inout)     :: elrhs(1,e%mnode),elmat(1,e%mnode,1,e%mnode)
       end subroutine         
       
       
    end interface
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   !Set Time Integrator
   call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)   
   !Physical Parameters
   call a%GetPhysicalParameters(imat,acden,acvis)   
   
   auxnbcs=a%ndofbc+a%ndofn
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_bouope_3rd')

   auxtens=(e%ndime-1)*(e%ndime-1)+2
   
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmat,'elmat','supf_bouope_1st')
   call a%Memor%alloc(1,e%mnode,elrhs,'elrhs','supf_bouope_1st')
   !Working arrays
   call a%Memor%alloc(1,e%mnode,welrhs,'welrhs','supf_bouope_1st')
   call a%Memor%alloc(1,e%mnode,1,e%mnode,welmat,'welmat','supf_bouope_1st')   
  
   write(*,*) 'supf_bouope_3rd: at this point, no statistics about wall law'
   
   !Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !Initialize
      welrhs=0.0_rp
      welmat=0.0_rp

      !Compute boundary terms
      call supf_bouope_3rdnormal(a,iboun,e,welrhs,welmat)    
      
      !At this point we do nothing with the pressure terms (either equation and unknown)

      raux = 1.0_rp/(acden*LHSdtinv)
      
      elmat(1,1:e%pnode,1,1:e%pnode) = raux*welmat(1,1:e%pnode,1,1:e%pnode)
      
      elrhs(1,1:e%pnode) = raux*welrhs(1,1:e%pnode) 

      !Dirichlet Boundary conditions
      if(a%kfl_exacs/=0.and.a%kfl_timei==1)then
         call supf_ExaPreElmdir(a,e,elmat,elrhs)
      else
         call supf_elmdir_press(a,e,elmat,elrhs)
      end if
      !Assembly
      call a%LinearSystemC%Assembly(e,elmat,elrhs)
      
   end do boundaries

   call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmat,'elmat','supf_bouope_1st')
   call a%Memor%dealloc(1,e%mnode,elrhs,'elrhs','supf_bouope_1st')
   call a%Memor%dealloc(1,e%mnode,welrhs,'welrhs','supf_bouope_1st')
   call a%Memor%dealloc(1,e%mnode,1,e%mnode,welmat,'welmat','supf_bouope_1st')     
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supf_bouope_3rd')


end subroutine supf_bouope_3rd
