module Mod_supf_bouope_2nd
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_BoundaryConditionsHooksAndSubroutines

  implicit none
contains
   subroutine SetBCPointers
      implicit none
      call ResetProcedureComposition
      
      !Set All Pointers To NULLSUB (in Mod_nsm_BaseElmope)
      call SetHooksToNULL
      call SetPointersBoundaryConditions(0)
      call SetPointersBoundaryConditions(1)
      call SetPointersBoundaryConditions(100)
   end subroutine
end module


subroutine supf_bouope_2nd(b)
   use typre
   use Mod_Memor
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_nsm_elmdir
   use Mod_php_elmdir
   use Mod_supf_bouope_2nd
   implicit none
   class(SUPFractionalStepProblem), pointer :: b
   
   real(rp), allocatable :: welmat(:,:,:,:)
   real(rp), allocatable :: welrhs(:,:)
   integer(ip) :: nboun,ndime
   
   integer(ip) :: currentbvess,auxbcstart,aux1bc,aux2bc   
   
   interface
       
       subroutine sup_bouope(a,iboun,e,elmat,elrhs,ndofn,bcstart)
          !This subroutine computes the boundary contribution, incompressible Navier-Stokes equations in Three Field Case for Stress
          !-----------------------------------------------------------------------
             use typre        
             use Mod_Element
             use Mod_ThreeField                
             implicit none
             class(ThreeFieldNSProblem) :: a
             integer(ip)                :: iboun
             class(FiniteElement)        :: e
             integer(ip), intent(in)    :: bcstart,ndofn
             real(rp),intent(inout)     :: elmat(ndofn,e%mnode,ndofn,e%mnode)              ! Element matrices
             real(rp),intent(inout)     :: elrhs(ndofn,e%mnode)
       end subroutine 
    end interface
    
    a=>b
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)  
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_bouope_2nd')
   
   auxtens=(e%ndime-1)*(e%ndime-1)+2      
   
   call a%Memor%alloc(auxtens,e%mnode,auxtens,e%mnode,elmat,'elmat','supf_bouope_2nd')
   call a%Memor%alloc(auxtens,e%mnode,elrhs,'elrhs','supf_bouope_2nd')
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,welmat,'welmat','supf_bouope_2nd')
   call a%Memor%alloc(a%ndofn,e%mnode,welrhs,'welrhs','supf_bouope_2nd')
   
   
   write(*,*) 'supf_bouope_2nd: at this point, no statistics about wall law'
   
   !Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !Initialize
      welmat=0.0_rp
      welrhs=0.0_rp
      bcstart=(e%ndime-1)*(e%ndime-1)+2
      !Compute boundary terms
      call sup_bouope(a,iboun,e,welmat,welrhs,a%ndofn,1_ip)
      
      !At this point we do nothing with the pressure terms (either equation and unknown)
      elmat(1:auxtens,1:e%pnode,1:auxtens,1:e%pnode) = welmat(1:auxtens,1:e%pnode,1:auxtens,1:e%pnode)
      elrhs(1:auxtens,1:e%pnode) = welrhs(1:auxtens,1:e%pnode) 
      
      do inode = 1,e%pnode
         do idime = 1,auxtens
            elrhs(idime,inode) = elrhs(idime,inode) - dot_product(welmat(idime,inode,e%ndime+1+auxtens,1:e%pnode),a%press(e%lnods(1:e%pnode),1))
         enddo
      enddo
      
      !Boundary conditions
      call nsm_rotdir(a,e,auxtens,elmat,elrhs)
      currentbvess=1
      auxbcstart=0_ip
      call php_elmdir(a,e,auxtens,auxtens,auxbcstart,currentbvess,elmat,elrhs)

      !Renumeration in case of periodic bc's
      !do inode=1,pnode
      !   lnodf(inode)=fictio(lnods(ispos+inode))
      !end do

      !Assembly
      call b%LinearSystemS%Assembly(e,elmat,elrhs)
      
   end do boundaries

   call a%Memor%dealloc(auxtens,e%mnode,auxtens,e%mnode,elmat,'elmat','supf_bouope_2nd')
   call a%Memor%dealloc(auxtens,e%mnode,elrhs,'elrhs','supf_bouope_2nd')
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,welmat,'welmat','supf_bouope_2nd')
   call a%Memor%dealloc(a%ndofn,e%mnode,welrhs,'welrhs','supf_bouope_2nd')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supf_bouope_2nd')
   
end subroutine supf_bouope_2nd
