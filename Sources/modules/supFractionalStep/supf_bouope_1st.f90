subroutine supf_bouope_1st(a)
   use typre
   use Mod_Memor
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_sup_elmdir
   use Mod_nsm_elmdir
   use Mod_php_elmdir
   implicit none
   class(SUPFractionalStepProblem), pointer :: a
   
   class(FiniteElement), pointer :: e => NULL()
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp), allocatable :: welmat(:,:,:,:)
   real(rp), allocatable :: welrhs(:,:)
   
   integer(ip) :: iboun,nboun,nelem,ndime,idime,inode,bcstart,auxtens,auxnbcs
   integer(ip) :: currentbvess,auxbcstart,aux1bc,aux2bc
   
   interface

       subroutine supf_bouope(a,iboun,e,elmat,elrhs,ndofn,bcstart)
          !This subroutine computes the boundary contribution, incompressible Navier-Stokes equations in Three Field Case for Stress
          !-----------------------------------------------------------------------
             use typre    
             use Mod_Element
             use Mod_SUPFractionalStep                
             implicit none
             class(SUPFractionalStepProblem) :: a
             integer(ip)                :: iboun
             integer(ip), intent(in)    :: bcstart,ndofn             
             class(FiniteElement)        :: e
             real(rp),intent(inout)     :: elmat(ndofn,e%mnode,ndofn,e%mnode)              ! Element matrices
             real(rp),intent(inout)     :: elrhs(ndofn,e%mnode)
       end subroutine         
       
       
    end interface
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   auxnbcs=a%ndofbc+a%ndofn
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sup_bouope')

   auxtens=(e%ndime-1)*(e%ndime-1)+2
   
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','supf_bouope_1st')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','supf_bouope_1st')
   !Working arrays
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,welmat,'welmat','supf_bouope_1st')
   call a%Memor%alloc(a%ndofn,e%mnode,welrhs,'welrhs','supf_bouope_1st')
   
  
   write(*,*) 'supf_bouope_1st: at this point, no statistics about wall law'
   
   !Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !Initialize
      welmat=0.0_rp
      welrhs=0.0_rp
      bcstart=(e%ndime-1)*(e%ndime-1)+2
      !Compute boundary terms
      call supf_bouope(a,iboun,e,welmat,welrhs,a%ndofn,bcstart)  
      
      !At this point we do nothing with the pressure terms (either equation and unknown)
      aux1bc=1+auxtens
      aux2bc=e%ndime+auxtens
      elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = welmat(aux1bc:aux2bc,1:e%pnode,aux1bc:aux2bc,1:e%pnode)
      elrhs(1:e%ndime,1:e%pnode) = welrhs(aux1bc:aux2bc,1:e%pnode) 
      do inode = 1,e%pnode
         do idime = 1,e%ndime
            elrhs(idime,inode) = elrhs(idime,inode) - dot_product(welmat(auxtens+idime,inode,e%ndime+1+auxtens,1:e%pnode),a%press(e%lnods(1:e%pnode),1))
         enddo
      enddo
      
      !Boundary conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
      currentbvess=auxtens+1
      auxbcstart=0_ip
      if(a%kfl_exacs/=0.and.a%kfl_timei==1)then
         call sup_Exaelmdir(a,e,e%ndime,e%ndime,auxbcstart,currentbvess,elmat,elrhs)
      elseif(a%kfl_bc_number>0.and.a%kfl_confi==1)then
         call sup_Exaelmdir(a,e,e%ndime,e%ndime,auxbcstart,currentbvess,elmat,elrhs)
      else
         call php_elmdir(a,e,e%ndime,e%ndime,auxbcstart,currentbvess,elmat,elrhs)
      end if 
            
      

      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   end do boundaries

   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','supf_bouope_1st')
   call a%Memor%dealloc(e%ndime,e%mnode,elrhs,'elrhs','supf_bouope_1st')
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,welmat,'welmat','supf_bouope_1st')
   call a%Memor%dealloc(a%ndofn,e%mnode,welrhs,'welrhs','supf_bouope_1st')   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sup_bouope')
end subroutine supf_bouope_1st
