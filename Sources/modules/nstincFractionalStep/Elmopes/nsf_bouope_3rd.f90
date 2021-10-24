subroutine nsf_bouope_3rd(a)
   use typre
   use Mod_Element
   use Mod_NSF_Element
   use Mod_nsm_elmdir
   use Mod_php_elmdir
   implicit none
   class(NSFractionalStepProblem) :: a
   
   class(FiniteElement), pointer :: e => NULL()
   
   integer(ip) :: nboun,inodb,jnodb,inode,jnode,iboun,igaub,idime
   real(rp)    :: acden,acvis,dsurf
   
   real(rp), allocatable      :: elmat(:,:,:,:)
   real(rp), allocatable      :: elrhs(:,:)
   real(rp)                   :: dvol,&
                                 gpvno,&
                                 hclen,&
                                 rdinv,&
                                 penal
   
   real(rp), allocatable      :: bopre(:,:), bovel(:,:,:), gpvelbo(:)
   real(rp), allocatable      :: bopreInc(:,:), bovelInc(:,:), gpvelInc(:) !Veloc and Press Increments
   real(rp) :: gppreInc(1,1)
   
   !todo multy materials
   integer(ip) :: imat=1

   
   !Initializations
   call a%Mesh%GetNboun(nboun)
   
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsm_bouop0')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','nsm_bouop0')
   call a%Memor%alloc(e%mnodb,2,bopre,'bopre','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,e%mnodb,2,bovel,'bovel','nsf_elmope_3rd')
   call a%Memor%alloc(e%mnodb,1,bopreInc,'bopreInc','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,e%mnodb,bovelInc,'bovelInc','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,gpvelInc,'gpvelInc','nsf_elmope_3rd')
   call a%Memor%alloc(e%ndime,gpvelbo,'gpvelbo','nsf_boupe_3rd')
   
   !Physical Parameters
   call a%GetPhysicalParameters(imat,acden,acvis)

   rdinv = 1.0_rp/real(e%ndime)
   
   !Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !Initialize
      elmat=0.0_rp
      elrhs=0.0_rp

      !Gathers
      call e%gatherb(1,bopre(:,1),a%press(:,1)) ! p_n+1
      call e%gatherb(1,bopre(:,2),a%press(:,3)) ! p_n
      call e%gatherb(e%ndime,bovel(:,:,1),a%veloc(:,:,1)) ! v_n+1
      call e%gatherb(e%ndime,bovel(:,:,2),a%veloc(:,:,3)) ! v_n+1
      
      !bopreInc is the difference between pressure at n and pressure at n+1
      bopreInc(:,1) = bopre(:,1) - bopre(:,2)
      bovelInc = bovel(:,:,1) - bovel(:,:,2)
      
      call e%elmdel
      
      !Gauss-Point Loop
      gauss: do igaub=1,e%pgaub
         e%igaub = igaub
         
         !Calculate exterior Normal
         call e%bounor
         
         !Derivatives at the boundary
         call e%elmderb
         
         dsurf=e%weigb(e%igaub)*e%eucta
         
      enddo gauss

      !Boundary conditions
      
      !Dirichlet in weak: Contribution to Rhs, U_int
      if(a%kfl_fixbo(iboun)==16)then
         do jnodb = 1,e%pnodb
            jnode = e%lboel(jnodb)
            do inodb = 1,e%pnodb
               inode = e%lboel(inodb)
               do idime = 1,e%ndime
                  elrhs(idime,inode) = elrhs(idime,inode) &
                  +dot_product(elmat(idime,inode,1:e%ndime,jnode),a%veloc(1:e%ndime,e%lnods(jnode),1))
               enddo
            enddo
         enddo
      endif

      !Boundary conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
      call php_elmdir(a,e,e%ndime,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   end do boundaries

   call a%Memor%dealloc(e%mnodb,2,bopre,'bopre','nsf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,e%mnodb,2,bovel,'bovel','nsf_elmope_3rd')
   call a%Memor%dealloc(e%mnodb,1,bopreInc,'bopreInc','nsf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsm_bouop0')
   call a%Memor%dealloc(e%ndime,e%mnode,elrhs,'elrhs','nsm_bouop0')
   call a%Memor%dealloc(e%ndime,e%mnodb,bovelInc,'bovelInc','nsf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,gpvelInc,'gpvelInc','nsf_elmope_3rd')
   call a%Memor%dealloc(e%ndime,gpvelbo,'gpvelbo','nsf_boupe_3rd')

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   
end subroutine
