subroutine supf_bouope_4rd(a)
   use typre
   use Mod_Element
   use Mod_SUPFractionalStep
   use Mod_Element
   use Mod_SUPF_Element
   use Mod_nsm_elmdir
   use Mod_php_elmdir   
   implicit none
   class(SUPFractionalStepProblem) :: a
   
   class(FiniteElement), pointer :: e => NULL()
   
   integer(ip) :: nelem,nboun,inodb,jnodb,inode,jnode,iboun,igaub,idime
   real(rp)    :: acden,acvis,dsurf
   integer(ip)           :: nsteps
   real(rp)              :: LHSdtinv
   
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp), allocatable :: elrhs(:,:)
   real(rp)              :: dvol,gppre(1)
   
   real(rp), allocatable :: bosig(:,:,:),gpsig(:),tract(:)
   real(rp), allocatable :: bopre(:,:)   
   integer(ip) :: currentbvess,bcstar,auxtens
   
   !todo multy materials
   integer(ip) :: imat=1
   

   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   
   !Memory allocation
  call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supf_bouope_3rd')
   
   auxtens=(e%ndime-1)*(e%ndime-1)+2
   
   call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsm_bouop0')
   call a%Memor%alloc(e%ndime,e%mnode,elrhs,'elrhs','nsm_bouop0')
   call a%Memor%alloc(auxtens,e%mnodb,1,bosig,'bopre','nsf_elmope_3rd')
   call a%Memor%alloc(auxtens,gpsig,'gpsig','nsf_elmope_3rd')   
   call a%Memor%alloc(e%ndime,tract,'tract','nsf_elmope_3rd')
   call a%Memor%alloc(e%mnodb,1,bopre,'bopre','nsf_elmope_3rd')   
   
   !Physical Parameters
   call a%GetPhysicalParameters(imat,acden,acvis)
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      !Initialize
      elmat=0.0_rp
      elrhs=0.0_rp

      !Gathers
      call e%gatherb(auxtens,bosig,a%sigma(:,:,3))
      call e%gatherb(1,bopre,a%press(:,3)) ! p_n      

      
      call e%elmdel
      
      
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         !Calculate exterior Normal
         call e%bounor
         
         !Derivatives at the boundary
         call e%elmderb
         
         dsurf=e%weigb(e%igaub)*e%eucta
         
         !Open flow, integration by parts of (v,grad p) => (n*v,p)_Gamma
         if (a%kfl_fixbo(iboun) == 14) then
            
            call e%interpb(auxtens,gpsig,bosig)
            call e%interpb(1,bopre,gppre)
            
            if(e%ndime==2)then
               tract(1) = e%baloc(1,e%ndime)*gpsig(1) + e%baloc(2,e%ndime)*gpsig(3)
               tract(2) = e%baloc(1,e%ndime)*gpsig(3) + e%baloc(2,e%ndime)*gpsig(2)
            elseif(e%ndime==3)then
               tract(1) = e%baloc(1,e%ndime)*gpsig(1) + e%baloc(2,e%ndime)*gpsig(6) + e%baloc(3,e%ndime)*gpsig(5) 
               tract(2) = e%baloc(1,e%ndime)*gpsig(6) + e%baloc(2,e%ndime)*gpsig(2) + e%baloc(3,e%ndime)*gpsig(4) 
               tract(3) = e%baloc(1,e%ndime)*gpsig(5) + e%baloc(2,e%ndime)*gpsig(4) + e%baloc(3,e%ndime)*gpsig(3)
            end if
         
         else if (a%kfl_fixbo(iboun) == 15)then
         
             if(e%ndime==2)then
               tract(1) = -e%baloc(1,e%ndime)*gppre(1)
               tract(2) = -e%baloc(2,e%ndime)*gppre(1)
            elseif(e%ndime==3)then
               tract(1) = -e%baloc(1,e%ndime)*gppre(1)
               tract(2) = -e%baloc(2,e%ndime)*gppre(1) 
               tract(3) = -e%baloc(3,e%ndime)*gppre(1)
            end if        
         
         
         endif
         
         do inodb = 1,e%pnodb
            inode = e%lboel(inodb)
            do idime = 1,e%ndime
               elrhs(idime,inode)=elrhs(idime,inode)+dsurf*(tract(idime)*e%shapb(inodb,e%igaub))
            enddo
         enddo
         
      enddo

      !Boundary conditions
      call nsm_rotdir(a,e,e%ndime,elmat,elrhs)
       ! php_elmdir(a,e,ndofn,ndofbc,ndofbcstart,currentbvess,elmat,elrhs)
      currentbvess=auxtens+1
      bcstar=0_ip
      call php_elmdir(a,e,e%ndime,e%ndime,bcstar,currentbvess,elmat,elrhs)
      
      ! Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
      
   end do boundaries

   call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmat,'elmat','nsm_bouop0')
   call a%Memor%dealloc(e%ndime,e%mnode,elrhs,'elrhs','nsm_bouop0')
   call a%Memor%dealloc(auxtens,e%mnodb,1,bosig,'bopre','nsf_elmope_3rd')
   call a%Memor%dealloc(auxtens,gpsig,'gpsig','nsf_elmope_3rd')  
   call a%Memor%dealloc(e%ndime,tract,'tract','nsf_elmope_3rd')  
   call a%Memor%dealloc(e%mnodb,1,bopre,'bopre','nsf_elmope_3rd')   
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   
   
end subroutine