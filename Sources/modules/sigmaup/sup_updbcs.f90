subroutine sup_updbcs(a)
!-----------------------------------------------------------------------
!****f* Nstinc/nsi_updbcs
! NAME 
!    nsi_updbcs
! DESCRIPTION
!    This routine updates the velocity boundary conditions:
!    1. Before a time step begins
!    2. Before a global iteration begins
!    3. Before an inner iteration begins
! USED BY
!    nsi_begste
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_ThreeField
   use Mod_Element  
   use Mod_Mesh
   use Mod_memor 
   use Mod_SupExacso     
   implicit none
   class(ThreeFieldNSProblem) :: a
   class(FiniteElement), pointer :: e => NULL()
   type(SupExacso)         :: exacso    
   real(rp), pointer       :: meshve(:,:,:) => NULL()
   
   integer(ip)             :: ibopo,idime,ipoin,ndime,npoin,ntens
   integer(ip)             :: iboun,nboun,ielem,inode,inodb
   !Exact Values
   real(rp), allocatable   :: exvel(:),exsig(:),exsigr(:,:),exveg(:,:),psi(:)
   !Nodal coordinates
   real(rp), pointer       :: coord(:)
   real(rp), pointer       :: exnor(:,:)   
   real (rp)               :: dummr
   integer(ip)             :: dummi
   logical                 :: isALE

   call a%Mesh%GetALE(isALE)
   call a%Mesh%GetMeshVeloc(meshve)
   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNdime(ndime)
      

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsi_updbcs')

   if (isALE) then
      do ipoin=1,npoin
         do idime=1,ndime
            if(a%kfl_fixno(idime,ipoin)==1 .or. a%kfl_fixno(idime,ipoin)==7) then
               if(a%kfl_funno(ipoin)==0) then
                 a%bvess(idime,ipoin,1)=a%bvess(idime,ipoin,2)+meshve(idime,ipoin,1)
               end if
            end if
         end do
      end do 

      do iboun=1,nboun
         call a%Mesh%BoundaryLoad(iboun,e)
         if (a%kfl_fixbo(iboun) == 4) then
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
               a%bvess(:,ipoin,1)=a%bvess(:,ipoin,2)*exnor(:,1)
            end do
         end if
      end do
   end if

   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_updbcs')
   
   if (a%kfl_exacs /= 0) then

      ntens=(ndime-1)*(ndime-1) + 2
      
      ! Allocate exact components
      call a%Memor%alloc(ndime,exvel,'exvel','sup_updbcs')  
      call a%Memor%alloc(ndime,ndime,exveg,'exveg','sup_updbcs')  
      call a%Memor%alloc(ntens,exsig,'exsig','sup_upbcs')
      call a%Memor%alloc(ntens,ndime,exsigr,'exsigr','sup_upbcs')
      
      do ipoin = 1,npoin
         call a%Mesh%GetPointCoord(ipoin,coord)    
         
         call exacso%sup_ComputeSolution(ndime,coord,a%ctime,a%LogFormulation,a)
         call exacso%sup_GetVelocity(ndime,exvel,exveg)        
         if (a%LogFormulation/=0) then
            call exacso%sup_GetPsi(ndime,exsig,exsigr)
         else   
            call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
         end if
                     
         do idime=1,ntens
            a%bvess(idime,ipoin,1)=exsig(idime)                 
            if(a%kfl_fixno(idime,ipoin)==1) then         
               a%sigma(idime,ipoin,1)=exsig(idime)
               a%sigma(idime,ipoin,2)=exsig(idime)
            end if
         end do  
      

         do idime=1,ndime
            a%bvess(ntens + idime,ipoin,1)=exvel(idime)                 
            if(a%kfl_fixno(ntens+idime,ipoin)==1) then         
               a%veloc(idime,ipoin,1)=exvel(idime)
               a%veloc(idime,ipoin,2)=exvel(idime)
            end if
         end do   
         
      end do
      
      !deallocate exact components
      call a%Memor%dealloc(ndime,exvel,'exvel','sup_updbcs')  
      call a%Memor%dealloc(ndime,ndime,exveg,'exveg','sup_updbcs')  
      call a%Memor%dealloc(ntens,exsig,'exsig','sup_upbcs')
      call a%Memor%dealloc(ntens,ndime,exsigr,'exsigr','sup_upbcs')    
   endif
   
end subroutine sup_updbcs
