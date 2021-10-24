module Mod_sup_elmdir 
   
contains
   
   subroutine sup_Exaelmdir(a,e,ndofn,ndofbc,ndofbcstart,currentbvess,elmat,elrhs)
      !------------------------------------------------------------------------
      ! This routine modifies the element stiffness to impose the correct 
      ! boundary conditions
      !------------------------------------------------------------------------
      use typre
      use Mod_Mesh
      use Mod_SupExacso       
      use Mod_ThreeField
      use Mod_Element
      implicit none
      class(ThreeFieldNSProblem)              :: a
      class(FiniteElement), intent(in)        :: e
      type(SupExacso)  :: exacso         
      
      integer(ip)                            :: ndofn,ndofbcstart,ndofbc
      real(rp), intent(inout)                :: elmat(ndofn,e%mnode,ndofn,e%mnode),elrhs(ndofn,e%mnode)
      
      integer(ip), intent(in)    :: currentbvess
      real(rp)                   :: adiag
      integer(ip)                :: inode,ipoin,idofn
      integer(ip)                :: iffix_aux
         
         
      !Exact Values
      real(rp)                   :: exvel(e%ndime),exveg(e%ndime,e%ndime),exprg(e%ndime), &
                                    exsig((e%ndime-1)*(e%ndime-1)+2),exsigr((e%ndime-1)*(e%ndime-1)+2,e%ndime)
      real(rp)                   :: expre,uex,vex,pi,wex
      real(rp), pointer          :: coord(:)
      
      
      !Dirichlet boundary conditions
      if(a%kfl_exacs/=0)then
      
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            
            call a%Mesh%GetPointCoord(ipoin,coord)
            
            ! Exact solution
            call exacso%sup_ComputeSolution(e%ndime,coord,a%ctime,a%LogFormulation,a)
            call exacso%sup_GetPressure(e%ndime,expre,exprg)         
            call exacso%sup_GetVelocity(e%ndime,exvel,exveg) 
            if (a%LogFormulation/=0) then
               call exacso%sup_GetPsi(e%ndime,exsig,exsigr)
            else   
               call exacso%sup_GetStress(e%ndime,a%LogFormulation,exsig,exsigr)
            end if      
            
            
            do idofn=1,ndofbc
               iffix_aux = a%kfl_fixno(idofn+currentbvess-1,ipoin)
               if(iffix_aux==1) then
                  
                  a%bvess(1,ipoin,1) = exsig(1)
                  a%bvess(2,ipoin,1) = exsig(2)
                  a%bvess(3,ipoin,1) = exsig(3)
                  a%bvess(4,ipoin,1) = exvel(1)
                  a%bvess(5,ipoin,1) = exvel(2)
               
                  elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                        - elmat(:,1:e%pnode,idofn+ndofbcstart,inode)*a%bvess(currentbvess+idofn-1,ipoin,1)
                  adiag=elmat(idofn+ndofbcstart,inode,idofn+ndofbcstart,inode)
                  elrhs(idofn+ndofbcstart,inode)=adiag*a%bvess(currentbvess+idofn-1,ipoin,1)
                  elmat(ndofbcstart+idofn,inode,:,1:e%pnode)=0.0_rp
                  elmat(:,1:e%pnode,ndofbcstart+idofn,inode)=0.0_rp
                  elmat(ndofbcstart+idofn,inode,ndofbcstart+idofn,inode)=adiag
               end if
            end do
            
            
         end do
         
      elseif(a%kfl_bc_number>0.and.a%kfl_confi==1)then 
            do inode=1,e%pnode
            ipoin=e%lnods(inode)
            
            call a%Mesh%GetPointCoord(ipoin,coord)
 
           
            do idofn=1,ndofbc
               iffix_aux = a%kfl_fixno(idofn+currentbvess-1,ipoin)
               if(iffix_aux==1) then
                  
                  
                  !boundary value definition
                  pi= 4.*atan(1.)
                  
                  if(e%ndime==2)then
                     if(a%kfl_bc_number==1.and.a%kfl_timei==1)then
                        uex=8.0_rp*(1.0_rp+tanh(8.0_rp*(0.5_rp)))*(coord(1)**2)*((1-coord(1))**2)*coord(2)*sin(pi*a%ctime)
                        vex=0.0_rp
                     elseif(a%kfl_bc_number==2.and.a%kfl_timei==1)then
                        uex=8.0_rp*(1.0_rp+tanh(8.0_rp*(a%ctime - 0.5_rp)))*(coord(1)**2)*((1-coord(1))**2)*coord(2)      
                        vex=0.0_rp
                     elseif(a%kfl_bc_number==2.and.a%kfl_timei==0)then
                        uex=16.0_rp*(coord(1)**2)*((1-coord(1))**2)*coord(2)      
                        vex=0.0_rp                  
                     end if
                     a%bvess(4,ipoin,1) = uex
                     a%bvess(5,ipoin,1) = vex
                  
                  elseif(e%ndime==3)then
                  
                     if(a%kfl_bc_number==1.and.a%kfl_timei==1)then
                        uex=256.0_rp*(coord(2)**2)*((1-coord(2))**2)*(coord(1)**2)*((1-coord(1))**2)*coord(3)      
                        vex=0.0_rp 
                        wex=0.0_rp
                     elseif(a%kfl_bc_number==2.and.a%kfl_timei==1)then
                        uex=256.0_rp*(coord(2)**2)*((1-coord(2))**2)*(coord(1)**2)*((1-coord(1))**2)*coord(3)*sin(pi*a%ctime)    
                        vex=0.0_rp
                        wex=0.0_rp
                     end if
                     a%bvess(7,ipoin,1) = uex
                     a%bvess(8,ipoin,1) = vex                     
                     a%bvess(9,ipoin,1) = wex 
                     
                  end if
               
                  elrhs(:,1:e%pnode)=elrhs(:,1:e%pnode) &
                        - elmat(:,1:e%pnode,idofn+ndofbcstart,inode)*a%bvess(currentbvess+idofn-1,ipoin,1)
                  adiag=elmat(idofn+ndofbcstart,inode,idofn+ndofbcstart,inode)
                  elrhs(idofn+ndofbcstart,inode)=adiag*a%bvess(currentbvess+idofn-1,ipoin,1)
                  elmat(ndofbcstart+idofn,inode,:,1:e%pnode)=0.0_rp
                  elmat(:,1:e%pnode,ndofbcstart+idofn,inode)=0.0_rp
                  elmat(ndofbcstart+idofn,inode,ndofbcstart+idofn,inode)=adiag
               end if
            end do
         end do
      
      end if
         
         
   end subroutine sup_Exaelmdir
   
end module

