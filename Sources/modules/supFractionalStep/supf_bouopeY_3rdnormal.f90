subroutine supf_bouope_3rdnormal(a,iboun,e,elrhs,elmat)
!This subroutine computes the boundary contribution, incompressible Navier-Stokes equations
!-----------------------------------------------------------------------
   use typre
   use Mod_SUPFractionalStep
   use Mod_Element
   implicit none
   class(SUPFractionalStepProblem) :: a
   integer(ip)                :: iboun,bcstart,ndofn
   class(FiniteElement)        :: e   
   real(rp),intent(inout)     :: elrhs(1,e%mnode),elmat(1,e%mnode,1,e%mnode)
   real(rp), external :: funcre

   
   
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: gpcod(e%ndime)   
   real(rp)    :: tract(e%ndime)
   real(rp)    :: updbcn,coorg
   real(rp)    :: bosig((e%ndime-1)*(e%ndime-1)+2,e%pnodb),gpsig((e%ndime-1)*(e%ndime-1)+2)
   real(rp)    :: bosig2((e%ndime-1)*(e%ndime-1)+2,e%pnodb),gpsig2((e%ndime-1)*(e%ndime-1)+2)
   real(rp)    :: grsig((e%ndime-1)*(e%ndime-1)+2,e%ndime),grsig2((e%ndime-1)*(e%ndime-1)+2,e%ndime)
   real(rp)    :: divsig(e%ndime),divsig2(e%ndime),bopre(e%pnodb,2)
   real(rp)    :: gppre(1),gppre2(1),grapre(1,e%ndime),grapre2(1,e%ndime)
   
   
   integer(ip) :: igaub,inodb,inode,idime,auxtens,jnodb,jnode

      
      call e%elmdel
      
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      
      !Old time step stress value
      call e%gatherb(auxtens,bosig,a%sigma(:,:,3))
      !Intermediate stress value
      call e%gatherb(auxtens,bosig2,a%sigma(:,:,2))    
      
      !Gathers
      call e%gatherb(1,bopre(:,1),a%press(:,1)) ! p_n+1
      call e%gatherb(1,bopre(:,2),a%press(:,3)) ! p_n 
      
      

      
      dsurf0 = 0.0_rp
      
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         tract=0.0_rp
   
         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         
         call e%interpb(e%ndime,e%bocod,gpcod)
         
         call e%interpb(auxtens,bosig,gpsig)  
         call e%interpb(auxtens,bosig2,gpsig2) 
         
         call e%interpb(1,bopre(:,1),gppre)  
         call e%interpb(1,bopre(:,2),gppre2) 
         
         
         call e%gradientb(auxtens,bosig,grsig)
         call e%gradientb(auxtens,bosig2,grsig2)   
         call e%gradientb(1,bopre(:,1),grapre)
         call e%gradientb(1,bopre(:,2),grapre2) 
         
         
         if(e%ndime==2)then
            divsig(1) = grsig(1,1) + grsig(3,2)
            divsig2(1)= grsig2(1,1) + grsig2(3,2)
            divsig(2) = grsig(3,1) + grsig(2,2)
            divsig2(2)= grsig2(3,1) + grsig2(2,2)
         elseif(e%ndime==3)then
            divsig(1) = grsig(1,1) + grsig(6,2) + grsig(5,3)
            divsig2(1)= grsig2(1,1) + grsig2(6,2) + grsig2(5,3)
            divsig(2) = grsig(6,1) + grsig(2,2) + grsig(4,3)
            divsig2(2)= grsig2(6,1) + grsig2(2,2) + grsig2(4,3)
            divsig(3) = grsig(5,1) + grsig(4,2) + grsig(3,3)
            divsig2(3)= grsig2(5,1) + grsig2(4,2) + grsig2(3,3)
         endif
         
         !Derivatives at the boundary
         call e%elmderb         

            
            if(e%ndime==2)then
            
               tract(1) = e%baloc(1,e%ndime)*(divsig2(1) - divsig(1))
               tract(2) = e%baloc(2,e%ndime)*(divsig2(2) - divsig(2))
                  
             elseif(e%ndime==3)then

               tract(1) = e%baloc(1,e%ndime)*(divsig2(1) - divsig(1))
               tract(2) = e%baloc(2,e%ndime)*(divsig2(2) - divsig(2))
               tract(3) = e%baloc(3,e%ndime)*(divsig2(3) - divsig(3))               
                  
             end if
             
             
       do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do jnodb=1,e%pnodb
            jnode = e%lboel(jnodb)
            do idime=1,e%ndime
               elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
                  + e%shapb(inodb,e%igaub)*(grapre(idime,e%ndime)*e%baloc(idime,e%ndime))*dsurf
            end do
         end do
      end do
             
         
        do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            do idime=1,e%ndime
                  elrhs(1,inode)=elrhs(1,inode)+dsurf*(grapre2(idime,1)*e%baloc(idime,e%ndime) + &
                     tract(idime)*e%shapb(inodb,e%igaub))
            end do   
              
        end do
      end do       


end subroutine supf_bouope_3rdnormal
