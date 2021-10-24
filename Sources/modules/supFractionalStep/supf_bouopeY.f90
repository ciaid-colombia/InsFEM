subroutine supf_bouopeY(a,iboun,e,elmat,elrhs,ndofn,bcstart)
!This subroutine computes the boundary contribution, incompressible Navier-Stokes equations
!-----------------------------------------------------------------------
   use typre
   use Mod_SupExacso 
   use Mod_SUPFractionalStep
   use Mod_Element
   implicit none
   class(SUPFractionalStepProblem) :: a
   integer(ip)                :: iboun,bcstart,ndofn
   class(FiniteElement)        :: e
   type(SupExacso)  :: exacso    
   
   real(rp),intent(inout)     :: elmat(ndofn,e%mnode,ndofn,e%mnode)              ! Element matrices
   real(rp),intent(inout)     :: elrhs(ndofn,e%mnode)
!    real(rp),intent(inout)     :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode)              ! Element matrices
!    real(rp),intent(inout)     :: elrhs(a%ndofn,e%mnode)
   
   !to do multy materials
   integer(ip) :: imat=1
   
   interface

   end interface
   
   real(rp), external :: funcre

   real(rp)    :: wmatr(e%ndime+1,e%mnode,e%ndime+1,e%mnode)
   real(rp)    :: wrhsi(e%ndime+1,e%mnode)
   
   
   real(rp)    :: bovel(e%ndime,e%pnodb),elvel(e%ndime,e%pnode)
   real(rp)    :: gpcod(e%ndime),gpvel(e%ndime),gpvno,grvel(e%ndime,e%ndime)
   
   real(rp)    :: exveg(e%ndime,e%ndime),exvel(e%ndime)
   
   real(rp)    :: acden,acvis,vista,dynpr
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: tract(e%ndime)
   real(rp)    :: updbcn,coorg
   real(rp)    :: bosig((e%ndime-1)*(e%ndime-1)+2,e%pnodb),gpsig((e%ndime-1)*(e%ndime-1)+2)
   real(rp)    :: bosig2((e%ndime-1)*(e%ndime-1)+2,e%pnodb),gpsig2((e%ndime-1)*(e%ndime-1)+2)
   
   
   integer(ip) :: igaub,inodb,inode,idime,auxtens


   if(    a%kfl_fixbo(iboun)==2 &   !Pressure
      ) then  

      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)
      
      call e%elmdel
      
      auxtens=(e%ndime-1)*(e%ndime-1)+2
      
      !Old time step stress value
      call e%gatherb(auxtens,bosig,a%sigma(:,:,3))
      !Intermediate stress value
      call e%gatherb(auxtens,bosig2,a%sigma(:,:,2))    

      
      dsurf0 = 0.0_rp
      
      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub
         
         wmatr=0.0_rp
         wrhsi=0.0_rp
         tract=0.0_rp
   
         !Calculate exterior Normal
         call e%bounor
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         
         call e%interpb(e%ndime,e%bocod,gpcod)
         call e%interpb(e%ndime,bovel,gpvel)
         
         call e%interpb(auxtens,bosig,gpsig)  
         call e%interpb(auxtens,bosig2,gpsig2) 
         
         !Velocity norm
         call vecnor(gpvel,e%ndime,gpvno,2)
         
         !Derivatives at the boundary
         call e%elmderb
         
         !Pressure: sig.n=-p n
         if(a%kfl_fixbo(iboun)==2) then
            if(a%kfl_funbo(iboun)/=0) then
               updbcn=funcre(a%funpa(a%kfl_funbo(iboun))%a,&
               a%kfl_funty(a%kfl_funbo(iboun),2),&
               a%kfl_funty(a%kfl_funbo(iboun),1),a%bctime)
            else
               updbcn=1.0_rp
            end if
            tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*updbcn*e%baloc(1:e%ndime,e%ndime)
            
            if(e%ndime==2)then
            
               tract(1) = tract(1) + e%baloc(1,e%ndime)*(gpsig2(1) - gpsig(1)) &
                  + e%baloc(2,e%ndime)*(gpsig2(3) - gpsig(3))
               tract(2) = tract(2) + e%baloc(1,e%ndime)*(gpsig2(3) - gpsig(3)) &
                  + e%baloc(2,e%ndime)*(gpsig2(2) - gpsig(2)) 
                  
             elseif(e%ndime==3)then

               tract(1) = tract(1) + e%baloc(1,e%ndime)*(gpsig2(1) - gpsig(1)) &
                  + e%baloc(2,e%ndime)*(gpsig2(6) - gpsig(6)) + e%baloc(3,e%ndime)*(gpsig2(5) - gpsig(5)) 
               tract(2) = tract(2) + e%baloc(1,e%ndime)*(gpsig2(6) - gpsig(6)) &
                  + e%baloc(2,e%ndime)*(gpsig2(2) - gpsig(2)) + e%baloc(3,e%ndime)*(gpsig2(4) - gpsig(4))              
               tract(3) = tract(3) + e%baloc(1,e%ndime)*(gpsig2(5) - gpsig(5)) &
                  + e%baloc(2,e%ndime)*(gpsig2(4) - gpsig(4)) + e%baloc(3,e%ndime)*(gpsig2(3) - gpsig(3))
                  
                  
             end if

         end if
         
          elmat((bcstart+1):(bcstart+e%ndime+1),:,(bcstart+1):(bcstart+e%ndime+1),:)= dsurf*wmatr + &
               elmat((bcstart+1):(bcstart+e%ndime+1),:,(bcstart+1):(bcstart+e%ndime+1),:)            
        
         do inodb=1,e%pnodb
            inode = e%lboel(inodb)
            do idime=1,e%ndime
                  elrhs(bcstart+idime,inode)=elrhs(bcstart+idime,inode)+dsurf*(wrhsi(idime,inode)+tract(idime)*e%shapb(inodb,e%igaub))
            end do   
              
        end do
      end do
   endif
       


end subroutine supf_bouopeY
