module Mod_SUPF_Element
   use typre
   use Mod_SUPFractionalStep
   use Mod_Element
   implicit none

contains
   
   !---------------------------------------------------------------------
   !First Step "First Approach" Split Momentum equation
   !---------------------------------------------------------------------
   
   subroutine supf_elmbuv1(e,dvolu,acvis,beta,acden,dtinv,gpadv,timom,gpvel,grvel,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpadv(e%pnode),gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,acden,dtinv,timom,acvis,beta
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime
      real(rp)                   :: aux,aux1,aux2

      aux  = dvolu*acden
      aux2 = timom*acden*acden*dvolu
      
      do inode=1,e%pnode
         do jnode=1,e%pnode         

            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
               + aux*e%shape(inode,e%igaus)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(1,1)) &
               + aux2*gpadv(inode)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(1,1))            
            elmat(1,inode,2,jnode) = elmat(1,inode,2,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(1,2))
               
            elmat(2,inode,1,jnode) = elmat(2,inode,1,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(2,1))               
            elmat(2,inode,2,jnode) = elmat(2,inode,2,jnode) &
               + aux*e%shape(inode,e%igaus)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(2,2)) &
               + aux2*gpadv(inode)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(2,2)) 
                 
               
               
         end do
      end do 

   end subroutine supf_elmbuv1 
   
   
   subroutine supf_elmbuv1Y(e,dvolu,acvis,beta,acden,dtinv,gpadv,timom,gpvel,grvel,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpadv(e%pnode),gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,acden,dtinv,timom,acvis,beta
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime
      real(rp)                   :: aux,aux1,aux2,gpconv(e%pnode),gpconv2(e%pnode)

      aux  = dvolu*acden
      aux2 = timom*acden*acden*dvolu
      
      do inode=1,e%pnode
      
         gpconv(inode)= (gpvel(1)*e%cartd(1,inode) + gpvel(2)*e%cartd(1,inode))
         
         do jnode=1,e%pnode   
         
            gpconv2(jnode)= (gpvel(1)*e%cartd(1,jnode) + gpvel(2)*e%cartd(1,jnode))
         
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
               + aux*e%shape(inode,e%igaus)*(dtinv*e%shape(jnode,e%igaus) + gpadv(jnode)) &
               + aux2*gpadv(inode)*(gpadv(jnode))            
            elmat(1,inode,2,jnode) = elmat(1,inode,2,jnode) + 0.0_rp               
            
            elmat(2,inode,1,jnode) = elmat(2,inode,1,jnode) + 0.0_rp            
            elmat(2,inode,2,jnode) = elmat(2,inode,2,jnode) &
               + aux*e%shape(inode,e%igaus)*(dtinv*e%shape(jnode,e%igaus) + gpadv(jnode)) &
               + aux2*gpadv(inode)*(gpadv(jnode)) 

           
         end do
      end do 

   end subroutine supf_elmbuv1Y    
   
   subroutine supf_elmbuv13d(e,dvolu,acvis,beta,acden,dtinv,gpadv,timom,gpvel,grvel,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpadv(e%pnode),gpvel(e%ndime),grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,acden,dtinv,timom,acvis,beta
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime
      real(rp)                   :: aux,aux1,aux2
      
      aux  = dvolu*acden
      aux2 = timom*acden*acden*dvolu
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
         
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
               + aux*e%shape(inode,e%igaus)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(1,1)) &
               + aux2*gpadv(inode)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(1,1))            
            elmat(1,inode,2,jnode) = elmat(1,inode,2,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(1,2))
            elmat(1,inode,3,jnode) = elmat(1,inode,3,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(1,3))               
               
            elmat(2,inode,1,jnode) = elmat(2,inode,1,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(2,1))               
            elmat(2,inode,2,jnode) = elmat(2,inode,2,jnode) &
               + aux*e%shape(inode,e%igaus)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(2,2)) &
               + aux2*gpadv(inode)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(2,2)) 
            elmat(2,inode,3,jnode) = elmat(2,inode,3,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(2,3))  
               
            elmat(3,inode,1,jnode) = elmat(3,inode,1,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(3,1))                
            elmat(3,inode,2,jnode) = elmat(3,inode,2,jnode) &
               + (aux*e%shape(inode,e%igaus) + aux2*gpadv(inode))*(e%shape(jnode,e%igaus)*grvel(3,2))
            elmat(3,inode,3,jnode) = elmat(3,inode,3,jnode) &
               + aux*e%shape(inode,e%igaus)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(3,3)) &
               + aux2*gpadv(inode)*(gpadv(jnode) + e%shape(jnode,e%igaus)*grvel(3,3))               
               
         end do
      end do 

   end subroutine supf_elmbuv13d    
   
   subroutine supf_elmbuv_estS(e,dvolu,acvis,acden,dtinv,tidiv,tisig,beta,grvel,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,tisig,beta,acvis,acden,dtinv,tidiv
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime
      real(rp)                   :: aux,aux1,aux2    
      

      aux = (2.0_rp*beta*acvis + (1.0_rp-beta)*tisig)*dvolu       
      aux1= acden*dvolu*dtinv
      aux2=tidiv*dvolu
        
      do inode=1,e%pnode
         do jnode=1,e%pnode
         
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + aux1*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) &
               + e%cartd(1,inode)*e%cartd(1,jnode)*aux2 &            
               + aux*(e%cartd(1,inode)*e%cartd(1,jnode) + 0.5_rp*e%cartd(2,inode)*e%cartd(2,jnode))
            elmat(1,inode,2,jnode) = elmat(1,inode,2,jnode) &
               + e%cartd(1,inode)*e%cartd(2,jnode)*aux2 &            
               + aux*(0.5_rp*e%cartd(2,inode)*e%cartd(1,jnode))
               
            elmat(2,inode,1,jnode) = elmat(2,inode,1,jnode) &
               + e%cartd(2,inode)*e%cartd(1,jnode)*aux2 &            
               + aux*(0.5_rp*e%cartd(1,inode)*e%cartd(2,jnode))               
            elmat(2,inode,2,jnode) = elmat(2,inode,2,jnode) + aux1*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) &
               + e%cartd(2,inode)*e%cartd(2,jnode)*aux2 &            
               + aux*(e%cartd(2,inode)*e%cartd(2,jnode) + 0.5_rp*e%cartd(1,inode)*e%cartd(1,jnode))
               
                
         end do
      end do 


   end subroutine supf_elmbuv_estS
   
   subroutine supf_elmbuv_estS3d(e,dvolu,acvis,acden,dtinv,tidiv,tisig,beta,grvel,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the lhs terms for ASGS for block U,V 
      !   (v, a·grad u) + s*(v,u) - tau1*(L*v, Lu) + (v, u/dt) - tau1*(L*v, u/dt)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,tisig,beta,acvis,acden,dtinv,tidiv
      real(rp),    intent(inout) :: elmat(e%ndime,e%mnode,e%ndime,e%mnode)
      integer(ip)                :: inode,jnode,idime,jdime
      real(rp)                   :: aux,aux1,aux2
      

      aux = (2.0_rp*beta*acvis + (1.0_rp-beta)*tisig)*dvolu      
      aux1=acden*dvolu*dtinv
      aux2=dvolu*tidiv
      
      do inode=1,e%pnode
         do jnode=1,e%pnode
         
            elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + aux1*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) &
               + e%cartd(1,inode)*e%cartd(1,jnode)*aux2 &
               + aux*(e%cartd(1,inode)*e%cartd(1,jnode)+0.5_rp*e%cartd(2,inode)*e%cartd(2,jnode)+0.5_rp*e%cartd(3,inode)*e%cartd(3,jnode))
            elmat(1,inode,2,jnode) = elmat(1,inode,2,jnode) &
               + e%cartd(1,inode)*e%cartd(2,jnode)*aux2 &
               + aux*(0.5_rp*e%cartd(2,inode)*e%cartd(1,jnode))
            elmat(1,inode,3,jnode) = elmat(1,inode,3,jnode) &
               + e%cartd(1,inode)*e%cartd(3,jnode)*aux2 &            
               + aux*(0.5_rp*e%cartd(3,inode)*e%cartd(1,jnode))                 
               
            elmat(2,inode,1,jnode) = elmat(2,inode,1,jnode) &
               + e%cartd(2,inode)*e%cartd(1,jnode)*aux2 &              
               + aux*(0.5_rp*e%cartd(1,inode)*e%cartd(2,jnode))               
            elmat(2,inode,2,jnode) = elmat(2,inode,2,jnode) + aux1*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) &
               + e%cartd(2,inode)*e%cartd(2,jnode)*aux2 &              
               + aux*(e%cartd(2,inode)*e%cartd(2,jnode)+0.5_rp*e%cartd(1,inode)*e%cartd(1,jnode)+0.5_rp*e%cartd(3,inode)*e%cartd(3,jnode)) 
            elmat(2,inode,3,jnode) = elmat(2,inode,3,jnode) &
               + e%cartd(2,inode)*e%cartd(3,jnode)*aux2 &              
               + aux*(0.5_rp*e%cartd(3,inode)*e%cartd(2,jnode))  
               
            elmat(3,inode,1,jnode) = elmat(3,inode,1,jnode) &
               + e%cartd(3,inode)*e%cartd(1,jnode)*aux2 &              
               + aux*(0.5_rp*e%cartd(1,inode)*e%cartd(3,jnode))                  
            elmat(3,inode,2,jnode) = elmat(3,inode,2,jnode) &
               + e%cartd(3,inode)*e%cartd(2,jnode)*aux2 &              
               + aux*(0.5_rp*e%cartd(2,inode)*e%cartd(3,jnode))
            elmat(3,inode,3,jnode) = elmat(3,inode,3,jnode) + aux1*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus) &
               + e%cartd(3,inode)*e%cartd(3,jnode)*aux2 &  
               + aux*(e%cartd(3,inode)*e%cartd(3,jnode)+0.5_rp*e%cartd(1,inode)*e%cartd(1,jnode)+0.5_rp*e%cartd(2,inode)*e%cartd(2,jnode))                
               
               
         end do
      end do 


   end subroutine supf_elmbuv_estS3d    
   
   subroutine supf_elmrhu1(e,acden,timom,dvolu,elext,gpvel,grvel,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elext(e%ndime),acden,gpadv(e%pnode)
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime),timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2


      do inode=1,e%pnode
         
         aux = e%shape(inode,e%igaus)*dvolu   
         
         elrhs(1,inode)= aux*(elext(1)) + elrhs(1,inode)            
         elrhs(2,inode)= aux*(elext(2)) + elrhs(2,inode)
          
       end do


   end subroutine supf_elmrhu1  
   
   
   subroutine supf_elmrhuConv(e,acden,timom,dvolu,elext,gpvel,grvel,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elext(e%ndime),acden,gpadv(e%pnode)
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime),timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2


      do inode=1,e%pnode
         
         aux = e%shape(inode,e%igaus)*dvolu       
         aux2= (timom*acden*gpadv(inode))*dvolu
         
         elrhs(1,inode)= aux*(acden*(gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2))) &
            + aux2*acden*(gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2)) + elrhs(1,inode)
            
         elrhs(2,inode)= aux*(acden*(gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2))) &
            + aux2*acden*(gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2)) + elrhs(2,inode)  
       end do


   end subroutine supf_elmrhuConv    
   
   
   
   subroutine supf_elmrhu1Y(e,acden,timom,dvolu,elext,gpvel,grvel,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elext(e%ndime),acden,gpadv(e%pnode)
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime),timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2


      do inode=1,e%pnode
         
         aux = e%shape(inode,e%igaus)*dvolu      
         
         elrhs(1,inode)= aux*(elext(1)) + elrhs(1,inode)
            
         elrhs(2,inode)= aux*(elext(2)) + elrhs(2,inode)     
            
            
       end do


   end subroutine supf_elmrhu1Y  
       
   
   subroutine supf_elmrhu13d(e,acden,timom,dvolu,elext,gpvel,grvel,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elext(e%ndime),acden,gpadv(e%pnode)
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime),timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2


      do inode=1,e%pnode
         
         aux = e%shape(inode,e%igaus)*dvolu       
         
         elrhs(1,inode)= aux*elext(1) + elrhs(1,inode)
            
         elrhs(2,inode)= aux*elext(2) + elrhs(2,inode)
            
         elrhs(3,inode)= aux*elext(3) + elrhs(3,inode)
            
            
       end do

   end subroutine supf_elmrhu13d    
   
   
   subroutine supf_elmrhuConv3d(e,acden,timom,dvolu,elext,gpvel,grvel,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: elext(e%ndime),acden,gpadv(e%pnode)
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),grvel(e%ndime,e%ndime),timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2


      do inode=1,e%pnode
         
         aux = e%shape(inode,e%igaus)*dvolu       
         aux2= (timom*acden*gpadv(inode))*dvolu
         
         elrhs(1,inode)= aux*(acden*(gpvel(1)*grvel(1,1)+gpvel(2)*grvel(1,2)+gpvel(3)*grvel(1,3))) &
            + aux2*acden*(gpvel(1)*grvel(1,1)+gpvel(2)*grvel(1,2)+gpvel(3)*grvel(1,3)) + elrhs(1,inode)
            
         elrhs(2,inode)= aux*(acden*(gpvel(1)*grvel(2,1)+gpvel(2)*grvel(2,2)+gpvel(3)*grvel(2,3))) &
            + aux2*acden*(gpvel(1)*grvel(2,1)+gpvel(2)*grvel(2,2)+gpvel(3)*grvel(2,3)) + elrhs(2,inode)
            
         elrhs(3,inode)= aux*(acden*(gpvel(1)*grvel(3,1)+gpvel(2)*grvel(3,2)+gpvel(3)*grvel(3,3))) &
            + aux2*acden*(gpvel(1)*grvel(3,1)+gpvel(2)*grvel(3,2)+gpvel(3)*grvel(3,3)) + elrhs(3,inode)
            
            
       end do

   end subroutine supf_elmrhuConv3d    
   
   
   subroutine supf_rhu_soss_step1(e,dvolu,acden,tidiv,gpdiv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: acden
      real(rp),    intent(in)    :: dvolu,tidiv,gpdiv(1)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux1,aux2
      
      
      aux1=tidiv*dvolu
      
      do inode=1,e%pnode      

         do idime=1,e%ndime

             elrhs(idime,inode) = aux1*gpdiv(1)*e%cartd(idime,inode) + elrhs(idime,inode)
               
         end do
       end do      

   end subroutine supf_rhu_soss_step1  
   
   
   subroutine supf_rhu_soss_conv(e,dvolu,acden,timom,gpadv,gpconv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: acden,gpconv(e%ndime),gpadv(e%pnode)
      real(rp),    intent(in)    :: dvolu,timom
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux2     

      
      do inode=1,e%pnode              
      
         aux2= (timom*acden*gpadv(inode))*dvolu      
         
         do idime=1,e%ndime

             elrhs(idime,inode) = aux2*gpconv(idime) + elrhs(idime,inode)
               
         end do
       end do      

   end subroutine supf_rhu_soss_conv   
   
   
   subroutine supf_elmrhuextra(e,dvolu,auxtens,gpsig,gppre,elrhs)
      !This subroutine computes the rhs in nsf1_elmope, additional terms involving the pressure and stress
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,gppre(1),gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip) :: inode
 
      do inode=1,e%pnode
      
         elrhs(1,inode) = elrhs(1,inode) + gppre(1)*e%cartd(1,inode)*dvolu &
            - (gpsig(1)*e%cartd(1,inode) + gpsig(3)*e%cartd(2,inode))*dvolu
         elrhs(2,inode) = elrhs(2,inode) + gppre(1)*e%cartd(2,inode)*dvolu &
            - (gpsig(3)*e%cartd(1,inode) + gpsig(2)*e%cartd(2,inode))*dvolu            
            
      end do
   end subroutine   

   subroutine supf_elmrhuextra3d(e,dvolu,auxtens,gpsig,gppre,elrhs)
      !This subroutine computes the rhs in nsf1_elmope, additional terms involving the pressure and stress
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,gppre(1),gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip) :: inode
 
      do inode=1,e%pnode
      
         elrhs(1,inode) = elrhs(1,inode) + gppre(1)*e%cartd(1,inode)*dvolu &
            - (gpsig(1)*e%cartd(1,inode) + gpsig(6)*e%cartd(2,inode) + gpsig(5)*e%cartd(3,inode))*dvolu
         elrhs(2,inode) = elrhs(2,inode) + gppre(1)*e%cartd(2,inode)*dvolu &
            - (gpsig(6)*e%cartd(1,inode) + gpsig(2)*e%cartd(2,inode) + gpsig(4)*e%cartd(3,inode))*dvolu
         elrhs(3,inode) = elrhs(3,inode) + gppre(1)*e%cartd(3,inode)*dvolu &
            - (gpsig(5)*e%cartd(1,inode) + gpsig(4)*e%cartd(2,inode) + gpsig(3)*e%cartd(3,inode))*dvolu             
            
      end do
   end subroutine 
   
   
   subroutine supf_elmrhu_oss(e,auxtens,tisig,dvolu,gprep,beta,elrhs)
      !This subroutine computes the rhs in nsf1_elmope, additional terms involving the pressure and stress
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,gprep(auxtens+e%ndime+1),tisig,beta
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      real(rp)    :: aux1
      integer(ip) :: inode
 
      aux1=-tisig*dvolu
      
      do inode=1,e%pnode
      
         elrhs(1,inode) = elrhs(1,inode) + (gprep(1)*e%cartd(1,inode) + gprep(3)*e%cartd(2,inode))*aux1
         elrhs(2,inode) = elrhs(2,inode) + (gprep(3)*e%cartd(1,inode) + gprep(2)*e%cartd(2,inode))*aux1           
            
      end do
   end subroutine   

   subroutine supf_elmrhu_oss3d(e,auxtens,tisig,dvolu,gprep,beta,elrhs)
      !This subroutine computes the rhs in nsf1_elmope, additional terms involving the pressure and stress
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,tisig,gprep(auxtens+e%ndime+1),beta
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      real(rp)    :: aux1
      integer(ip) :: inode
 
      aux1=-tisig*dvolu
      
      do inode=1,e%pnode
      
         elrhs(1,inode) = elrhs(1,inode) &
            + (gprep(1)*e%cartd(1,inode) + gprep(6)*e%cartd(2,inode) + gprep(5)*e%cartd(3,inode))*aux1
         elrhs(2,inode) = elrhs(2,inode) &
            + (gprep(6)*e%cartd(1,inode) + gprep(2)*e%cartd(2,inode) + gprep(4)*e%cartd(3,inode))*aux1 
         elrhs(3,inode) = elrhs(3,inode) &
            + (gprep(5)*e%cartd(1,inode) + gprep(4)*e%cartd(2,inode) + gprep(3)*e%cartd(3,inode))*aux1            
            
      end do
   end subroutine    
   
   
   subroutine supf_stressdivergence(e,ntens,elvar,divvar)   
   implicit none
   
   class(FiniteElement) :: e
   integer(ip)             :: ntens
   real(rp), intent(in)    :: elvar(ntens,e%mnode)
   real(rp), intent(out)   :: divvar(e%ndime)
   integer(ip)             :: idime,inode
   
   
   divvar = 0.0_rp
   
   if(e%ndime==2)then
      do inode=1,e%pnode
       divvar(1) =  divvar(1) + e%cartd(1,inode)*elvar(1,e%mnode) + e%cartd(2,inode)*elvar(3,e%mnode)
       divvar(2) =  divvar(2) + e%cartd(1,inode)*elvar(3,e%mnode) + e%cartd(2,inode)*elvar(2,e%mnode)
      end do
   elseif(e%ndime==3)then
      do inode=1,e%pnode
       divvar(1) =  divvar(1) + e%cartd(1,inode)*elvar(1,e%mnode) + e%cartd(2,inode)*elvar(6,e%mnode) &
            + e%cartd(3,inode)*elvar(5,e%mnode)
       divvar(2) =  divvar(2) + e%cartd(1,inode)*elvar(6,e%mnode) + e%cartd(2,inode)*elvar(2,e%mnode) &
            + e%cartd(3,inode)*elvar(4,e%mnode)
       divvar(3) =  divvar(3) + e%cartd(1,inode)*elvar(5,e%mnode) + e%cartd(2,inode)*elvar(4,e%mnode) &
            + e%cartd(3,inode)*elvar(3,e%mnode)            
      end do   
   end if
   
   end subroutine
   
   
   subroutine supf_stressdivergence_boundary(e,ntens,elvar,divvar)   
   implicit none
   
   class(FiniteElement) :: e
   integer(ip)             :: ntens
   real(rp), intent(in)    :: elvar(ntens,e%pnodb)
   real(rp), intent(out)   :: divvar(e%ndime)
   integer(ip)             :: idime,inode
   
   
   divvar = 0.0_rp
   
   if(e%ndime==2)then
      do inode=1,e%pnodb
       divvar(1) =  divvar(1) + e%cartb(1,inode)*elvar(1,e%mnode) + e%cartb(2,inode)*elvar(3,e%mnode)
       divvar(2) =  divvar(2) + e%cartb(1,inode)*elvar(3,e%mnode) + e%cartb(2,inode)*elvar(2,e%mnode)
      end do
   elseif(e%ndime==3)then
      do inode=1,e%pnodb
       divvar(1) =  divvar(1) + e%cartb(1,inode)*elvar(1,e%mnode) + e%cartb(2,inode)*elvar(6,e%mnode) &
            + e%cartb(3,inode)*elvar(5,e%mnode)
       divvar(2) =  divvar(2) + e%cartb(1,inode)*elvar(6,e%mnode) + e%cartb(2,inode)*elvar(2,e%mnode) &
            + e%cartb(3,inode)*elvar(4,e%mnode)
       divvar(3) =  divvar(3) + e%cartb(1,inode)*elvar(5,e%mnode) + e%cartb(2,inode)*elvar(4,e%mnode) &
            + e%cartb(3,inode)*elvar(3,e%mnode)            
      end do   
   end if
   
   end subroutine   
   
   
   !---------------------------------------------------------------------
   !Second Step "First Approach" Split Momentum equation
   !---------------------------------------------------------------------     
   
   !Galerkin terms
   subroutine supf_elmbstGal(e,auxconv,auxtrac,acvis,auxVE,auxG,auxPTT,dtinv,dvolu,auxtens,gpadv,grvel,gpsig,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmst matrix
    !    1/2mu(S,T)+tau1*(div(T),div(S))-tau3*(1/2mu*(T),1/2mu*(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxconv,auxtrac
      real(rp),    intent(in)    :: dvolu,acvis,auxG,auxVE,dtinv
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime),gpsig(auxtens)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1
      real(rp)                   :: Res11s11,Res11s12,Res11s22,Res12s11,Res12s12,Res12S22, &
                                    Res22s11,Res22s12,Res22S22                                    

 
      do jnode=1,e%pnode      
      
         Res11s11 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac*2.0_rp*grvel(1,1)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(1))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res11s12 = (auxVE*auxtrac)*(-2.0_rp)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(3))*e%shape(jnode,e%igaus)
         Res11s22 = 0.0_rp
         
         Res12s11 = -(auxVE*auxtrac)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(3)*e%shape(jnode,e%igaus)
         Res12s12 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               + auxVE*(auxG/acvis)*(gpsig(1) + gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res12s22 = -(auxVE*auxtrac)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(3)*e%shape(jnode,e%igaus)          
         
         Res22s11 = 0.0_rp
         Res22s12 = (auxVE*auxtrac)*(-2.0_rp)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(3))*e%shape(jnode,e%igaus)
         Res22s22 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac*2.0_rp*grvel(2,2)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)      
      
         do inode=1,e%pnode
         
            aux1=e%shape(inode,e%igaus)*dvolu
      
            elmat(1,inode,1,jnode) = aux1*Res11s11 + elmat(1,inode,1,jnode)
            elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
            elmat(1,inode,3,jnode) = aux1*Res11s12 + elmat(1,inode,3,jnode)
            
            elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
            elmat(2,inode,2,jnode) = aux1*Res22s22 + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = aux1*Res22s12 + elmat(2,inode,3,jnode)
            
            elmat(3,inode,1,jnode) = aux1*(2.0_rp*Res12s11) + elmat(3,inode,1,jnode)
            elmat(3,inode,2,jnode) = aux1*(2.0_rp*Res12s22) + elmat(3,inode,2,jnode)                                         
            elmat(3,inode,3,jnode) = aux1*(2.0_rp*Res12s12) + elmat(3,inode,3,jnode)                                  

         end do
      end do
   end subroutine
   
   subroutine supf_elmbstGal3d(e,auxconv,auxtrac1,acvis,auxVE,auxG,auxPTT,dtinv,dvolu,auxtens,gpadv,grvel,gpsig,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmst matrix
    !    1/2mu(S,T)+tau1*(div(T),div(S))-tau3*(1/2mu*(T),1/2mu*(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxconv,auxtrac1
      real(rp),    intent(in)    :: dvolu,acvis,auxG,auxVE,dtinv
      real(rp),    intent(in)    :: gpadv(e%pnode),grvel(e%ndime,e%ndime),gpsig(auxtens)
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1
      real(rp)                   :: Res11s11,Res11s12,Res11s13,Res11s22,Res11s23,Res11s33, &
                                    Res22s11,Res22s12,Res22s13,Res22s22,Res22s23,Res22s33, &
                                    Res33s11,Res33s12,Res33s13,Res33s22,Res33s23,Res33s33, &
                                    Res12s11,Res12s12,Res12s13,Res12s22,Res12s23,Res12s33, &
                                    Res13s11,Res13s12,Res13s13,Res13s22,Res13s23,Res13s33, &
                                    Res23s11,Res23s12,Res23s13,Res23s22,Res23s23,Res23s33
 
      do jnode=1,e%pnode      
      
         Res11s11 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(1,1)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(1))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res11s12 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(6))*e%shape(jnode,e%igaus)
         Res11s13 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(1,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(5))*e%shape(jnode,e%igaus)
         Res11s22 = 0.0_rp
         Res11s23 = 0.0_rp
         Res11s33 = 0.0_rp
         
         Res22s11 = 0.0_rp
         Res22s12 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(6))*e%shape(jnode,e%igaus)
         Res22s13 = 0.0_rp         
         Res22s22 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(2,2)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)   
         Res22s23 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(2,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(4))*e%shape(jnode,e%igaus)
         Res22s33 = 0.0_rp
         
         Res33s11 = 0.0_rp
         Res33s12 = 0.0_rp
         Res33s13 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(3,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(5))*e%shape(jnode,e%igaus)
         Res33s22 = 0.0_rp
         Res33s23 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(3,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(4))*e%shape(jnode,e%igaus)
         Res33s33 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(3,3)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(3))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
            
         Res12s11 = -(auxVE*auxtrac1)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)
         Res12s12 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               - auxVE*(auxtrac1*(grvel(2,2) + grvel(1,1))*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(gpsig(1) + gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res12s13 = -(auxVE*auxtrac1)*grvel(2,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)
         Res12s22 = -(auxVE*auxtrac1)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)          
         Res12s23 = -(auxVE*auxtrac1)*grvel(1,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)
         Res12s33 = 0.0_rp
         
         Res13s11 = -(auxVE*auxtrac1)*grvel(3,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)
         Res13s12 = -(auxVE*auxtrac1)*grvel(3,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)
         Res13s13 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               - auxVE*(auxtrac1*(grvel(3,3) + grvel(1,1))*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(gpsig(1) + gpsig(3))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res13s22 = 0.0_rp        
         Res13s23 = -(auxVE*auxtrac1)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)
         Res13s33 = -(auxVE*auxtrac1)*grvel(1,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)         

         Res23s11 = 0.0_rp
         Res23s12 = -(auxVE*auxtrac1)*grvel(3,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)
         Res23s13 = -(auxVE*auxtrac1)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)
         Res23s22 = -(auxVE*auxtrac1)*grvel(3,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)         
         Res23s23 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               - auxVE*(auxtrac1*(grvel(3,3) + grvel(2,2))*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(gpsig(2) + gpsig(3))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)       
         Res23s33 = -(auxVE*auxtrac1)*grvel(2,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)         
         
      
         do inode=1,e%pnode
         
            aux1=e%shape(inode,e%igaus)*dvolu
      
            elmat(1,inode,1,jnode) = aux1*Res11s11 + elmat(1,inode,1,jnode)
            elmat(1,inode,2,jnode) = aux1*Res11s22 + elmat(1,inode,2,jnode)
            elmat(1,inode,3,jnode) = aux1*Res11s33 + elmat(1,inode,3,jnode)
            elmat(1,inode,4,jnode) = aux1*Res11s23 + elmat(1,inode,4,jnode)
            elmat(1,inode,5,jnode) = aux1*Res11s13 + elmat(1,inode,5,jnode)            
            elmat(1,inode,6,jnode) = aux1*Res11s12 + elmat(1,inode,6,jnode)
            
            elmat(2,inode,1,jnode) = aux1*Res22s11 + elmat(2,inode,1,jnode)
            elmat(2,inode,2,jnode) = aux1*Res22s22 + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = aux1*Res22s33 + elmat(2,inode,3,jnode)
            elmat(2,inode,4,jnode) = aux1*Res22s23 + elmat(2,inode,4,jnode)
            elmat(2,inode,5,jnode) = aux1*Res22s13 + elmat(2,inode,5,jnode)            
            elmat(2,inode,6,jnode) = aux1*Res22s12 + elmat(2,inode,6,jnode)
            
            elmat(3,inode,1,jnode) = aux1*Res33s11 + elmat(3,inode,1,jnode)
            elmat(3,inode,2,jnode) = aux1*Res33s22 + elmat(3,inode,2,jnode)
            elmat(3,inode,3,jnode) = aux1*Res33s33 + elmat(3,inode,3,jnode)
            elmat(3,inode,4,jnode) = aux1*Res33s23 + elmat(3,inode,4,jnode)
            elmat(3,inode,5,jnode) = aux1*Res33s13 + elmat(3,inode,5,jnode)            
            elmat(3,inode,6,jnode) = aux1*Res33s12 + elmat(3,inode,6,jnode)
            
            elmat(4,inode,1,jnode) = aux1*(2.0_rp*Res23s11) + elmat(4,inode,1,jnode)
            elmat(4,inode,2,jnode) = aux1*(2.0_rp*Res23s22) + elmat(4,inode,2,jnode)                                         
            elmat(4,inode,3,jnode) = aux1*(2.0_rp*Res23s33) + elmat(4,inode,3,jnode)               
            elmat(4,inode,4,jnode) = aux1*(2.0_rp*Res23s23) + elmat(4,inode,4,jnode)
            elmat(4,inode,5,jnode) = aux1*(2.0_rp*Res23s13) + elmat(4,inode,5,jnode)                                         
            elmat(4,inode,6,jnode) = aux1*(2.0_rp*Res23s12) + elmat(4,inode,6,jnode)  

            elmat(5,inode,1,jnode) = aux1*(2.0_rp*Res13s11) + elmat(5,inode,1,jnode)
            elmat(5,inode,2,jnode) = aux1*(2.0_rp*Res13s22) + elmat(5,inode,2,jnode)                                         
            elmat(5,inode,3,jnode) = aux1*(2.0_rp*Res13s33) + elmat(5,inode,3,jnode)               
            elmat(5,inode,4,jnode) = aux1*(2.0_rp*Res13s23) + elmat(5,inode,4,jnode)
            elmat(5,inode,5,jnode) = aux1*(2.0_rp*Res13s13) + elmat(5,inode,5,jnode)                                         
            elmat(5,inode,6,jnode) = aux1*(2.0_rp*Res13s12) + elmat(5,inode,6,jnode)              

            elmat(6,inode,1,jnode) = aux1*(2.0_rp*Res12s11) + elmat(6,inode,1,jnode)
            elmat(6,inode,2,jnode) = aux1*(2.0_rp*Res12s22) + elmat(6,inode,2,jnode)                                         
            elmat(6,inode,3,jnode) = aux1*(2.0_rp*Res12s33) + elmat(6,inode,3,jnode)               
            elmat(6,inode,4,jnode) = aux1*(2.0_rp*Res12s23) + elmat(6,inode,4,jnode)
            elmat(6,inode,5,jnode) = aux1*(2.0_rp*Res12s13) + elmat(6,inode,5,jnode)                                         
            elmat(6,inode,6,jnode) = aux1*(2.0_rp*Res12s12) + elmat(6,inode,6,jnode)              
            
         end do
      end do    

   end subroutine supf_elmbstGal3d 
   
   subroutine supf_elmbst_splitoss(e,dvolu,timom,auxtens,beta,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmst matrix
    !    1/2mu(S,T)+tau1*(div(T),div(S))-tau3*(1/2mu*(T),1/2mu*(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,timom,beta
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1    
      
      aux1=(1.0_rp-beta)*timom*dvolu 
      do jnode=1,e%pnode
         do inode=1,e%pnode
         
            elmat(1,inode,1,jnode) = aux1*e%cartd(1,jnode)*e%cartd(1,inode) + elmat(1,inode,1,jnode)
            elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
            elmat(1,inode,3,jnode) = aux1*e%cartd(2,jnode)*e%cartd(1,inode) + elmat(1,inode,3,jnode)
            
            elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
            elmat(2,inode,2,jnode) = aux1*e%cartd(2,jnode)*e%cartd(2,inode) + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = aux1*e%cartd(1,jnode)*e%cartd(2,inode) + elmat(2,inode,3,jnode)
            
            elmat(3,inode,1,jnode) = aux1*e%cartd(1,jnode)*e%cartd(2,inode) + elmat(3,inode,1,jnode)
            elmat(3,inode,2,jnode) = aux1*e%cartd(2,jnode)*e%cartd(1,inode) + elmat(3,inode,2,jnode)                                         
            elmat(3,inode,3,jnode) = aux1*(e%cartd(1,jnode)*e%cartd(1,inode)+e%cartd(2,jnode)*e%cartd(2,inode)) &
                  + elmat(3,inode,3,jnode)                                  

         end do
      end do   

   end subroutine supf_elmbst_splitoss
   
   subroutine supf_elmbst_splitoss3d(e,dvolu,timom,auxtens,beta,elmat)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Galerkin component of elmst matrix
    !    1/2mu(S,T)+tau1*(div(T),div(S))-tau3*(1/2mu*(T),1/2mu*(S))
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,timom,beta
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode
      real(rp)                   :: aux1   
 
      do jnode=1,e%pnode      
         do inode=1,e%pnode
         
            aux1=(1.0_rp-beta)*timom*dvolu
      
            elmat(1,inode,1,jnode) = aux1*e%cartd(1,jnode)*e%cartd(1,inode) + elmat(1,inode,1,jnode)
            elmat(1,inode,2,jnode) = 0.0_rp + elmat(1,inode,2,jnode)
            elmat(1,inode,3,jnode) = 0.0_rp + elmat(1,inode,3,jnode)
            elmat(1,inode,4,jnode) = 0.0_rp + elmat(1,inode,4,jnode)
            elmat(1,inode,5,jnode) = aux1*e%cartd(3,jnode)*e%cartd(1,inode) + elmat(1,inode,5,jnode)            
            elmat(1,inode,6,jnode) = aux1*e%cartd(2,jnode)*e%cartd(1,inode) + elmat(1,inode,6,jnode)
            
            elmat(2,inode,1,jnode) = 0.0_rp + elmat(2,inode,1,jnode)
            elmat(2,inode,2,jnode) = aux1*e%cartd(2,jnode)*e%cartd(2,inode) + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = 0.0_rp + elmat(2,inode,3,jnode)
            elmat(2,inode,4,jnode) = aux1*e%cartd(3,jnode)*e%cartd(2,inode) + elmat(2,inode,4,jnode)
            elmat(2,inode,5,jnode) = 0.0_rp + elmat(2,inode,5,jnode)            
            elmat(2,inode,6,jnode) = aux1*e%cartd(1,jnode)*e%cartd(2,inode) + elmat(2,inode,6,jnode)
            
            elmat(3,inode,1,jnode) = 0.0_rp + elmat(3,inode,1,jnode)
            elmat(3,inode,2,jnode) = 0.0_rp + elmat(3,inode,2,jnode)
            elmat(3,inode,3,jnode) = aux1*e%cartd(3,jnode)*e%cartd(3,inode) + elmat(3,inode,3,jnode)
            elmat(3,inode,4,jnode) = aux1*e%cartd(2,jnode)*e%cartd(3,inode) + elmat(3,inode,4,jnode)
            elmat(3,inode,5,jnode) = aux1*e%cartd(1,jnode)*e%cartd(3,inode) + elmat(3,inode,5,jnode)            
            elmat(3,inode,6,jnode) = 0.0_rp + elmat(3,inode,6,jnode)
            
            elmat(4,inode,1,jnode) = 0.0_rp + elmat(4,inode,1,jnode)
            elmat(4,inode,2,jnode) = aux1*e%cartd(2,jnode)*e%cartd(3,inode) + elmat(4,inode,2,jnode)                                         
            elmat(4,inode,3,jnode) = aux1*e%cartd(3,jnode)*e%cartd(2,inode) + elmat(4,inode,3,jnode)               
            elmat(4,inode,4,jnode) = aux1*(e%cartd(3,jnode)*e%cartd(3,inode) + e%cartd(2,jnode)*e%cartd(2,inode)) &
                  + elmat(4,inode,4,jnode)
            elmat(4,inode,5,jnode) = aux1*e%cartd(1,jnode)*e%cartd(2,inode) + elmat(4,inode,5,jnode)                                         
            elmat(4,inode,6,jnode) = aux1*e%cartd(1,jnode)*e%cartd(3,inode) + elmat(4,inode,6,jnode)  

            elmat(5,inode,1,jnode) = aux1*e%cartd(1,jnode)*e%cartd(3,inode) + elmat(5,inode,1,jnode)
            elmat(5,inode,2,jnode) = 0.0_rp + elmat(5,inode,2,jnode)                                         
            elmat(5,inode,3,jnode) = aux1*e%cartd(3,jnode)*e%cartd(1,inode) + elmat(5,inode,3,jnode)               
            elmat(5,inode,4,jnode) = aux1*e%cartd(2,jnode)*e%cartd(1,inode) + elmat(5,inode,4,jnode)
            elmat(5,inode,5,jnode) = aux1*(e%cartd(1,jnode)*e%cartd(1,inode) + e%cartd(3,jnode)*e%cartd(3,inode)) &
                  + elmat(5,inode,5,jnode)                                         
            elmat(5,inode,6,jnode) = aux1*e%cartd(2,jnode)*e%cartd(3,inode) + elmat(5,inode,6,jnode)              

            elmat(6,inode,1,jnode) = aux1*e%cartd(1,jnode)*e%cartd(2,inode) + elmat(6,inode,1,jnode)
            elmat(6,inode,2,jnode) = aux1*e%cartd(2,jnode)*e%cartd(1,inode) + elmat(6,inode,2,jnode)                                         
            elmat(6,inode,3,jnode) = 0.0_rp + elmat(6,inode,3,jnode)               
            elmat(6,inode,4,jnode) = aux1*e%cartd(3,jnode)*e%cartd(1,inode) + elmat(6,inode,4,jnode)
            elmat(6,inode,5,jnode) = aux1*e%cartd(3,jnode)*e%cartd(2,inode) + elmat(6,inode,5,jnode)                                         
            elmat(6,inode,6,jnode) = aux1*(e%cartd(1,jnode)*e%cartd(1,inode) + e%cartd(2,jnode)*e%cartd(2,inode)) &
                  + elmat(6,inode,6,jnode)              
            
         end do
      end do   

   end subroutine supf_elmbst_splitoss3d    
   
   subroutine supf_elmbstEst(e,auxGie,auxconv,auxtrac1,auxtrac2,auxVE,auxG,auxPTT,acvis,tisig,dtinv,gpadv,gpsig,dvolu,auxtens,grvel,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in constitutive equation
    !    (-(1/2eta)*tau + auxVE*a*grad(tau),R_constituitva)
    !   +(auxVE*(2*grad(a)*tau -(auxG/eta)*sigma*tau),R_constituitva)
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxconv,auxtrac1,auxtrac2,auxGie
      real(rp),    intent(in)    :: acvis,auxVE,auxG,dvolu,gpsig(auxtens),gpadv(e%pnode),dtinv
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),tisig     
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode,auxstab,aux1,auxstab2,aux4
      real(rp)                   :: aux2,aux3
      real(rp)                   :: Res11s11,Res11s12,Res11s22,Res12s11,Res12s12,Res12S22, &
                                    Res22s11,Res22s12,Res22S22                                     
      
      do jnode=1,e%pnode  
         
         Res11s11 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(1,1)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(1))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res11s12 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(3))*e%shape(jnode,e%igaus)
         Res11s22 = 0.0_rp
         
         Res12s11 = -(auxVE*auxtrac1)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(3)*e%shape(jnode,e%igaus)
         Res12s12 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               + auxVE*(auxG/acvis)*(gpsig(1) + gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res12s22 = -(auxVE*auxtrac1)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(3)*e%shape(jnode,e%igaus)          
         
         Res22s11 = 0.0_rp
         Res22s12 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(3))*e%shape(jnode,e%igaus)
         Res22s22 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(2,2)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)      
             
         do inode=1,e%pnode
            
            aux2= (e%shape(inode,e%igaus)*(-tisig/(2.0_rp*acvis)) + (tisig*auxVE)*gpadv(inode))*dvolu          

            aux3= ((tisig*auxVE)*e%shape(inode,e%igaus))*dvolu         
               
            elmat(1,inode,1,jnode) = aux2*Res11s11 &
                  + aux3*(auxtrac2*2.0_rp*grvel(1,1))*Res11s11 - aux3*((auxGie*auxG/acvis)*gpsig(1))*Res11s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s11 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res12s11 &                
                  + elmat(1,inode,1,jnode)
                  
            elmat(1,inode,2,jnode) = 0.0_rp &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s22 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res12s22 &                      
                  + elmat(1,inode,2,jnode)
                  
            elmat(1,inode,3,jnode) = aux2*Res11s12 &
                  + aux3*(auxtrac2*2.0_rp*grvel(1,1))*Res11s12 - aux3*((auxGie*auxG/acvis)*gpsig(1))*Res11s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s12 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res12s12 &                     
                  + elmat(1,inode,3,jnode)
            
            elmat(2,inode,1,jnode) = 0.0_rp &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s11 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res12s11 &             
                  + elmat(2,inode,1,jnode)
            
            elmat(2,inode,2,jnode) = aux2*Res22s22 &
                  + aux3*(auxtrac2*2.0_rp*grvel(2,2))*Res22s22 - aux3*((auxGie*auxG/acvis)*gpsig(2))*Res22s22 &              
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s22 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res12s22 &                     
                  + elmat(2,inode,2,jnode)
                  
            elmat(2,inode,3,jnode) = aux2*Res22s12 &
                  + aux3*(auxtrac2*2.0_rp*grvel(2,2))*Res22s12 - aux3*((auxGie*auxG/acvis)*gpsig(2))*Res22s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s12 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res12s12 &                      
                  + elmat(2,inode,3,jnode)            
                 
            elmat(3,inode,1,jnode) = aux2*(2.0_rp*Res12s11) & 
                  + aux3*(auxtrac2*2.0_rp*grvel(2,1))*Res11s11 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res11s11 &              
                  - aux3*((auxGie*auxG/acvis)*(gpsig(1)+gpsig(2)))*(Res12s11) &                            
                  + elmat(3,inode,1,jnode)
                  
                  
            elmat(3,inode,2,jnode) = aux2*(2.0_rp*Res12s22) &
                  - aux3*((auxGie*auxG/acvis)*(gpsig(1)+gpsig(2)))*(Res12s22) &                
                  + aux3*(auxtrac2*2.0_rp*grvel(1,2))*Res22s22 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res22s22 &                  
                  + elmat(3,inode,2,jnode)           
            
                  
            elmat(3,inode,3,jnode) = aux2*(2.0_rp*Res12s12) &
                  + aux3*(auxtrac2*2.0_rp*grvel(2,1))*Res11s12 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res11s12 &               
                  + aux3*(auxtrac2*2.0_rp*grvel(1,2))*Res22s12 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res22s12 &              
                  - aux3*((auxGie*auxG/acvis)*(gpsig(1)+gpsig(2)))*(Res12s12) &                         
                  + elmat(3,inode,3,jnode)                    

         end do
      end do   

   end subroutine supf_elmbstEst
   
   subroutine supf_elmbstEst3d(e,auxconv,auxtrac1,auxtrac2,auxVE,auxG,auxPTT,acvis,tisig,dtinv,gpadv,gpsig,dvolu,auxtens,grvel,elmat) 
    !-----------------------------------------------------------------------
    !
    ! This routine computes the first lhs term for ASGS in constitutive equation
    !    (-(1/2eta)*tau + auxVE*a*grad(tau),R_constituitva)
    !   +(auxVE*(2*grad(a)*tau -(auxG/eta)*sigma*tau),R_constituitva)
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxconv,auxtrac1,auxtrac2
      real(rp),    intent(in)    :: acvis,auxVE,auxG,dvolu,gpsig(auxtens),gpadv(e%pnode),dtinv
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),tisig     
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,jnode,aux1,aux4
      real(rp)                   :: aux2,aux3
      real(rp)                   :: Res11s11,Res11s12,Res11s13,Res11s22,Res11s23,Res11S33, &
                                    Res22s11,Res22s12,Res22s13,Res22s22,Res22s23,Res22S33, & 
                                    Res33s11,Res33s12,Res33s13,Res33s22,Res33s23,Res33S33, & 
                                    Res12s11,Res12s12,Res12s13,Res12s22,Res12s23,Res12S33, & 
                                    Res13s11,Res13s12,Res13s13,Res13s22,Res13s23,Res13S33, & 
                                    Res23s11,Res23s12,Res23s13,Res23s22,Res23s23,Res23S33 
      
      do jnode=1,e%pnode  
         
         Res11s11 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(1,1)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(1))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res11s12 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(6))*e%shape(jnode,e%igaus)
         Res11s13 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(1,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(5))*e%shape(jnode,e%igaus)
         Res11s22 = 0.0_rp
         Res11s23 = 0.0_rp
         Res11s33 = 0.0_rp
         
         Res22s11 = 0.0_rp
         Res22s12 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(6))*e%shape(jnode,e%igaus)
         Res22s13 = 0.0_rp         
         Res22s22 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(2,2)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)   
         Res22s23 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(2,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(4))*e%shape(jnode,e%igaus)
         Res22s33 = 0.0_rp
         
         Res33s11 = 0.0_rp
         Res33s12 = 0.0_rp
         Res33s13 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(3,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(5))*e%shape(jnode,e%igaus)
         Res33s22 = 0.0_rp
         Res33s23 = (auxVE*auxtrac1)*(-2.0_rp)*grvel(3,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*(2.0_rp*gpsig(4))*e%shape(jnode,e%igaus)
         Res33s33 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode) - auxtrac1*2.0_rp*grvel(3,3)*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(2.0_rp*gpsig(3))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
            
         Res12s11 = -(auxVE*auxtrac1)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)
         Res12s12 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               - auxVE*(auxtrac1*(grvel(2,2) + grvel(1,1))*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(gpsig(1) + gpsig(2))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res12s13 = -(auxVE*auxtrac1)*grvel(2,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)
         Res12s22 = -(auxVE*auxtrac1)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)          
         Res12s23 = -(auxVE*auxtrac1)*grvel(1,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)
         Res12s33 = 0.0_rp
         
         Res13s11 = -(auxVE*auxtrac1)*grvel(3,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)
         Res13s12 = -(auxVE*auxtrac1)*grvel(3,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)
         Res13s13 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               - auxVE*(auxtrac1*(grvel(3,3) + grvel(1,1))*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(gpsig(1) + gpsig(3))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)
         Res13s22 = 0.0_rp        
         Res13s23 = -(auxVE*auxtrac1)*grvel(1,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)
         Res13s33 = -(auxVE*auxtrac1)*grvel(1,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)         

         Res23s11 = 0.0_rp
         Res23s12 = -(auxVE*auxtrac1)*grvel(3,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(5)*e%shape(jnode,e%igaus)
         Res23s13 = -(auxVE*auxtrac1)*grvel(2,1)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(6)*e%shape(jnode,e%igaus)
         Res23s22 = -(auxVE*auxtrac1)*grvel(3,2)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)         
         Res23s23 = e%shape(jnode,e%igaus)/(2.0_rp*acvis) + auxVE*(auxconv*gpadv(jnode)) &
               - auxVE*(auxtrac1*(grvel(3,3) + grvel(2,2))*e%shape(jnode,e%igaus)) &
               + auxVE*(auxG/acvis)*(gpsig(2) + gpsig(3))*e%shape(jnode,e%igaus) + auxVE*dtinv*e%shape(jnode,e%igaus)       
         Res23s33 = -(auxVE*auxtrac1)*grvel(2,3)*e%shape(jnode,e%igaus) + auxVE*(auxG/acvis)*gpsig(4)*e%shape(jnode,e%igaus)         
         
         do inode=1,e%pnode
         
            aux2= (e%shape(inode,e%igaus)*(-tisig/(2.0_rp*acvis)) + (tisig*auxVE)*gpadv(inode))*dvolu          

            aux3= ((tisig*auxVE)*e%shape(inode,e%igaus))*dvolu         
               
            elmat(1,inode,1,jnode) = aux2*Res11s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11s11 - aux3*((auxG/acvis)*gpsig(1))*Res11s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s11 - aux3*((auxG/acvis)*gpsig(6))*Res12s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13s11 - aux3*((auxG/acvis)*gpsig(5))*Res13s11 &
                  + elmat(1,inode,1,jnode)               
            elmat(1,inode,2,jnode) = aux2*Res11s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11s22 - aux3*((auxG/acvis)*gpsig(1))*Res11s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s22 - aux3*((auxG/acvis)*gpsig(6))*Res12s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13s22 - aux3*((auxG/acvis)*gpsig(5))*Res13s22 &                    
                  + elmat(1,inode,2,jnode)               
            elmat(1,inode,3,jnode) = aux2*Res11s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11s33 - aux3*((auxG/acvis)*gpsig(1))*Res11s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s33 - aux3*((auxG/acvis)*gpsig(6))*Res12s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13s33 - aux3*((auxG/acvis)*gpsig(5))*Res13s33 &                    
                  + elmat(1,inode,3,jnode)
            elmat(1,inode,4,jnode) = aux2*Res11s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11s23 - aux3*((auxG/acvis)*gpsig(1))*Res11s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s23 - aux3*((auxG/acvis)*gpsig(6))*Res12s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13s23 - aux3*((auxG/acvis)*gpsig(5))*Res13s23 &                    
                  + elmat(1,inode,4,jnode)
            elmat(1,inode,5,jnode) = aux2*Res11s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11s13 - aux3*((auxG/acvis)*gpsig(1))*Res11s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s13 - aux3*((auxG/acvis)*gpsig(6))*Res12s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13s13 - aux3*((auxG/acvis)*gpsig(5))*Res13s13 &                    
                  + elmat(1,inode,5,jnode)               
            elmat(1,inode,6,jnode) = aux2*Res11s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11s12 - aux3*((auxG/acvis)*gpsig(1))*Res11s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12s12 - aux3*((auxG/acvis)*gpsig(6))*Res12s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13s12 - aux3*((auxG/acvis)*gpsig(5))*Res13s12 &                    
                  + elmat(1,inode,6,jnode)
            
            elmat(2,inode,1,jnode) = aux2*Res22s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22s11 - aux3*((auxG/acvis)*gpsig(2))*Res22s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s11 - aux3*((auxG/acvis)*gpsig(6))*Res12s11 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23s11 - aux3*((auxG/acvis)*gpsig(4))*Res23s11 &                
                  + elmat(2,inode,1,jnode)         
            elmat(2,inode,2,jnode) = aux2*Res22s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22s22 - aux3*((auxG/acvis)*gpsig(2))*Res22s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s22 - aux3*((auxG/acvis)*gpsig(6))*Res12s22 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23s22 - aux3*((auxG/acvis)*gpsig(4))*Res23s22 &                    
                  + elmat(2,inode,2,jnode)
            elmat(2,inode,3,jnode) = aux2*Res22s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22s33 - aux3*((auxG/acvis)*gpsig(2))*Res22s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s33 - aux3*((auxG/acvis)*gpsig(6))*Res12s33 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23s33 - aux3*((auxG/acvis)*gpsig(4))*Res23s33 &                    
                  + elmat(2,inode,3,jnode)               
            elmat(2,inode,4,jnode) = aux2*Res22s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22s23 - aux3*((auxG/acvis)*gpsig(2))*Res22s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s23 - aux3*((auxG/acvis)*gpsig(6))*Res12s23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23s23 - aux3*((auxG/acvis)*gpsig(4))*Res23s23 &                    
                  + elmat(2,inode,4,jnode)
            elmat(2,inode,5,jnode) = aux2*Res22s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22s13 - aux3*((auxG/acvis)*gpsig(2))*Res22s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s13 - aux3*((auxG/acvis)*gpsig(6))*Res12s13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23s13 - aux3*((auxG/acvis)*gpsig(4))*Res23s13 &                    
                  + elmat(2,inode,5,jnode)               
            elmat(2,inode,6,jnode) = aux2*Res22s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22s12 - aux3*((auxG/acvis)*gpsig(2))*Res22s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12s12 - aux3*((auxG/acvis)*gpsig(6))*Res12s12 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23s12 - aux3*((auxG/acvis)*gpsig(4))*Res23s12 &                    
                  + elmat(2,inode,6,jnode)
                  
            elmat(3,inode,1,jnode) = aux2*Res33s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33s11 - aux3*((auxG/acvis)*gpsig(3))*Res33s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13s11 - aux3*((auxG/acvis)*gpsig(5))*Res13s11 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23s11 - aux3*((auxG/acvis)*gpsig(4))*Res23s11 &                
                  + elmat(3,inode,1,jnode)         
            elmat(3,inode,2,jnode) = aux2*Res33s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33s22 - aux3*((auxG/acvis)*gpsig(3))*Res33s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13s22 - aux3*((auxG/acvis)*gpsig(5))*Res13s22 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23s22 - aux3*((auxG/acvis)*gpsig(4))*Res23s22 &                
                  + elmat(3,inode,2,jnode)
            elmat(3,inode,3,jnode) = aux2*Res33s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33s33 - aux3*((auxG/acvis)*gpsig(3))*Res33s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13s33 - aux3*((auxG/acvis)*gpsig(5))*Res13s33 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23s33 - aux3*((auxG/acvis)*gpsig(4))*Res23s33 &                
                  + elmat(3,inode,3,jnode)                
            elmat(3,inode,4,jnode) = aux2*Res33s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33s23 - aux3*((auxG/acvis)*gpsig(3))*Res33s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13s23 - aux3*((auxG/acvis)*gpsig(5))*Res13s23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23s23 - aux3*((auxG/acvis)*gpsig(4))*Res23s23 &                
                  + elmat(3,inode,4,jnode) 
            elmat(3,inode,5,jnode) = aux2*Res33s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33s13 - aux3*((auxG/acvis)*gpsig(3))*Res33s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13s13 - aux3*((auxG/acvis)*gpsig(5))*Res13s13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23s13 - aux3*((auxG/acvis)*gpsig(4))*Res23s13 &                
                  + elmat(3,inode,5,jnode) 
            elmat(3,inode,6,jnode) = aux2*Res33s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33s12 - aux3*((auxG/acvis)*gpsig(3))*Res33s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13s12 - aux3*((auxG/acvis)*gpsig(5))*Res13s12 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23s12 - aux3*((auxG/acvis)*gpsig(4))*Res23s12 &                
                  + elmat(3,inode,6,jnode)              
                  
            elmat(4,inode,1,jnode) = (2.0_rp*aux2)*Res23s11 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23s11 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12s11 - aux3*((auxG/acvis)*gpsig(5))*Res12s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13s11 - aux3*((auxG/acvis)*gpsig(6))*Res13s11 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22s11 - aux3*((auxG/acvis)*gpsig(4))*Res22s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33s11 - aux3*((auxG/acvis)*gpsig(4))*Res33s11 &              
                  + elmat(4,inode,1,jnode)                
            elmat(4,inode,2,jnode) = (2.0_rp*aux2)*Res23s22 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23s22 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12s22 - aux3*((auxG/acvis)*gpsig(5))*Res12s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13s22 - aux3*((auxG/acvis)*gpsig(6))*Res13s22 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22s22 - aux3*((auxG/acvis)*gpsig(4))*Res22s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33s22 - aux3*((auxG/acvis)*gpsig(4))*Res33s22 &              
                  + elmat(4,inode,2,jnode)
            elmat(4,inode,3,jnode) = (2.0_rp*aux2)*Res23s33 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23s33 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12s33 - aux3*((auxG/acvis)*gpsig(5))*Res12s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13s33 - aux3*((auxG/acvis)*gpsig(6))*Res13s33 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22s33 - aux3*((auxG/acvis)*gpsig(4))*Res22s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33s33 - aux3*((auxG/acvis)*gpsig(4))*Res33s33 &              
                  + elmat(4,inode,3,jnode)               
            elmat(4,inode,4,jnode) = (2.0_rp*aux2)*Res23s23 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23s23 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12s23 - aux3*((auxG/acvis)*gpsig(5))*Res12s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13s23 - aux3*((auxG/acvis)*gpsig(6))*Res13s23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22s23 - aux3*((auxG/acvis)*gpsig(4))*Res22s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33s23 - aux3*((auxG/acvis)*gpsig(4))*Res33s23 &              
                  + elmat(4,inode,4,jnode)               
            elmat(4,inode,5,jnode) = (2.0_rp*aux2)*Res23s13 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23s13 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12s13 - aux3*((auxG/acvis)*gpsig(5))*Res12s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13s13 - aux3*((auxG/acvis)*gpsig(6))*Res13s13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22s13 - aux3*((auxG/acvis)*gpsig(4))*Res22s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33s13 - aux3*((auxG/acvis)*gpsig(4))*Res33s13 &              
                  + elmat(4,inode,5,jnode)
            elmat(4,inode,6,jnode) = (2.0_rp*aux2)*Res23s12 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23s12 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12s12 - aux3*((auxG/acvis)*gpsig(5))*Res12s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13s12 - aux3*((auxG/acvis)*gpsig(6))*Res13s12 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22s12 - aux3*((auxG/acvis)*gpsig(4))*Res22s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33s12 - aux3*((auxG/acvis)*gpsig(4))*Res33s12 &              
                  + elmat(4,inode,6,jnode)               
                  
            elmat(5,inode,1,jnode) = (2.0_rp*aux2)*Res13s11 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13s11 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12s11 - aux3*((auxG/acvis)*gpsig(4))*Res12s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23s11 - aux3*((auxG/acvis)*gpsig(6))*Res23s11 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11s11 - aux3*((auxG/acvis)*gpsig(5))*Res11s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33s11 - aux3*((auxG/acvis)*gpsig(5))*Res33s11 &              
                  + elmat(5,inode,1,jnode)
            elmat(5,inode,2,jnode) = (2.0_rp*aux2)*Res13s22 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13s22 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12s22 - aux3*((auxG/acvis)*gpsig(4))*Res12s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23s22 - aux3*((auxG/acvis)*gpsig(6))*Res23s22 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11s22 - aux3*((auxG/acvis)*gpsig(5))*Res11s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33s22 - aux3*((auxG/acvis)*gpsig(5))*Res33s22 &              
                  + elmat(5,inode,2,jnode)
            elmat(5,inode,3,jnode) = (2.0_rp*aux2)*Res13s33 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13s33 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12s33 - aux3*((auxG/acvis)*gpsig(4))*Res12s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23s33 - aux3*((auxG/acvis)*gpsig(6))*Res23s33 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11s33 - aux3*((auxG/acvis)*gpsig(5))*Res11s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33s33 - aux3*((auxG/acvis)*gpsig(5))*Res33s33 &              
                  + elmat(5,inode,3,jnode)
            elmat(5,inode,4,jnode) = (2.0_rp*aux2)*Res13s23 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13s23 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12s23 - aux3*((auxG/acvis)*gpsig(4))*Res12s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23s23 - aux3*((auxG/acvis)*gpsig(6))*Res23s23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11s23 - aux3*((auxG/acvis)*gpsig(5))*Res11s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33s23 - aux3*((auxG/acvis)*gpsig(5))*Res33s23 &              
                  + elmat(5,inode,4,jnode)
            elmat(5,inode,5,jnode) = (2.0_rp*aux2)*Res13s13 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13s13 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12s13 - aux3*((auxG/acvis)*gpsig(4))*Res12s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23s13 - aux3*((auxG/acvis)*gpsig(6))*Res23s13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11s13 - aux3*((auxG/acvis)*gpsig(5))*Res11s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33s13 - aux3*((auxG/acvis)*gpsig(5))*Res33s13 &              
                  + elmat(5,inode,5,jnode)               
            elmat(5,inode,6,jnode) = (2.0_rp*aux2)*Res13s12 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13s12 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12s12 - aux3*((auxG/acvis)*gpsig(4))*Res12s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23s12 - aux3*((auxG/acvis)*gpsig(6))*Res23s12 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11s12 - aux3*((auxG/acvis)*gpsig(5))*Res11s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33s12 - aux3*((auxG/acvis)*gpsig(5))*Res33s12 &              
                  + elmat(5,inode,6,jnode)
                  
            elmat(6,inode,1,jnode) = (2.0_rp*aux2)*Res12s11 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12s11 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12s11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13s11 - aux3*((auxG/acvis)*gpsig(4))*Res13s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23s11 - aux3*((auxG/acvis)*gpsig(5))*Res23s11 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11s11 - aux3*((auxG/acvis)*gpsig(6))*Res11s11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22s11 - aux3*((auxG/acvis)*gpsig(6))*Res22s11 &              
                  + elmat(6,inode,1,jnode)                
            elmat(6,inode,2,jnode) = (2.0_rp*aux2)*Res12s22 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12s22 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12s22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13s22 - aux3*((auxG/acvis)*gpsig(4))*Res13s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23s22 - aux3*((auxG/acvis)*gpsig(5))*Res23s22 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11s22 - aux3*((auxG/acvis)*gpsig(6))*Res11s22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22s22 - aux3*((auxG/acvis)*gpsig(6))*Res22s22 &              
                  + elmat(6,inode,2,jnode)
            elmat(6,inode,3,jnode) = (2.0_rp*aux2)*Res12s33 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12s33 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12s33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13s33 - aux3*((auxG/acvis)*gpsig(4))*Res13s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23s33 - aux3*((auxG/acvis)*gpsig(5))*Res23s33 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11s33 - aux3*((auxG/acvis)*gpsig(6))*Res11s33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22s33 - aux3*((auxG/acvis)*gpsig(6))*Res22s33 &              
                  + elmat(6,inode,3,jnode)               
            elmat(6,inode,4,jnode) = (2.0_rp*aux2)*Res12s23 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12s23 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12s23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13s23 - aux3*((auxG/acvis)*gpsig(4))*Res13s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23s23 - aux3*((auxG/acvis)*gpsig(5))*Res23s23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11s23 - aux3*((auxG/acvis)*gpsig(6))*Res11s23 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22s23 - aux3*((auxG/acvis)*gpsig(6))*Res22s23 &              
                  + elmat(6,inode,4,jnode)                                  
            elmat(6,inode,5,jnode) = (2.0_rp*aux2)*Res12s13 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12s13 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12s13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13s13 - aux3*((auxG/acvis)*gpsig(4))*Res13s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23s13 - aux3*((auxG/acvis)*gpsig(5))*Res23s13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11s13 - aux3*((auxG/acvis)*gpsig(6))*Res11s13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22s13 - aux3*((auxG/acvis)*gpsig(6))*Res22s13 &              
                  + elmat(6,inode,5,jnode)
            elmat(6,inode,6,jnode) = (2.0_rp*aux2)*Res12s12 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12s12 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12s12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13s12 - aux3*((auxG/acvis)*gpsig(4))*Res13s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23s12 - aux3*((auxG/acvis)*gpsig(5))*Res23s12 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11s12 - aux3*((auxG/acvis)*gpsig(6))*Res11s12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22s12 - aux3*((auxG/acvis)*gpsig(6))*Res22s12 &              
                  + elmat(6,inode,6,jnode)
               
         end do
      end do  

   end subroutine supf_elmbstEst3d   
   
   subroutine supf_elmrhcGal(e,auxVE,beta,auxG,auxPTT,acvis,dvolu,auxtens,elextS,gpvel,grvel,gpsig,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: elextS(auxtens),gpvel(e%ndime),grvel(e%ndime,e%ndime),beta
      real(rp),    intent(in)    :: dvolu,acvis,auxG,auxVE,gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux1
      real(rp)                   :: Res11,Res12,Res22
      
      Res11 = elextS(1) + (1.0_rp-beta)*grvel(1,1) + (auxG/acvis)*(gpsig(1)*gpsig(1) + gpsig(3)*gpsig(3))*auxVE
      Res12 = elextS(3) + 0.5_rp*(1.0_rp-beta)*(grvel(1,2) + grvel(2,1)) &
         + (auxG/acvis)*(gpsig(1) + gpsig(2))*gpsig(3)*auxVE
      Res22 = elextS(2) + (1.0_rp-beta)*grvel(2,2) + (auxG/acvis)*(gpsig(3)*gpsig(3) + gpsig(2)*gpsig(2))*auxVE         
     
      do inode=1,e%pnode
      
         aux1=e%shape(inode,e%igaus)*dvolu
         
         elrhs(1,inode) = aux1*Res11 + elrhs(1,inode) 
         elrhs(2,inode) = aux1*Res22 + elrhs(2,inode)
         elrhs(3,inode) = aux1*(2.0_rp*Res12) + elrhs(3,inode)
  
      end do

   end subroutine supf_elmrhcGal
   
   subroutine supf_elmrhcGal3d(e,auxVE,beta,auxG,auxPTT,acvis,dvolu,auxtens,elextS,gpvel,grvel,gpsig,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT
      real(rp),    intent(in)    :: elextS(auxtens),gpvel(e%ndime),grvel(e%ndime,e%ndime),beta
      real(rp),    intent(in)    :: dvolu,acvis,auxG,auxVE,gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux1
      real(rp)                   :: Res11,Res12,Res13,Res22,Res23,Res33
      
      Res11 = elextS(1) + (auxVE)*((auxG/acvis)*(gpsig(1)*gpsig(1) + gpsig(6)*gpsig(6) + gpsig(5)*gpsig(5))) &
         + (1.0_rp-beta)*grvel(1,1) 
      Res22 = elextS(2) + (auxVE)*((auxG/acvis)*(gpsig(2)*gpsig(2) + gpsig(6)*gpsig(6) + gpsig(4)*gpsig(4))) &
         + (1.0_rp-beta)*grvel(2,2)
      Res33 = elextS(3) + (auxVE)*((auxG/acvis)*(gpsig(5)*gpsig(5) + gpsig(4)*gpsig(4) + gpsig(3)*gpsig(3))) &
         + (1.0_rp-beta)*grvel(3,3)
      Res12 = elextS(6) + (auxVE)*((auxG/acvis)*((gpsig(1) + gpsig(2))*gpsig(6) + gpsig(4)*gpsig(5))) &
         + 0.5_rp*(1.0_rp-beta)*(grvel(1,2) + grvel(2,1))
      Res13 = elextS(5) + (auxVE)*((auxG/acvis)*((gpsig(1) + gpsig(3))*gpsig(5) + gpsig(4)*gpsig(6))) &
         + 0.5_rp*(1.0_rp-beta)*(grvel(1,3) + grvel(3,1))
      Res23 = elextS(4) + (auxVE)*((auxG/acvis)*((gpsig(2) + gpsig(3))*gpsig(4) + gpsig(6)*gpsig(5))) &
         + 0.5_rp*(1.0_rp-beta)*(grvel(2,3) + grvel(3,2))

      do inode=1,e%pnode
      
         aux1=e%shape(inode,e%igaus)*dvolu
         
         elrhs(1,inode) = aux1*Res11 + elrhs(1,inode) 
         elrhs(2,inode) = aux1*Res22 + elrhs(2,inode)
         elrhs(3,inode) = aux1*Res33 + elrhs(3,inode) 
         elrhs(4,inode) = aux1*(2.0_rp*Res23) + elrhs(4,inode)         
         elrhs(5,inode) = aux1*(2.0_rp*Res13) + elrhs(5,inode)         
         elrhs(6,inode) = aux1*(2.0_rp*Res12) + elrhs(6,inode)
  
      end do

   end subroutine supf_elmrhcGal3d  
   
   subroutine supf_elmrhc_splitoss(e,auxtens,beta,timom,dvolu,gpdivs,elrhs)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the rhs terms in constitutive equation for OSS
   !    
   ! -tau1(div(tau),resid_momentum)- tau3/2mu(tau,resid_constitutive)
   !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: gpdivs(e%ndime)
      real(rp),    intent(in)    :: dvolu,beta,timom
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,aux2
      real(rp)                   :: aux1,aux3

      aux1 = (1.0_rp-beta)*timom*dvolu  
      do inode=1,e%pnode

         elrhs(1,inode) = e%cartd(1,inode)*aux1*gpdivs(1) + elrhs(1,inode)
         elrhs(2,inode) = e%cartd(2,inode)*aux1*gpdivs(2) + elrhs(2,inode)
         elrhs(3,inode) = aux1*(e%cartd(2,inode)*gpdivs(1) + e%cartd(1,inode)*gpdivs(2)) + elrhs(3,inode)        

      end do
   end subroutine supf_elmrhc_splitoss  
   
   subroutine supf_elmrhc_splitoss3d(e,auxtens,beta,timom,dvolu,gpdivs,elrhs)
   !-----------------------------------------------------------------------
   !
   ! This routine computes the rhs terms in constitutive equation for OSS
   !    
   ! -tau1(div(tau),resid_momentum)- tau3/2mu(tau,resid_constitutive)
   !-----------------------------------------------------------------------
      use typre
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: gpdivs(e%ndime)
      real(rp),    intent(in)    :: dvolu,beta,timom
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,aux2
      real(rp)                   :: aux1,aux3

      aux1 = (1.0_rp-beta)*timom*dvolu   
      do inode=1,e%pnode

         elrhs(1,inode) = e%cartd(1,inode)*aux1*gpdivs(1) + elrhs(1,inode)
         elrhs(2,inode) = e%cartd(2,inode)*aux1*gpdivs(2) + elrhs(2,inode)
         elrhs(3,inode) = e%cartd(3,inode)*aux1*gpdivs(3) + elrhs(3,inode)
         elrhs(4,inode) = aux1*(e%cartd(3,inode)*gpdivs(2) + e%cartd(2,inode)*gpdivs(3)) + elrhs(4,inode)            
         elrhs(5,inode) = aux1*(e%cartd(3,inode)*gpdivs(1) + e%cartd(1,inode)*gpdivs(3)) + elrhs(5,inode)            
         elrhs(6,inode) = aux1*(e%cartd(2,inode)*gpdivs(1) + e%cartd(1,inode)*gpdivs(2)) + elrhs(6,inode)        

      end do
   end subroutine supf_elmrhc_splitoss3d  
  
   subroutine supf_elmrhcVES_oss(e,auxGie,auxVE,auxG,auxPTT,tisig,acvis,dvolu,auxtens,gpadv,grvel,gpsig,gprep,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxGie
      real(rp),    intent(in)    :: tisig,gprep(e%ndime+1+auxtens),gpadv(e%pnode),grvel(e%ndime,e%ndime) 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux,aux2,aux3,aux5
      real(rp)                   :: Res11,Res12,Res22
         
         
      Res11 = gprep(1)
      Res12 = gprep(3)
      Res22 = gprep(2)

      do inode=1,e%pnode

     
         aux2= (e%shape(inode,e%igaus)*(-tisig/(2.0_rp*acvis)) + (tisig*auxVE)*gpadv(inode))*dvolu          
         aux5= (tisig*auxVE)*(e%shape(inode,e%igaus))*dvolu
         
              
         elrhs(1,inode) = aux2*Res11 &
               + aux5*(2.0_rp*grvel(1,1) - (auxGie*auxG/acvis)*gpsig(1))*(Res11) &              
               + aux5*(2.0_rp*grvel(1,2) - (auxGie*auxG/acvis)*gpsig(3))*(Res12) &               
               + elrhs(1,inode) 
               
         elrhs(2,inode) = aux2*Res22 &
               + aux5*(2.0_rp*grvel(2,2) - (auxGie*auxG/acvis)*gpsig(2))*(Res22) &                
               + aux5*(2.0_rp*grvel(2,1) - (auxGie*auxG/acvis)*gpsig(3))*(Res12) &                    
               + elrhs(2,inode)
               
         elrhs(3,inode) = aux2*(2.0_rp*Res12)  &
               + aux5*(2.0_rp*grvel(2,1) - (auxGie*auxG/acvis)*gpsig(3))*(Res11) &
               + aux5*(2.0_rp*grvel(1,2) - (auxGie*auxG/acvis)*gpsig(3))*(Res22) &
               - aux5*(auxGie*auxG/acvis)*(gpsig(1) + gpsig(2))*(Res12) &                                   
               + elrhs(3,inode)              
  
      end do

   end subroutine supf_elmrhcVES_oss   

   subroutine supf_elmrhcVES_oss3d(e,auxGie,auxVE,auxG,auxPTT,tisig,acvis,dvolu,auxtens,gpadv,grvel,gpsig,gprep,elrhs)   
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxGie
      real(rp),    intent(in)    :: tisig,gprep(e%ndime+1+auxtens),gpadv(e%pnode),grvel(e%ndime,e%ndime) 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpsig(auxtens)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux,aux2,aux3,aux5
      real(rp)                   :: Res11,Res12,Res13,Res22,Res23,Res33
      
      
      
      Res11 = gprep(1)            
      Res12 = gprep(6)
      Res13 = gprep(5)
      Res22 = gprep(2)  
      Res23 = gprep(4)
      Res33 = gprep(3)   
 
      do inode=1,e%pnode 
      
         aux2= (e%shape(inode,e%igaus)*(-tisig/(2.0_rp*acvis)) + (tisig*auxVE)*gpadv(inode))*dvolu          
         aux3= (tisig*auxVE)*(e%shape(inode,e%igaus))*dvolu                  
               
 
         elrhs(1,inode) = aux2*Res11 &
                  + aux3*(2.0_rp*grvel(1,1))*Res11 - aux3*((auxGie*auxG/acvis)*gpsig(1))*Res11 &
                  + aux3*(2.0_rp*grvel(1,2))*Res12 - aux3*((auxGie*auxG/acvis)*gpsig(6))*Res12 &
                  + aux3*(2.0_rp*grvel(1,3))*Res13 - aux3*((auxGie*auxG/acvis)*gpsig(5))*Res13 &   
                  + elrhs(1,inode) 
               
         elrhs(2,inode) =  aux2*Res22 &
                  + aux3*(2.0_rp*grvel(2,2))*Res22 - aux3*((auxGie*auxG/acvis)*gpsig(2))*Res22 &           
                  + aux3*(2.0_rp*grvel(2,1))*Res12 - aux3*((auxGie*auxG/acvis)*gpsig(6))*Res12 & 
                  + aux3*(2.0_rp*grvel(2,3))*Res23 - aux3*((auxGie*auxG/acvis)*gpsig(4))*Res23 &                
                  + elrhs(2,inode)
               
         elrhs(3,inode) =  aux2*Res33 &
                  + aux3*(2.0_rp*grvel(3,3))*Res33 - aux3*((auxGie*auxG/acvis)*gpsig(3))*Res33 &           
                  + aux3*(2.0_rp*grvel(3,1))*Res13 - aux3*((auxGie*auxG/acvis)*gpsig(5))*Res13 & 
                  + aux3*(2.0_rp*grvel(3,2))*Res23 - aux3*((auxGie*auxG/acvis)*gpsig(4))*Res23 &                   
                  + elrhs(3,inode)
                  
         elrhs(4,inode) = (2.0_rp*aux2)*Res23 &
                  + aux3*(2.0_rp*(grvel(3,3) + grvel(2,2)))*Res23 - aux3*((auxGie*auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23 &
                  + aux3*(2.0_rp*grvel(3,1))*Res12 - aux3*((auxGie*auxG/acvis)*gpsig(5))*Res12 &           
                  + aux3*(2.0_rp*grvel(2,1))*Res13 - aux3*((auxGie*auxG/acvis)*gpsig(6))*Res13 & 
                  + aux3*(2.0_rp*grvel(3,2))*Res22 - aux3*((auxGie*auxG/acvis)*gpsig(4))*Res22 &           
                  + aux3*(2.0_rp*grvel(2,3))*Res33 - aux3*((auxGie*auxG/acvis)*gpsig(4))*Res33 &              
                  + elrhs(4,inode)
                  
         elrhs(5,inode) = (2.0_rp*aux2)*Res13 &
                  + aux3*(2.0_rp*(grvel(3,3) + grvel(1,1)))*Res13 - aux3*((auxGie*auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13 &
                  + aux3*(2.0_rp*grvel(3,2))*Res12 - aux3*((auxGie*auxG/acvis)*gpsig(4))*Res12 &           
                  + aux3*(2.0_rp*grvel(1,2))*Res23 - aux3*((auxGie*auxG/acvis)*gpsig(6))*Res23 & 
                  + aux3*(2.0_rp*grvel(3,1))*Res11 - aux3*((auxGie*auxG/acvis)*gpsig(5))*Res11 &           
                  + aux3*(2.0_rp*grvel(1,3))*Res33 - aux3*((auxGie*auxG/acvis)*gpsig(5))*Res33 &            
                  + elrhs(5,inode)
                 
         elrhs(6,inode) = (2.0_rp*aux2)*Res12 &
                  + aux3*(2.0_rp*(grvel(2,2) + grvel(1,1)))*Res12 - aux3*((auxGie*auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12 &
                  + aux3*(2.0_rp*grvel(2,3))*Res13 - aux3*((auxGie*auxG/acvis)*gpsig(4))*Res13 &           
                  + aux3*(2.0_rp*grvel(1,3))*Res23 - aux3*((auxGie*auxG/acvis)*gpsig(5))*Res23 & 
                  + aux3*(2.0_rp*grvel(2,1))*Res11 - aux3*((auxGie*auxG/acvis)*gpsig(6))*Res11 &           
                  + aux3*(2.0_rp*grvel(1,2))*Res22 - aux3*((auxGie*auxG/acvis)*gpsig(6))*Res22 &          
                  + elrhs(6,inode) 
  
      end do

   end subroutine supf_elmrhcVES_oss3d    

   
   subroutine supf_elmrhcEst(e,auxGie,auxconv,auxtrac1,auxtrac2,auxVE,beta,auxG,auxPTT,tisig,acvis,gpsig,gpvel,gpadv,grvel &
         ,dvolu,auxtens,elextS,elrhs)         
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxconv,auxtrac1,auxtrac2,auxGie
      real(rp),    intent(in)    :: gpsig(auxtens),tisig,elextS(auxtens),beta 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpvel(e%ndime),gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux2,aux5
      real(rp)                   :: Res11,Res12,Res22
      
      
!       Res11 = elextS(1) + (auxG/acvis)*(gpsig(1)*gpsig(1) + gpsig(3)*gpsig(3))*auxVE
!       Res12 = elextS(3) + (auxG/acvis)*(gpsig(1) + gpsig(2))*gpsig(3)*auxVE
!       Res22 = elextS(2) + (auxG/acvis)*(gpsig(3)*gpsig(3) + gpsig(2)*gpsig(2))*auxVE                   
                     
      Res11 = elextS(1) + (1.0_rp-beta)*grvel(1,1) + (auxVE)*((auxG/acvis)*(gpsig(1)*gpsig(1) + gpsig(3)*gpsig(3)))
      Res12 = elextS(3) + 0.5_rp*(1.0_rp-beta)*(grvel(1,2) + grvel(2,1)) &
         + (auxVE)*((auxG/acvis)*(gpsig(1) + gpsig(2))*gpsig(3)) 
      Res22 = elextS(2) + (1.0_rp-beta)*grvel(2,2) + (auxVE)*((auxG/acvis)*(gpsig(3)*gpsig(3) + gpsig(2)*gpsig(2)))

      do inode=1,e%pnode      

         aux2=(e%shape(inode,e%igaus)*(-(tisig/(2.0_rp*acvis)))  + (tisig*auxVE)*(gpadv(inode)))*dvolu  
         aux5=(tisig*auxVE)*(e%shape(inode,e%igaus))*dvolu
         
         elrhs(1,inode) = aux2*Res11 &
               + aux5*(auxtrac2*2.0_rp*grvel(1,1) - (auxGie*auxG/acvis)*gpsig(1))*(Res11) &              
               + aux5*(auxtrac2*2.0_rp*grvel(1,2) - (auxGie*auxG/acvis)*gpsig(3))*(Res12) &              
               + elrhs(1,inode) 
               
         elrhs(2,inode) = aux2*Res22 &
               + aux5*(auxtrac2*2.0_rp*grvel(2,2) - (auxGie*auxG/acvis)*gpsig(2))*(Res22) &                 
               + aux5*(auxtrac2*2.0_rp*grvel(2,1) - (auxGie*auxG/acvis)*gpsig(3))*(Res12) &                   
               + elrhs(2,inode)
               
         elrhs(3,inode) = aux2*(2.0_rp*Res12)  &
               + aux5*(auxtrac2*2.0_rp*grvel(2,1) - (auxGie*auxG/acvis)*gpsig(3))*(Res11) &
               + aux5*(auxtrac2*2.0_rp*grvel(1,2) - (auxGie*auxG/acvis)*gpsig(3))*(Res22) &
               + aux5*(-(auxGie*auxG/acvis)*gpsig(1))*(Res12) &
               + aux5*(-(auxGie*auxG/acvis)*gpsig(2))*(Res12) &                 
               + elrhs(3,inode)               
 
      end do

   end subroutine supf_elmrhcEst   

   
   subroutine supf_elmrhcEst3d(e,auxconv,auxtrac1,auxtrac2,auxVE,beta,auxG,auxPTT,tisig,acvis,gpsig,gpvel,gpadv,grvel &
         ,dvolu,auxtens,elextS,elrhs)         
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      integer(ip), intent(in)    :: auxtens,auxPTT,auxconv,auxtrac1,auxtrac2
      real(rp),    intent(in)    :: gpsig(auxtens),tisig,elextS(auxtens),beta 
      real(rp),    intent(in)    :: dvolu,auxVE,auxG,acvis,gpvel(e%ndime),gpadv(e%pnode),grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode
      real(rp)                   :: aux2,aux3
      real(rp)                   :: Res11,Res12,Res22,Res13,Res23,Res33
      
      
      Res11 = elextS(1) + (auxVE)*((auxG/acvis)*(gpsig(1)*gpsig(1) + gpsig(6)*gpsig(6) + gpsig(5)*gpsig(5))) &
         + (1.0_rp-beta)*grvel(1,1) 
      Res22 = elextS(2) + (auxVE)*((auxG/acvis)*(gpsig(2)*gpsig(2) + gpsig(6)*gpsig(6) + gpsig(4)*gpsig(4))) &
         + (1.0_rp-beta)*grvel(2,2)
      Res33 = elextS(3) + (auxVE)*((auxG/acvis)*(gpsig(5)*gpsig(5) + gpsig(4)*gpsig(4) + gpsig(3)*gpsig(3))) &
         + (1.0_rp-beta)*grvel(3,3)
      Res12 = elextS(6) + (auxVE)*((auxG/acvis)*((gpsig(1) + gpsig(2))*gpsig(6) + gpsig(4)*gpsig(5))) &
         + 0.5_rp*(1.0_rp-beta)*(grvel(1,2) + grvel(2,1))
      Res13 = elextS(5) + (auxVE)*((auxG/acvis)*((gpsig(1) + gpsig(3))*gpsig(5) + gpsig(4)*gpsig(6))) &
         + 0.5_rp*(1.0_rp-beta)*(grvel(1,3) + grvel(3,1))
      Res23 = elextS(4) + (auxVE)*((auxG/acvis)*((gpsig(2) + gpsig(3))*gpsig(4) + gpsig(6)*gpsig(5))) &
         + 0.5_rp*(1.0_rp-beta)*(grvel(2,3) + grvel(3,2))                 
                     

      do inode=1,e%pnode      

         aux2= (e%shape(inode,e%igaus)*(-(tisig/(2.0_rp*acvis))) + (tisig*auxVE)*gpadv(inode))*dvolu  
         aux3= ((tisig*auxVE)*e%shape(inode,e%igaus))*dvolu      
               
         elrhs(1,inode) = aux2*Res11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,1))*Res11 - aux3*((auxG/acvis)*gpsig(1))*Res11 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res12 - aux3*((auxG/acvis)*gpsig(6))*Res12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res13 - aux3*((auxG/acvis)*gpsig(5))*Res13 &   
                  + elrhs(1,inode) 
               
         elrhs(2,inode) =  aux2*Res22 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,2))*Res22 - aux3*((auxG/acvis)*gpsig(2))*Res22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res12 - aux3*((auxG/acvis)*gpsig(6))*Res12 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res23 - aux3*((auxG/acvis)*gpsig(4))*Res23 &                
                  + elrhs(2,inode)
               
         elrhs(3,inode) =  aux2*Res33 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,3))*Res33 - aux3*((auxG/acvis)*gpsig(3))*Res33 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res13 - aux3*((auxG/acvis)*gpsig(5))*Res13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res23 - aux3*((auxG/acvis)*gpsig(4))*Res23 &                   
                  + elrhs(3,inode)
                  
         elrhs(4,inode) = (2.0_rp*aux2)*Res23 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(2,2)))*Res23 - aux3*((auxG/acvis)*(gpsig(2)+gpsig(3)))*Res23 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res12 - aux3*((auxG/acvis)*gpsig(5))*Res12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res13 - aux3*((auxG/acvis)*gpsig(6))*Res13 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res22 - aux3*((auxG/acvis)*gpsig(4))*Res22 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res33 - aux3*((auxG/acvis)*gpsig(4))*Res33 &              
                  + elrhs(4,inode)
                  
         elrhs(5,inode) = (2.0_rp*aux2)*Res13 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(3,3) + grvel(1,1)))*Res13 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(3)))*Res13 &
                  + aux3*(2.0_rp*auxtrac2*grvel(3,2))*Res12 - aux3*((auxG/acvis)*gpsig(4))*Res12 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res23 - aux3*((auxG/acvis)*gpsig(6))*Res23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(3,1))*Res11 - aux3*((auxG/acvis)*gpsig(5))*Res11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res33 - aux3*((auxG/acvis)*gpsig(5))*Res33 &            
                  + elrhs(5,inode)
                 
         elrhs(6,inode) = (2.0_rp*aux2)*Res12 &
                  + aux3*(2.0_rp*auxtrac2*(grvel(2,2) + grvel(1,1)))*Res12 - aux3*((auxG/acvis)*(gpsig(1)+gpsig(2)))*Res12 &
                  + aux3*(2.0_rp*auxtrac2*grvel(2,3))*Res13 - aux3*((auxG/acvis)*gpsig(4))*Res13 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,3))*Res23 - aux3*((auxG/acvis)*gpsig(5))*Res23 & 
                  + aux3*(2.0_rp*auxtrac2*grvel(2,1))*Res11 - aux3*((auxG/acvis)*gpsig(6))*Res11 &           
                  + aux3*(2.0_rp*auxtrac2*grvel(1,2))*Res22 - aux3*((auxG/acvis)*gpsig(6))*Res22 &          
                  + elrhs(6,inode)                 
  
      end do

   end subroutine supf_elmrhcEst3d      
   
   
   !---------------------------------------------------------------------
   !Third Step "First Approach" Split continuity equation
   !---------------------------------------------------------------------    

   subroutine supf_embpq_step3(e,dvolu,acden,timom,dtinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for p_n+1 for Fractional step
      ! -(div u, q)-tau*(u · grad u, grad q)+tau*(proj, grad q)+(1/rho)*dt*(grad p_n, grad q)
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: acden
      real(rp),    intent(in)    :: dvolu,timom,dtinv
      real(rp),    intent(inout) :: elmat(1,e%mnode,1,e%mnode)
      integer(ip)                :: inode,jnode,idime
      real(rp)                   :: aux1,aux2

     
      aux1 = dtinv*acden
      aux2= -(timom + (1.0_rp/aux1))*dvolu

      do inode=1,e%pnode
         do jnode=1,e%pnode
            do idime=1,e%ndime
               elmat(1,inode,1,jnode) = (e%cartd(idime,inode)*e%cartd(idime,jnode))*aux2 &
                     + elmat(1,inode,1,jnode)
            end do
         end do
      end do

   end subroutine supf_embpq_step3    
   
   subroutine supf_elmrhp(e,dvolu,acden,timom,dtinv,divvel,divsig1,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for p_n+1 for Fractional step
      ! -(div u, q)-tau*(u · grad u, grad q)+tau*(proj, grad q)+(1/rho)*dt*(grad p_n, grad q)
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: acden
      real(rp),    intent(in)    :: dvolu,timom,dtinv,divvel,divsig1(e%ndime)
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: tmp,tmp2
      real(rp)                   :: aux1

      tmp = timom*dvolu      
      aux1 = dtinv*acden
      tmp2= dvolu*(1.0_rp/aux1)

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = -(tmp2*divsig1(idime))*e%cartd(idime,inode) &
                  + elrhs(1,inode)
         end do         
         elrhs(1,inode) = dvolu*divvel*e%shape(inode,e%igaus) + elrhs(1,inode)
      end do

   end subroutine supf_elmrhp 
   
   subroutine supf_elmrhp_splitoss(e,timom,dvolu,gpgrap,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for OSS
    !    tau1*(grad q, resid_continuity)
    !
    !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e      
      real(rp),    intent(in)    :: gpgrap(e%ndime),timom,dvolu
      real(rp),    intent(inout) :: elrhs(1,e%mnode)
      integer(ip)                :: inode,idime,aux1
      real(rp)                   :: tmp1


      tmp1 = -(dvolu*timom)

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = e%cartd(idime,inode)*gpgrap(idime)*tmp1 + elrhs(1,inode)
         end do
      end do

   end subroutine supf_elmrhp_splitoss    
   
   subroutine supf_elmrhpextra(e,dvolu,acden,dtinv,divsig2,grpre,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for p_n+1 for Fractional step
      ! -(div u, q)-tau*(u · grad u, grad q)+tau*(proj, grad q)+(1/rho)*dt*(grad p_n, grad q)
      !
      !-----------------------------------------------------------------------
      use typre
      implicit none

      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: acden
      real(rp),    intent(in)    :: grpre(e%ndime)
      real(rp),    intent(in)    :: dvolu,dtinv,divsig2(e%ndime)
      real(rp),    intent(inout) :: elrhs(1,e%mnode)

      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: tmp,tmp2
      real(rp)                   :: aux1
      
      aux1 = dtinv*acden
      tmp2= dvolu*(1.0_rp/aux1)

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(1,inode) = (tmp2*(-grpre(idime) + divsig2(idime)))*e%cartd(idime,inode) &
                  + elrhs(1,inode)
         end do         
      end do

   end subroutine supf_elmrhpextra     
   
   
   subroutine supf_elmdir_press(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement) :: e
      class(SUPFractionalStepProblem) :: a
      real(rp) :: elmat(1,e%mnode,1,e%mnode), elrhs(1,e%mnode)
      
      real(rp) :: adiag
      integer(ip) :: inode,ipoin,jnode
      real(rp)    :: prepr
      
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         
         if (ipoin == a%nodpr .or. a%kfl_fixpr(ipoin) /= 0) then
            if (a%kfl_fixpr(ipoin) == 2) then
               prepr = a%bvpress(ipoin)
            else
               prepr = a%prepr
            endif
            adiag=elmat(1,inode,1,inode)
            elmat(1,inode,1,1:e%pnode)=0.0_rp
            do jnode = 1,e%pnode
               elrhs(1,jnode) = elrhs(1,jnode)- elmat(1,jnode,1,inode)*prepr
            enddo
            elmat(1,1:e%pnode,1,inode) = 0.0_rp
            elmat(1,inode,1,inode) = adiag
            elrhs(1,inode) = adiag*prepr
         endif
      enddo
   end subroutine 
   
   !---------------------------------------------------------------------
   !Fourth Step "First Approach" Correction momentum equation
   !---------------------------------------------------------------------    
   
   subroutine supf_elmrhu_corr(e,dvolu,denac,dtinv,gpveln,gppre,gpsig,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for Fractional step : end-of-step u
      !    (v, u_int/dt) + (p_n+1-p_n, div v)
      !
      !-----------------------------------------------------------------------
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: gpveln(e%ndime),gpsig((e%ndime-1)*(e%ndime-1)+2)
      real(rp),    intent(in)    :: dvolu,dtinv,gppre,denac
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)

      integer(ip)                :: inode,idime
      real(rp)                   :: aux


      
      if(e%ndime==2)then
         do inode=1,e%pnode
            aux = e%shape(inode,e%igaus)*dvolu*denac*dtinv
            elrhs(1,inode) = (gppre*e%cartd(1,inode)-(gpsig(1)*e%cartd(1,inode) + gpsig(3)*e%cartd(2,inode)))*dvolu &
                  + gpveln(1)*aux + elrhs(1,inode)
            elrhs(2,inode) = (gppre*e%cartd(2,inode)-(gpsig(3)*e%cartd(1,inode) + gpsig(2)*e%cartd(2,inode)))*dvolu &
                  + gpveln(2)*aux + elrhs(2,inode)
         end do
      elseif(e%ndime==3)then
         do inode=1,e%pnode
            aux = e%shape(inode,e%igaus)*dvolu*denac*dtinv
            elrhs(1,inode) = (gppre*e%cartd(1,inode)-(gpsig(1)*e%cartd(1,inode)+gpsig(6)*e%cartd(2,inode)+gpsig(5)*e%cartd(3,inode)))*dvolu &
               + gpveln(1)*aux + elrhs(1,inode)
            elrhs(2,inode) = (gppre*e%cartd(2,inode)-(gpsig(6)*e%cartd(1,inode)+gpsig(2)*e%cartd(2,inode)+gpsig(4)*e%cartd(3,inode)))*dvolu &
               + gpveln(2)*aux + elrhs(2,inode)
            elrhs(3,inode) = (gppre*e%cartd(3,inode)-(gpsig(5)*e%cartd(1,inode)+gpsig(4)*e%cartd(2,inode)+gpsig(3)*e%cartd(3,inode)))*dvolu &
               + gpveln(3)*aux + elrhs(3,inode) 
         end do
      end if
      
   end subroutine supf_elmrhu_corr
   
   subroutine supf_elmrhu_corr_extra(e,dvolu,gppre,gpsig,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for Fractional step : end-of-step u
      !    (v, u_int/dt) + (p_n+1-p_n, div v)
      !
      !-----------------------------------------------------------------------
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: gpsig((e%ndime-1)*(e%ndime-1)+2)
      real(rp),    intent(in)    :: dvolu,gppre
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime
      real(rp)                   :: aux

      if(e%ndime==2)then
         do inode=1,e%pnode
            elrhs(1,inode) = (-gppre*e%cartd(1,inode) + gpsig(1)*e%cartd(1,inode) + gpsig(3)*e%cartd(2,inode))*dvolu &
               + elrhs(1,inode)
            elrhs(2,inode) = (-gppre*e%cartd(2,inode) + gpsig(3)*e%cartd(1,inode) + gpsig(2)*e%cartd(2,inode))*dvolu + elrhs(2,inode)
         end do
      elseif(e%ndime==3)then
         do inode=1,e%pnode
            elrhs(1,inode) = (-gppre*e%cartd(1,inode)+gpsig(1)*e%cartd(1,inode)+gpsig(6)*e%cartd(2,inode)+gpsig(5)*e%cartd(3,inode))*dvolu &
               + elrhs(1,inode)
            elrhs(2,inode) = (-gppre*e%cartd(2,inode)+gpsig(6)*e%cartd(1,inode)+gpsig(2)*e%cartd(2,inode)+gpsig(4)*e%cartd(3,inode))*dvolu &
               + elrhs(2,inode)
            elrhs(3,inode) = (-gppre*e%cartd(3,inode)+gpsig(5)*e%cartd(1,inode)+gpsig(4)*e%cartd(2,inode)+gpsig(3)*e%cartd(3,inode))*dvolu &
               + elrhs(3,inode) 
         end do
      end if
      
   end subroutine supf_elmrhu_corr_extra 
   
   
   subroutine supf_elmrhuY(e,acden,timom,tisig,beta,acvis,dvolu,gpvel,gpvel2,grvel,grvel2,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: acden,gpadv(e%pnode),grvel2(e%ndime,e%ndime),tisig
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),gpvel2(e%ndime),grvel(e%ndime,e%ndime),timom,beta,acvis
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2,aux1,aux3


      do inode=1,e%pnode
         
         aux1 = (2.0_rp*acvis*beta + (1.0_rp - beta)*tisig)*dvolu   
           
           
         elrhs(1,inode)= aux1*(grvel2(1,1)*e%cartd(1,inode) + 0.5_rp*(grvel2(1,2) + grvel2(2,1))*e%cartd(2,inode)) &        
            + elrhs(1,inode)
            
         elrhs(2,inode)= aux1*(grvel2(2,2)*e%cartd(2,inode) + 0.5_rp*(grvel2(1,2) + grvel2(2,1))*e%cartd(1,inode)) &                    
            + elrhs(2,inode)              
           
            
    
       end do


   end subroutine supf_elmrhuY  
   
   
   subroutine supf_elmrhuYConv(e,acden,timom,tisig,beta,acvis,dvolu,gpvel,gpvel2,grvel,grvel2,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: acden,gpadv(e%pnode),grvel2(e%ndime,e%ndime),tisig
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),gpvel2(e%ndime),grvel(e%ndime,e%ndime),timom,beta,acvis
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2,aux1,aux3


      do inode=1,e%pnode
         
         aux  = acden*e%shape(inode,e%igaus)*dvolu
         aux2= (timom*acden*gpadv(inode))*dvolu
         aux3= timom*acden*dvolu
         
         elrhs(1,inode)= aux*(gpvel2(1)*grvel2(1,1) + gpvel2(2)*grvel2(1,2))  &     
            + (aux2 + aux)*acden*(gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2)) &
            + aux3*(gpvel2(1)*e%cartd(1,inode) + gpvel2(2)*e%cartd(2,inode))*(gpvel2(1)*grvel2(1,1) + gpvel2(2)*grvel2(1,2)) &
            + elrhs(1,inode)
            
         elrhs(2,inode)= aux*(gpvel2(1)*grvel2(2,1) + gpvel2(2)*grvel2(2,2))  &      
            + (aux2 + aux)*acden*(gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2)) &
            + aux3*(gpvel2(1)*e%cartd(1,inode) + gpvel2(2)*e%cartd(2,inode))*(gpvel2(1)*grvel2(2,1) + gpvel2(2)*grvel2(2,2)) &            
            + elrhs(2,inode)       
            
    
       end do


   end subroutine supf_elmrhuYConv    
   
   
   
   subroutine supf_elmrhuY3d(e,acden,timom,tisig,beta,acvis,dvolu,gpvel,gpvel2,grvel,grvel2,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: acden,gpadv(e%pnode),grvel2(e%ndime,e%ndime),tisig
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),gpvel2(e%ndime),grvel(e%ndime,e%ndime),timom,beta,acvis
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2,aux1,aux3


      do inode=1,e%pnode         

         aux1 = (2.0_rp*acvis*beta + (1.0_rp-beta)*tisig)*dvolu
         
         elrhs(1,inode)= aux1*(grvel2(1,1)*e%cartd(1,inode) + 0.5_rp*(grvel2(1,2) + grvel2(2,1))*e%cartd(2,inode) &
            + 0.5_rp*(grvel2(1,3) + grvel2(3,1))*e%cartd(3,inode)) + elrhs(1,inode)
            
         elrhs(2,inode)= aux1*(grvel2(2,2)*e%cartd(2,inode) + 0.5_rp*(grvel2(1,2) + grvel2(2,1))*e%cartd(1,inode) &
            + 0.5_rp*(grvel2(2,3) + grvel2(3,2))*e%cartd(3,inode)) + elrhs(2,inode)
            
         elrhs(3,inode)= aux1*(grvel2(3,3)*e%cartd(3,inode) + 0.5_rp*(grvel2(1,3) + grvel2(3,1))*e%cartd(1,inode) &
            + 0.5_rp*(grvel2(2,3) + grvel2(3,2))*e%cartd(2,inode)) + elrhs(3,inode)            
            
    
       end do


   end subroutine supf_elmrhuY3d     
   
   
   subroutine supf_elmrhuYConv3d(e,acden,timom,tisig,beta,acvis,dvolu,gpvel,gpvel2,grvel,grvel2,gpadv,elrhs)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS in constitutive equation
    !    -tau1*(div(T), f) -tau1(div(T),rho*u/dt)
    !
    !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e

      real(rp),    intent(in)    :: acden,gpadv(e%pnode),grvel2(e%ndime,e%ndime),tisig
      real(rp),    intent(in)    :: dvolu,gpvel(e%ndime),gpvel2(e%ndime),grvel(e%ndime,e%ndime),timom,beta,acvis
      real(rp),    intent(inout) :: elrhs(e%ndime,e%mnode)
      integer(ip)                :: inode,idime,jdime
      real(rp)                   :: aux,aux2,aux1,aux3


      do inode=1,e%pnode
         
         aux  = acden*e%shape(inode,e%igaus)*dvolu
         aux2= (timom*acden*gpadv(inode))*dvolu
         aux3= timom*dvolu*acden
         
         elrhs(1,inode)= aux*(gpvel2(1)*grvel2(1,1) + gpvel2(2)*grvel2(1,2) + gpvel2(3)*grvel2(1,3))  &         
            + (aux2 + aux)*acden*(gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2) + gpvel(3)*grvel(1,3)) &
            + aux3*(gpvel2(1)*e%cartd(1,inode) + gpvel2(2)*e%cartd(2,inode) + gpvel2(3)*e%cartd(3,inode))* &
               (gpvel2(1)*grvel2(1,1) + gpvel2(2)*grvel2(1,2) + gpvel2(3)*grvel2(1,3)) &           
            + elrhs(1,inode)
            
         elrhs(2,inode)= aux*(gpvel2(1)*grvel2(2,1) + gpvel2(2)*grvel2(2,2) + gpvel2(3)*grvel2(2,3))  &         
            + (aux2 + aux)*acden*(gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2) + gpvel(3)*grvel(2,3)) &
            + aux3*(gpvel2(1)*e%cartd(1,inode) + gpvel2(2)*e%cartd(2,inode) + gpvel2(3)*e%cartd(3,inode))* &
               (gpvel2(1)*grvel2(2,1) + gpvel2(2)*grvel2(2,2) + gpvel2(3)*grvel2(2,3)) &               
            + elrhs(2,inode)
            
         elrhs(3,inode)= aux*(gpvel2(1)*grvel2(3,1) + gpvel2(2)*grvel2(3,2) + gpvel2(3)*grvel2(3,3))  &         
            + (aux2 + aux)*acden*(gpvel(1)*grvel(3,1) + gpvel(2)*grvel(3,2) + gpvel(3)*grvel(3,3)) &
            + aux3*(gpvel2(1)*e%cartd(1,inode) + gpvel2(2)*e%cartd(2,inode) + gpvel2(3)*e%cartd(3,inode))* &
               (gpvel2(1)*grvel2(3,1) + gpvel2(2)*grvel2(3,2) + gpvel2(3)*grvel2(3,3)) &             
            + elrhs(3,inode)            
            
    
       end do


   end subroutine supf_elmrhuYConv3d      
   
   
   !---------------------------------------------------------------------
   !Fifth Step "First Approach" Correction constitutive equation
   !---------------------------------------------------------------------    
   subroutine supf_elmbst_5th(e,auxtens,dvolu,auxVE,dtinv,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for Fractional step : end-of-step u
      !    (v, u_int/dt) + (p_n+1-p_n, div v)
      !
      !-----------------------------------------------------------------------
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: dvolu,dtinv,auxVE
      real(rp),    intent(inout) :: elmat(auxtens,e%mnode,auxtens,e%mnode)
      integer(ip)                :: inode,idime,jnode
      real(rp)                   :: aux,aux1

      
      
         do inode=1,e%pnode
            do jnode=1,e%pnode
               aux1= dtinv*dvolu*auxVE
            
               elmat(1,inode,1,jnode) = e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*aux1  + elmat(1,inode,1,jnode)
               elmat(1,inode,2,jnode) = 0.0_rp  + elmat(1,inode,2,jnode)
               elmat(1,inode,3,jnode) = 0.0_rp  + elmat(1,inode,3,jnode)
               
               elmat(2,inode,1,jnode) = 0.0_rp  + elmat(2,inode,1,jnode)
               elmat(2,inode,2,jnode) = e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*aux1  + elmat(2,inode,2,jnode)
               elmat(2,inode,3,jnode) = 0.0_rp  + elmat(2,inode,3,jnode)
               
               elmat(3,inode,1,jnode) = 0.0_rp  + elmat(3,inode,1,jnode)
               elmat(3,inode,2,jnode) = 0.0_rp  + elmat(3,inode,2,jnode)
               elmat(3,inode,3,jnode) = e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*aux1  + elmat(3,inode,3,jnode)
            
               
            end do
         end do

      
   end subroutine supf_elmbst_5th   
   
   
   
   
   subroutine supf_elmrhs_corr(e,auxtens,dvolu,beta,auxVE,dtinv,grvel,gpsig,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for Fractional step : end-of-step u
      !    (v, u_int/dt) + (p_n+1-p_n, div v)
      !
      !-----------------------------------------------------------------------
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime),gpsig(auxtens)
      real(rp),    intent(in)    :: dvolu,dtinv,auxVE,beta
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,idime
      real(rp)                   :: aux,aux1

      
      
      if(e%ndime==2)then
         do inode=1,e%pnode
            aux1=(1.0_rp-beta)*dvolu*e%shape(inode,e%igaus)
            aux = e%shape(inode,e%igaus)*dvolu*auxVE*dtinv
            
            elrhs(1,inode) = gpsig(1)*aux + grvel(1,1)*aux1 + elrhs(1,inode)
            elrhs(2,inode) = gpsig(2)*aux + grvel(2,2)*aux1 + elrhs(2,inode)
            elrhs(3,inode) = gpsig(3)*aux + (grvel(1,2) + grvel(2,1))*aux1 + elrhs(3,inode)
         end do
      elseif(e%ndime==3)then
         do inode=1,e%pnode
            aux1=(1.0_rp-beta)*dvolu*e%shape(inode,e%igaus)
            aux = e%shape(inode,e%igaus)*dvolu*auxVE*dtinv
            
            elrhs(1,inode) = gpsig(1)*aux + grvel(1,1)*aux1 + elrhs(1,inode)
            elrhs(2,inode) = gpsig(2)*aux + grvel(2,2)*aux1 + elrhs(2,inode)
            elrhs(3,inode) = gpsig(3)*aux + grvel(3,3)*aux1 + elrhs(3,inode)
            elrhs(4,inode) = gpsig(4)*aux + (grvel(2,3) + grvel(3,2))*aux1 + elrhs(4,inode)
            elrhs(5,inode) = gpsig(5)*aux + (grvel(1,3) + grvel(3,1))*aux1 + elrhs(5,inode)
            elrhs(6,inode) = gpsig(6)*aux + (grvel(1,2) + grvel(2,1))*aux1 + elrhs(6,inode)
            
         end do
      end if
      
   end subroutine supf_elmrhs_corr
   
   subroutine supf_elmrhs_corr_extra(e,auxtens,dvolu,beta,grvel2,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms for Fractional step : end-of-step u
      !    (v, u_int/dt) + (p_n+1-p_n, div v)
      !
      !-----------------------------------------------------------------------
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      integer(ip), intent(in)    :: auxtens
      real(rp),    intent(in)    :: grvel2(e%ndime,e%ndime)
      real(rp),    intent(in)    :: dvolu,beta
      real(rp),    intent(inout) :: elrhs(auxtens,e%mnode)
      integer(ip)                :: inode,idime
      real(rp)                   :: aux

      
      
      if(e%ndime==2)then
         do inode=1,e%pnode
            aux=-(1.0_rp-beta)*dvolu*e%shape(inode,e%igaus)
            elrhs(1,inode) = grvel2(1,1)*aux + elrhs(1,inode)
            elrhs(2,inode) = grvel2(2,2)*aux + elrhs(2,inode)
            elrhs(3,inode) = (grvel2(1,2) + grvel2(2,1))*aux + elrhs(3,inode)
         end do
      elseif(e%ndime==3)then
         do inode=1,e%pnode
            aux=-(1.0_rp-beta)*dvolu*e%shape(inode,e%igaus)
            elrhs(1,inode) = grvel2(1,1)*aux + elrhs(1,inode)
            elrhs(2,inode) = grvel2(2,2)*aux + elrhs(2,inode)
            elrhs(3,inode) = grvel2(3,3)*aux + elrhs(3,inode)
            elrhs(4,inode) = (grvel2(2,3) + grvel2(3,2))*aux + elrhs(4,inode)
            elrhs(5,inode) = (grvel2(1,3) + grvel2(3,1))*aux + elrhs(5,inode)
            elrhs(6,inode) = (grvel2(1,2) + grvel2(2,1))*aux + elrhs(6,inode)
            
         end do
      end if
      
   end subroutine supf_elmrhs_corr_extra      


end module