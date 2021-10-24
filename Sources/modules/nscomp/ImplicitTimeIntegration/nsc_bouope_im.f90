subroutine nsc_bouope_im(NSCompImplicitProblem)
   use typre
   use Mod_nsc_BaseElmope
   use Mod_nsc_elmdir_im
   use Mod_nsc_HangingNodes
   use Mod_nsc_bouback_im
   implicit none
   class(NSCompressibleImplicitProblem), target :: NSCompImplicitProblem
   
   real(rp), allocatable :: wmatr(:,:,:,:)
   real(rp), allocatable :: wrhsi(:,:)

   integer(ip) :: iboun,nboun,ndime
   integer(ip) :: igaub,inodb,inode,idofn
   integer(ip) :: kfl_HangingNodes

   real(rp), allocatable :: boden(:),elbde(:)
   real(rp), allocatable :: bomom(:,:),elbmo(:,:),gpbmo(:)
   real(rp), allocatable :: boene(:),elben(:)
   real(rp), allocatable :: tract(:)

   real(rp)    :: dsurf

   a=>NSCompImplicitProblem
   

   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
  
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_bouope_im')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','nsc_bouope_im')
   
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,wmatr,'wmatr','nsc_bouope_im')
   call a%Memor%alloc(a%ndofn,e%mnode,wrhsi,'wrhsi','nsc_bouope_im')

   call a%Mesh%GetHanging(kfl_HangingNodes)

   !Variables!
   call a%Memor%alloc(a%ndofn,tract,'tract','nsc_bouope_im')
   call a%Memor%alloc(e%pnodb,boden,'boden','nsc_bouope_im')
   call a%Memor%alloc(e%ndime,e%pnodb,bomom,'bomom','nsc_bouope_im')
   call a%Memor%alloc(e%pnodb,boene,'boene','nsc_bouope_im')
   call a%Memor%alloc(e%ndime,gpbmo,'gpbmo','nsc_bouope_im')
   call a%Memor%alloc(e%mnode,elbde,'elbde','nsc_bouope_im')
   call a%Memor%alloc(e%ndime,e%mnode,elbmo,'elbmo','nsc_bouope_im')
   call a%Memor%alloc(e%mnode,elben,'elben','nsc_bouope_im')

   ! Loop over boundaries
   boundaries: do iboun=1,nboun

      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      elmat=0.0_rp
      elrhs=0.0_rp

      if(    a%kfl_fixbo(iboun)==3.or.&   !Mom=Neu  Ene=Free
         & a%kfl_fixbo(iboun)==4.or.&     !Mom=Neu  Ene=Neu
         & a%kfl_fixbo(iboun)==5.or.&     !Mom=Neu  Ene=Rob
         & a%kfl_fixbo(iboun)==6.or.&     !Mom=Free  Ene=Neu
         & a%kfl_fixbo(iboun)==7.or.&     !Mom=Free  Ene=Rob
         & a%kfl_fixbo(iboun)==8.or.&     !Mom=Wall_law Ene=Free
         & a%kfl_fixbo(iboun)==9 &         !Mom= backlow penalty
         ) then  

         !Physical Parameters
         call a%GetPhysicalParameters(acvis,actco,accph,accvh)
         
         call e%elmdel
         
         !Gather operations
         call e%gatherb(1_ip,boden,a%densf(:,1))
         call e%gatherb(e%ndime,bomom,a%momen(:,:,1))
         call e%gatherb(1_ip,boene,a%energ(:,1))
         call e%gather(1_ip,elbde,a%densf(:,1))       
         call e%gather(e%ndime,elbmo,a%momen(:,:,1))      
         call e%gather(1_ip,elben,a%energ(:,1))       
         
         dsurf = 0.0_rp
         
         !Gauss-Point Loop
         gauss_points : do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Initialize
            wmatr=0.0_rp
            wrhsi=0.0_rp
            tract=0.0_rp
      
            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            
            !Derivatives at the boundary
            call e%elmderb
   
            !Backflow penalty
            if(a%kfl_fixbo(iboun)==9) then
               call e%interpb(e%ndime,bomom,gpbmo)
               call nsc_bouback_im(e,gpbmo,a%bvnat(iboun)%a(1),wmatr)

            end if

            elmat = elmat + dsurf*wmatr            
           
            do inodb=1,e%pnodb
               inode = e%lboel(inodb)
               do idofn=1,a%ndofn
                     elrhs(idofn,inode)=elrhs(idofn,inode)+dsurf*(wrhsi(idofn,inode)+tract(idofn)*e%shapb(inodb,e%igaub))
               end do   
              
           end do
   
         end do gauss_points

      end if 

      if (kfl_HangingNodes == 1) call ModifyMatricesHanging
      
      !Boundary conditions
      call nsc_elmdir(a,e,elmat,elrhs)

      
      ! Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)

   end do boundaries

   call a%Memor%dealloc(a%ndofn,tract,'tract','nsc_bouope_im')
   call a%Memor%dealloc(e%pnodb,boden,'boden','nsc_bouope_im')
   call a%Memor%dealloc(e%ndime,e%pnodb,bomom,'bomom','nsc_bouope_im')
   call a%Memor%dealloc(e%pnodb,boene,'boene','nsc_bouope_im')
   call a%Memor%dealloc(e%ndime,gpbmo,'gpbmo','nsc_bouope_im')
   call a%Memor%dealloc(e%mnode,elbde,'elbde','nsc_bouope_im')
   call a%Memor%dealloc(e%ndime,e%mnode,elbmo,'elbmo','nsc_bouope_im')
   call a%Memor%dealloc(e%mnode,elben,'elben','nsc_bouope_im')
   
   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','nsc_bouope_im')
   call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','nsc_bouope_im')

   call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,wmatr,'wmatr','nsc_bouope_im')
   call a%Memor%dealloc(a%ndofn,e%mnode,wrhsi,'wrhsi','nsc_bouope_im')

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   
end subroutine nsc_bouope_im
