subroutine blbopo(a) 
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a

   integer(ip) :: iaux,iaux2,iboun,ielem,inodb,inode,ipoin,ispos,ispos2,pblty

   call a%Timer%Blbopo%Tic
   
   pblty=1
   call a%Memor%alloc(a%npoin+1,a%pbopo,'pbopo','blbopo')
   !first loop to count
   do iboun = 1,a%nboun
      ispos2 = (a%nnodb(1)+1)*(iboun-1)                      
      if (a%nelty > 1) then 
         ispos2 = a%pboel(iboun)-1
         pblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
      endif
      ielem = a%lboel(ispos2+a%nnodb(pblty)+1) 
      ispos = a%nnode(1)*(ielem-1)                               
      if (a%nelty > 1) ispos = a%pnods(ielem)-1
      do inodb = 1,a%nnodb(pblty)
         inode = a%lboel(ispos2+inodb)
         ipoin = a%lnods(ispos+inode)
         a%pbopo(ipoin) = a%pbopo(ipoin)+1
      enddo
   enddo

   !second loop to fill
   iaux = a%pbopo(1)
   a%pbopo(1) = 1
   do ipoin = 2,a%npoin
      iaux2 = a%pbopo(ipoin)
      a%pbopo(ipoin) = a%pbopo(ipoin-1) + iaux
      iaux = iaux2
   enddo
   if (a%npoin > 0) then
      a%pbopo(a%npoin+1) = a%pbopo(a%npoin) +iaux
   endif
   
   call a%Memor%alloc(a%pbopo(a%npoin+1),a%lbopo,'lbopo','blbopo')
   
   do iboun = 1,a%nboun
      ispos2 = (a%nnodb(1)+1)*(iboun-1)                               
      if (a%nelty > 1) then 
         ispos2 = a%pboel(iboun)-1
         pblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
      endif
      ielem = a%lboel(ispos2+a%nnodb(pblty)+1) 
      ispos = a%nnode(1)*(ielem-1)                               
      if (a%nelty > 1) ispos = a%pnods(ielem)-1
      do inodb = 1,a%nnodb(pblty)
         inode = a%lboel(ispos2+inodb)
         ipoin = a%lnods(ispos+inode)
         a%lbopo(a%pbopo(ipoin)) = iboun
         a%pbopo(ipoin) = a%pbopo(ipoin)+1
      enddo
   enddo
   do ipoin = a%npoin,2,-1
      a%pbopo(ipoin) = a%pbopo(ipoin-1)
   enddo
   a%pbopo(1) = 1
   
   call a%Timer%Blbopo%Toc

end subroutine
