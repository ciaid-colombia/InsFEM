subroutine getlnods(ielem,pnode,ispos,pnods,lnods,nelty,nnode)
   use typre
   implicit none
   
   integer(ip) :: ielem,pnode,pnods(*),nelty,nnode(*)
   integer(ip) :: ispos
   integer(ip), target  :: lnods(*)
   
   if (nelty == 1) then
      pnode = nnode(1)
      ispos = pnode*(ielem-1)
   else
      pnode = pnods(ielem+1)-pnods(ielem)
      ispos = pnods(ielem)-1
   endif   

end subroutine
