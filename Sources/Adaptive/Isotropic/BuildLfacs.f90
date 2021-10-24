subroutine BuildLfacs(a)
   use typre
   use Mod_IsotropicAdaptiveRefiner
   use Mod_AdaptiveElementDependentSubroutines
   implicit none
   class(IsotropicAdaptiveRefiner) :: a
   
   integer(ip), allocatable :: pelpo(:),peaux(:)
   integer(ip), allocatable :: lelpo(:)

   integer(ip) :: ielem,inode,iface,cface,iv,jelem,jface,isequal,ipoin
   integer(ip) :: ispos,jspos,ifpos,jfpos
   integer(ip) :: facnod(a%mnodb,a%mface), facnod2(a%mnodb,a%mface)
   
   integer(ip) :: nchild, nface, nchild2,nface2,pnodb,pnodb2,pnode,pnode2,nchildface,nchildface2
   integer(ip) :: ielty,ielty2

   !FACES
   !build list of neighbour elements
   !It is stored in pfacs and lfacs (CSR). 
   !pfacs: starting point for each elements
   !lfacs: (1): corresponding element
   !lfacs: (2): corresponding face in the element

   
   !Count dimension of lfacs and allocate
   allocate(a%pfacs(a%nelem+1))
   a%pfacs(1) = 1
   do ielem = 1,a%nelem
      pnode = a%pnods(ielem+1)-a%pnods(ielem)
      ielty = a%ElementTypeList(ielem)%ielty
      call ElementGetDimensions(ielty,pnodb,nchild,nface,nchildface)
      a%pfacs(ielem+1) = a%pfacs(ielem) + nface
   enddo
   allocate(a%lfacs(2,a%pfacs(a%nelem+1)-1))
   a%lfacs(:,:) = 0
   
   !Create Pelpo and Lelpo
   allocate(pelpo(a%npoin+1))
   pelpo(:) = 0
   pelpo(1) = 1
   do ielem = 1,a%nelem
      ispos = a%pnods(ielem)-1
      pnode = a%pnods(ielem+1)-a%pnods(ielem)
      do inode = 1,pnode
         ipoin = a%lnods(ispos+inode)
         pelpo(ipoin+1) = pelpo(ipoin+1)+1
      enddo
   enddo
   
   ! Compress pelpo
   do ipoin=1,a%npoin
      pelpo(ipoin+1)=pelpo(ipoin+1)+pelpo(ipoin)
   end do
   !Create lelpo
   allocate(peaux(a%npoin))
   peaux(:)=-1
   allocate(lelpo(pelpo(a%npoin+1)))
   
   do ielem = 1,a%nelem
      ispos = a%pnods(ielem)-1
      pnode = a%pnods(ielem+1)-a%pnods(ielem)
      do inode=1,pnode
         ipoin=a%lnods(ispos+inode)
         peaux(ipoin)=peaux(ipoin)+1
         lelpo(pelpo(ipoin)+peaux(ipoin))=ielem
      end do
   end do
   
   
   ielemdo: do ielem = 1,a%nelem
      !build faces of ielem
      ispos = a%pnods(ielem)-1
      pnode = a%pnods(ielem+1)-a%pnods(ielem)
      ielty = a%ElementTypeList(ielem)%ielty
      call face_from_el(ielty,nface,pnodb,a%lnods(ispos+1),0,1,facnod)
      cface = 0
      
      ifpos = a%pfacs(ielem)-1
      do iface = 1,nface
         if (a%lfacs(1,ifpos+iface) /= 0) cface = cface+1  !faces already done
      enddo        
      
      if (cface < nface) then
         inodedo: do inode = 1,pnode
            ipoin = a%lnods(ispos+inode)
            do iv = pelpo(ipoin),pelpo(ipoin+1)-1
                  jelem = lelpo(iv)
                  if (jelem > ielem) then
                     
                     jspos = a%pnods(jelem)-1
                     pnode2 = a%pnods(jelem+1) - a%pnods(jelem)
                     ielty2 = a%ElementTypeList(jelem)%ielty
                     call ElementGetDimensions(ielty2,pnodb2,nchild2,nface2,nchildface2)
                     !compare faces and fill neighbour elements list
                     call face_from_el(ielty2,nface2,pnodb2,a%lnods(jspos+1),0,1,facnod2)
                     ifacedo: do iface = 1,nface
                        if (a%lfacs(1,ifpos+iface) == 0 .and. pnodb == pnodb2) then
                              jfacedo: do jface = 1,nface2
                                 jfpos = a%pfacs(jelem)-1
                                 if(a%lfacs(1,jfpos+jface) == 0) then
                                    call face_compare(facnod(1,iface),facnod2(1,jface),pnodb,isequal)
                                    if (isequal == 1) then
                                          a%lfacs(1,ifpos+iface) = jelem
                                          a%lfacs(2,ifpos+iface) = jface
                                          a%lfacs(1,jfpos+jface) = ielem
                                          a%lfacs(2,jfpos+jface) = iface
                                          cface = cface+1
                                          exit jfacedo
                                    endif
                                 endif
                              enddo jfacedo
                        endif
                        if (cface == nface) exit inodedo
                     enddo ifacedo
                     
                  endif
            enddo
         enddo inodedo  
      endif
   enddo ielemdo    
   
   deallocate(pelpo,peaux)
   deallocate(lelpo)
 end subroutine

   