subroutine dombou(a)
   !-----------------------------------------------------------------------
   !****f* Domain/dombou
   ! NAME
   !    dombou
   ! DESCRIPTION
   !    Compute a%nboun, lnodb, a%lboel and a%ltypb automatically
   !-----------------------------------------------------------------------
   use typre
   use Mod_Mesh
   implicit none 
   
   class(FemMesh), target :: a
   
   integer(ip), allocatable :: faces(:,:,:),bomsh(:,:)
   integer(ip)              :: istat,ielty,iface,ilist
   integer(ip)              :: jelem,jface,ndofb,jelty,ielem,iboun,inodb,inode,ipoin,ispos,ispos2,jpoin,jspos
   logical(lg)              :: equal_faces
   integer(ip), pointer     :: tlboel(:) => NULL()
   logical :: auxcount
   
   call a%Timer%Dombou%Tic
   
   !Construct cfael
   do ielty=1,a%nelty
      call dombou_nface(a%ndime,a%nnode(ielty),a%nface(ielty))             !Set NFACE
   end do
   a%mface=maxval(a%nface)
   
   call a%Memor%alloc(a%mnodb,a%mface,a%nelem,faces,'faces','dombou')
   call a%Memor%alloc(a%mface,a%nelem,bomsh,'bomsh','dombou')
   call a%Memor%alloc(a%mnodb,a%mface,a%nelty,a%cfael,'cfael','dombou')
   do ielty=1,a%nelty
      call dombou_cfael(&                                                !Set CFAEL, array
            a%mnodb,a%nnodb(ielty),a%nface(ielty),a%ndime,&                 !conecting el & fa 
            a%nnode(ielty),a%cfael(1,1,ielty)) 
   end do
   
   !Construct faces
   ielty = 1
   do ielem=1,a%nelem   
      ispos = a%nnode(1)*(ielem-1)
      if (a%nelty > 1) then 
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
         ispos = a%pnods(ielem)-1
      endif
      do iface=1,a%nface(ielty)
         do inodb=1,a%nnodb(ielty)
            inode=a%cfael(inodb,iface,ielty)
            ipoin = a%lnods(ispos+inode)
            faces(inodb,iface,ielem)=a%lnods(ispos+inode)
         end do
         call dombou_sortb(a%nnodb(ielty),&                              !Sorte FACES
               faces(1,iface,ielem))
      end do
   end do
   
   !Compute a%nboun and fill in bomsh
   a%nboun=0
   ielty = 1
   jelty = 1
   do ielem=1,a%nelem                                                    !Compare and eliminate repited faces
      ispos = a%nnode(1)*(ielem-1)                               
      if (a%nelty > 1) then 
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
         ispos = a%pnods(ielem)-1
      endif
      do iface=1,a%nface(ielty)
         !Check if all the nodes in the face belong to another domain, if so eliminate
         auxcount = .false.
         do inodb = 1,a%nnodb(ielty)
            ipoin = faces(inodb,iface,ielem)
            if (ipoin <= a%npoinLocal) auxcount = .true.
            if (a%kfl_perio == 1) then
               if (ipoin > 0) then
                  if (a%MasterSlave(ipoin) <= a%npoinLocal) auxcount = .true.
               endif
            endif   
         enddo
         if (.not. auxcount) then
               !eliminate face
               faces(1,iface,ielem) = 0
         else
               ipoin=faces(1,iface,ielem)
               if(ipoin/=0) then
                  ilist = a%pelpo(ipoin)
                  do while (ilist < a%pelpo(ipoin+1)) 
                     jelem=a%lelpo(ilist)
                     if(jelem/=ielem) then
                        jspos = a%nnode(1)*(ielem-1)                               !eliminate the repited faces
                        if (a%nelty > 1) jelty = a%ltype(a%pnods(jelem+1)-a%pnods(jelem))
                        jface=0
                        do while(jface/=a%nface(jelty))
                           jface=jface+1
                           if(faces(1,jface,jelem)/=0) then
                              equal_faces=.true.
                              inodb=0
                              do while(equal_faces.and.&
                                    &   inodb/=a%nnodb(ielty))
                                 inodb=inodb+1
                                 if(faces(inodb,iface,ielem)&
                                       & /=faces(inodb,jface,jelem))&
                                       equal_faces=.false.
                              end do
                              if(equal_faces) then
                                 faces(1,iface,ielem)=0                 !IFACE and JFACE
                                 faces(1,jface,jelem)=0                 !are eliminated
                                 jface=a%nface(jelty)                     !Exit JFACE do
                                 ilist=a%pelpo(ipoin+1)                 !Exit JELEM do
                              end if
                           end if
                        end do
                     end if
                     ilist=ilist+1
                  end do
                  if(faces(1,iface,ielem)/=0) then
                     a%nboun=a%nboun+1
                     do inodb=1,a%nnodb(ielty)
                        inode=a%cfael(inodb,iface,ielty)
                        faces(inodb,iface,ielem)=a%lnods(ispos+inode)
                     end do
                     bomsh(iface,ielem)=a%nboun
                  end if
               end if
         endif
      end do
   end do
   
   
   if (a%nelty > 1) then 
      call a%Memor%palloc(a%nboun*(a%mnodb+1),tlboel,'tlboel','dombou')
   else
      call a%Memor%alloc(a%nboun*(a%mnodb+1),a%lboel,'lboel','dombou')
      tlboel => a%lboel
   endif
   call a%Memor%alloc(a%nboun+1,a%pboel,'pboel','dombou')
   
   !Compute lnodb and a%lboel
   ielty = 1
   a%pboel(1) = 1
   do ielem=1,a%nelem
      ispos = a%nnode(1)*(ielem-1)                               
      if (a%nelty > 1) then 
         ielty = a%ltype(a%pnods(ielem+1)-a%pnods(ielem))
         ispos = a%pnods(ielem)-1
      endif
      do iface=1,a%nface(ielty)
         if(faces(1,iface,ielem)/=0) then
            iboun=bomsh(iface,ielem)
            a%pboel(iboun+1) = a%pboel(iboun) + a%nnodb(ielty)+1
            ispos2 = a%pboel(iboun)-1
            tlboel(ispos2+a%nnodb(ielty)+1) = ielem

            do inodb=1,a%nnodb(ielty)
               ipoin= faces(inodb,iface,ielem)
               nodes: do inode=1,a%nnode(ielty)
                  jpoin=a%lnods(ispos+inode)
                  if(ipoin==jpoin) then
                     tlboel(ispos2+inodb)=inode
                     exit nodes
                  end if
               end do nodes
            end do

         end if
      end do
   end do

   if (a%nelty > 1) then
      call a%Memor%alloc(a%pboel(a%nboun+1)-1,a%lboel,'lboel','dombou')
      a%lboel = tlboel(1:(a%pboel(a%nboun+1)-1))
      call a%Memor%pdealloc(a%nboun*(a%mnodb+1),tlboel,'tlboel','dombou')
   else
      call a%Memor%dealloc(a%nboun+1,a%pboel,'pboel','dombou')
      call a%Memor%alloc(1,a%pboel,'pboel','dombou')
      a%pboel(1) = a%mnodb
   endif

   call a%Memor%dealloc(a%mnodb,a%mface,a%nelem,faces,'faces','dombou')
   call a%Memor%dealloc(a%mface,a%nelem,bomsh,'bomsh','dombou')

   call a%Timer%Dombou%Toc
   
end subroutine dombou

subroutine dombou_nface(ndime,nnode,nface)

!-----------------------------------------------------------------------
!
!    Compute nface
!
!-----------------------------------------------------------------------
  use      typre
  implicit none 

  integer(ip), intent(in)  ::  ndime,nnode
  integer(ip), intent(out) ::  nface
      
  if (ndime==2) then 
     if (nnode<=4) then 
        nface=nnode
     else
        if (nnode==6.or.nnode==7.or.nnode==10.or.nnode.eq.15) then !   p2, ?,p3,p4
           nface=3
        else
           nface=4
        end if
     end if
  else
     if (nnode==4) then
        nface=4
     else if (nnode==10) then
        nface=4
     else if (nnode==15) then
        nface=4
     else if (nnode==8) then
        nface=6
     else if (nnode==20) then
        nface=6
     else if (nnode==27) then
        nface=6
     else if (nnode==6) then
        nface=0        !!! For prisms there is not only one nface value
     end if
  end if
end subroutine dombou_nface

subroutine dombou_sortb(n,a)

!-----------------------------------------------------------------------
!
!    Sort a vector
!
!-----------------------------------------------------------------------

  use      typre
  implicit none
  integer(ip)   n,a(n),i,j,t
  
  do i=1,n-1
     do j=i+1,n
        if (a(i)>a(j)) then
           t   =a(i)
           a(i)=a(j) 
           a(j)=t
        end if
     end do
  end do

end subroutine dombou_sortb

subroutine dombou_cfael(mnodb,nnodb,nface,ndime,nnode,cfael)

!-----------------------------------------------------------------------
!
!   Construct CFAEL that contains the conectivity between 
!   faces and elements
!
!-----------------------------------------------------------------------

  use      typre
  implicit none

  integer(ip), intent(in)  :: mnodb,nnodb,nface,ndime,nnode
  integer(ip), intent(out) :: cfael(mnodb,nface)
  integer(ip)              :: iface

  iface=0
  if (ndime==2) then
     if (nnode==3) then       !P1
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=1
     else if (nnode==4) then  !Q1
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=4
        iface=iface+1
        cfael(1,iface)=4
        cfael(2,iface)=1
     else if (nnode==6.or.nnode==7) then  !P2
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=4
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=5
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=1
        cfael(3,iface)=6
     else if (nnode==8.or.nnode==9) then  !Q2
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=5
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=6
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=4
        cfael(3,iface)=7
        iface=iface+1
        cfael(1,iface)=4
        cfael(2,iface)=1
        cfael(3,iface)=8
     else if (nnode.eq.16) then  ! Q3
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=5
        cfael(4,iface)=6
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=7
        cfael(4,iface)=8        
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=4
        cfael(3,iface)=9
        cfael(4,iface)=10
        iface=iface+1
        cfael(1,iface)=4
        cfael(2,iface)=1
        cfael(3,iface)=11
        cfael(4,iface)=12    
      else if (nnode.eq.10) then ! P3
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=4
        cfael(4,iface)=5
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=6
        cfael(4,iface)=7        
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=1
        cfael(3,iface)=8
        cfael(4,iface)=9

!  Añadido 
      else if (nnode==15) then ! P4
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=4
        cfael(4,iface)=5
        cfael(5,iface)=6
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=7
        cfael(4,iface)=8
        cfael(5,iface)=9        
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=1
        cfael(3,iface)=10
        cfael(4,iface)=11
        cfael(5,iface)=12
     else if (nnode==25) then  ! Q4
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=5
        cfael(4,iface)=6
        cfael(5,iface)=7
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=8
        cfael(4,iface)=9
        cfael(5,iface)=10        
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=4
        cfael(3,iface)=11
        cfael(4,iface)=12
        cfael(5,iface)=13
        iface=iface+1
        cfael(1,iface)=4
        cfael(2,iface)=1
        cfael(3,iface)=14
        cfael(4,iface)=15
        cfael(5,iface)=16  
! Hasta Aqui 
     end if
  else if (ndime==3) then
     if (nnode==4.or.nnode==5) then
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=3
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=4
        cfael(3,iface)=3
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=4
        cfael(3,iface)=2
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=3
        cfael(3,iface)=4
     else if (nnode==8) then
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=3
        cfael(4,iface)=4
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=7
        cfael(4,iface)=6
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=4
        cfael(3,iface)=8
        cfael(4,iface)=7
        iface=iface+1
        cfael(1,iface)=4
        cfael(2,iface)=1
        cfael(3,iface)=5
        cfael(4,iface)=8
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=6
        cfael(4,iface)=5
        iface=iface+1
        cfael(1,iface)=5
        cfael(2,iface)=6
        cfael(3,iface)=7
        cfael(4,iface)=8
     else if (nnode==10.or.nnode==15) then
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=4
        cfael(4,iface)=5
        cfael(5,iface)=9
        cfael(6,iface)=8
        iface=iface+1
        cfael(1,iface)=2
        cfael(2,iface)=3
        cfael(3,iface)=4
        cfael(4,iface)=6
        cfael(5,iface)=10
        cfael(6,iface)=9
        iface=iface+1
        cfael(1,iface)=3
        cfael(2,iface)=1
        cfael(3,iface)=4
        cfael(4,iface)=7
        cfael(5,iface)=8
        cfael(6,iface)=10
        iface=iface+1
        cfael(1,iface)=1
        cfael(2,iface)=2
        cfael(3,iface)=3
        cfael(4,iface)=5
        cfael(5,iface)=6
        cfael(6,iface)=7
        if (nnode==15) then                             ! OJO revisar
           cfael(7,1)=12
           cfael(7,2)=13
           cfael(7,3)=14
           cfael(7,4)=11
        end if
     else if (nnode==20.or.nnode==27) then           ! OJO completar
     endif
  end if

end subroutine dombou_cfael

