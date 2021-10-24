subroutine Finbou(a,pnodb,knodb,iboun)
!-----------------------------------------------------------------------
!    This routine looks for iboun in in the list knodb using
!    a%pbopo and a%lbopo
!-----------------------------------------------------------------------
  use typre
  use Mod_Mesh
  implicit none
  
  class(FemMesh) :: a
  integer(ip) :: pnodb,knodb(pnodb)
  integer(ip) :: iboun
  
  integer(ip)              :: onodb,inook,ibook,ibopo,ielem,inodb,inode,ispos,ispos2,jnodb
  integer(ip)              :: qnodb,jpoin,pblty

  
  
  iboun=0
  ibook=0
  pblty=1
  !loop_iboun: do while(iboun<nboun)
  !   iboun=iboun+1
  loop_iboun: do ibopo = a%pbopo(knodb(1)),a%pbopo(knodb(1)+1)-1
	iboun = a%lbopo(ibopo)
   
     ispos2 = (a%nnodb(1)+1)*(iboun-1)                               
     if (a%nelty > 1) then 
       ispos2 = a%pboel(iboun)-1
       pblty = a%ltypb(a%pboel(iboun+1)-a%pboel(iboun)-1)
     endif
     ielem = a%lboel(ispos2+a%nnodb(pblty)+1) 
     ispos = a%nnode(1)*(ielem-1)                               
     if (a%nelty > 1) ispos = a%pnods(ielem)-1
     
     qnodb=a%nnodb(pblty)
     if(pnodb==qnodb) then
        onodb=0
        jnodb=0
        inook=0
        loop_jnodb: do while(jnodb<pnodb)
           jnodb=jnodb+1
           jpoin=knodb(jnodb)
           inodb=0
           loop_inodb: do while(inodb<qnodb)
              inodb=inodb+1
              inode = a%lboel(ispos2+inodb)
              if(jpoin==a%lnods(ispos+inode)) then
                 onodb=onodb+1
                 inook=1
                 exit loop_inodb
              end if
           end do loop_inodb
           if(inook/=1) exit loop_jnodb
        end do loop_jnodb
        if(onodb==pnodb) then
           ibook=1
           exit loop_iboun
        end if
     end if
  end do loop_iboun
  if(ibook/=1) iboun=0

end subroutine Finbou
