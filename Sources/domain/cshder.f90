subroutine cshder(a)
!-----------------------------------------------------------------------
!****f* Domain/cshder
! NAME
!    cshder
! DESCRIPTION
!    Calculates a%shape functions and a%derivatives in the reference domain
!    and call a%Memor%allocs arrays calculations on the physical domain (performed
!    at each loop over elements). Calculates the lumped mass matrix.
! OUTPUT
!    SHAPE(a%mnode,a%mgaus,a%nelty)
!    DERIV(a%ndime,a%mnode,a%mgaus,a%nelty)
!    HESLO(a%ntens,a%mnode,a%mgaus,a%nelty)
!    if(a%nboun>0) :
!       WEIGP(a%mgaus,a%nelty)
!       SHAPB(a%mnodb,a%mgaub,a%nelty)
!       DERIB(a%ndimb,a%mnodb,a%mgaub,a%nelty)
!       WEIGB(a%mgaub,a%nelty)
!       SHAGA(a%mgaus,a%mnode,a%nelty)
!       SHAGB(a%mgaub,a%mnodb,a%nelty)
!-----------------------------------------------------------------------
  use typre
  use def_parame
  use Mod_Mesh
  implicit none
  
  class(FemMesh), intent(inout) :: a


  integer(ip) :: ielty,ielem,igaub,igaus,inodb,inode,ipoin,ispos,linea,pnode
  real(rp)    :: posgp(a%ndime,a%mgaus,a%nelty)         !Domain integration rules (defined in setdim)
  real(rp)    :: posgb(a%ndimb,a%mgaub,a%nelty)         !Boundary integration rules
  real(rp)    :: heslb(a%ntens,a%mnodb,a%mgaub,a%nelty) !Used as dummy at the boundary

  real(rp)    :: poscg(a%ndime)                         !Rules to calculate functions at the cog
  real(rp)    :: hescg(a%ntens,a%mnode,a%nelty)         !Used as dummy at cog

  real(rp)    :: posgc(a%ndime,a%mnode,a%nelty)         !Closed integration rules
  real(rp)    :: shapc(a%mnode,a%mnode,a%nelty)
  real(rp)    :: deric(a%ndime,a%mnode,a%mnode,a%nelty)
  real(rp)    :: heslc(a%ntens,a%mnode,a%mnode,a%nelty)
  real(rp)    :: weigc(a%mnode,a%nelty)

  real(rp)    :: posgw(a%ndime,a%mnode,a%nelty)         !Rules to calculate a%shaga y a%shagb
  real(rp)    :: weigw(a%mnode,a%nelty)
  
  real(rp)    :: elcod(a%ndime,a%mnode),cartd(a%ndime,a%mnode),xjaci(a%ndime,a%ndime),xjacm(a%ndime,a%ndime)
  real(rp)    :: detjm,dvolu
  

  integer(ip) :: istat,prule
  
  integer(ip) :: llapl,ltopo,lrule,nnodb,ngaub
  real(rp)    :: hnatu

  !Element a%shape function and a%derivatives (a%shape,a%deriv,a%heslo,a%weigp) 
  call a%Memor%alloc(a%mnode,a%mgaus,a%nelty,a%shape,'shape','cshder')
  call a%Memor%alloc(a%ndime,a%mnode,a%mgaus,a%nelty,a%deriv,'deriv','cshder')
  call a%Memor%alloc(a%ntens,a%mnode,a%mgaus,a%nelty,a%heslo,'heslo','cshder')
  call a%Memor%alloc(a%mgaus,a%nelty,a%weigp,'weigp','cshder') 

  !Evaluate a%shapes using the integration rule defined in setdim
  do ielty=1,a%nelty
     call rulepw(a%ndime,a%ngaus(ielty),a%lrule(ielty),&
          a%posgp(1,1,ielty),a%weigp(1,ielty))
     do igaus=1,a%ngaus(ielty)
        call shafun(&
             a%posgp(1,igaus,ielty),a%ndime,a%nnode(ielty),a%ntens,&
             a%shape(1,igaus,ielty),a%deriv(1,1,igaus,ielty),&
             a%heslo(1,1,igaus,ielty))
     end do
  end do

  !Element Center of gravity a%shape function and a%derivatives (a%shacg,a%dercg,a%weicg) 
  call a%Memor%alloc(a%mnode,a%nelty,a%shacg,'shacg','cshder')
  call a%Memor%alloc(a%ndime,a%mnode,a%nelty,a%dercg,'dercg','cshder')
  call a%Memor%alloc(a%nelty,a%weicg,'weicg','cshder')  
  
  do ielty=1,a%nelty
     if(a%ltopo(ielty)==0) then      !Quad/Hexa
        prule=1
     else if(a%ltopo(ielty)==1) then !Tria/Tetr
        prule=3
     else if(a%ltopo(ielty)==2) then !Prisms
        prule=5
     end if
     call rulepw(a%ndime,one,prule,poscg,a%weicg(ielty))
     call shafun(&
&         poscg,a%ndime,a%nnode(ielty),a%ntens,&
&         a%shacg(1,ielty),a%dercg(1,1,ielty),hescg(1,1,ielty))
  end do
  
  !Element shape functions and derivatives for closed integration rule (vmass, projections)
  call a%Memor%alloc(a%mnode,a%mnode,a%nelty,a%shape_clo,'shape_clo','cshder')
  call a%Memor%alloc(a%ndime,a%mnode,a%mnode,a%nelty,a%deriv_clo,'deriv_clo','cshder')
  call a%Memor%alloc(a%ntens,a%mnode,a%mnode,a%nelty,a%heslo_clo,'heslo_clo','cshder')
  call a%Memor%alloc(a%mnode,a%nelty,a%weigp_clo,'weigp_clo','cshder') 
  
  do ielty=1,a%nelty
     !This is a fake setdim, only for closed quadrature. What I need here is lrule 
     call setdim(&
     a%nnode(ielty),a%nnode(ielty),1_ip,a%ndime,&
     llapl,ltopo,lrule,nnodb,ngaub,hnatu) 
  
     call rulepw(a%ndime,a%nnode(ielty),lrule,&
         posgw(1,1,ielty),a%weigp_clo(1,ielty))
     do igaus=1,a%nnode(ielty)
        call shafun(&
             posgw(1,igaus,ielty),a%ndime,a%nnode(ielty),a%ntens,&
             a%shape_clo(1,igaus,ielty),a%deriv_clo(1,1,igaus,ielty),&
             a%heslo_clo(1,1,igaus,ielty))
     end do
  end do
  
  
  
  
  
  
  

  !Boundary a%shape function and a%derivatives (a%shapb,a%derib,heslb,a%weigb)
  call a%Memor%alloc(a%mnodb,a%mgaub,a%nelty,a%shapb,'shapb','cshder')
  call a%Memor%alloc(a%ndimb,a%mnodb,a%mgaub,a%nelty,a%derib,'derib','cshder')        
  call a%Memor%alloc(a%mgaub,a%nelty,a%weigb,'weigb','cshder')
  
  do ielty=1,a%nelty
     call rulepw(a%ndimb,a%ngaub(ielty),a%lrule(ielty),&
          posgb(1,1,ielty),a%weigb(1,ielty))
     do igaub=1,a%ngaub(ielty)
        call shafun(&
             posgb(1,igaub,ielty),a%ndimb,a%nnodb(ielty),a%ntens,&
             a%shapb(1,igaub,ielty),&
             a%derib(1,1,igaub,ielty),&
             heslb(1,1,igaub,ielty))
     end do
  end do

  !Computes the interpolation functions associated to the integration
  !points in order to extrapolate values from these integration points
  !to the nodes
  call a%Memor%alloc(a%mgaus,a%mnode,a%nelty,a%shaga,'shaga','cshder')
  do ielty=1,a%nelty
     call rulepw(a%ndime,a%nnode(ielty),(a%ltopo(ielty)+1)*2,&
          posgw(1,1,ielty),weigw (1,ielty))
     do inode=1,a%nnode(ielty)   
        call shafga(&
             posgw(1,inode,ielty),a%ndime,a%ltopo(ielty),a%ngaus(ielty),&
             a%shaga(1,inode,ielty))
     end do
  end do

  call a%Memor%alloc(a%mgaub,a%mnodb,a%nelty,a%shagb,'shagb','cshder')
  do ielty=1,a%nelty
     call rulepw(a%ndimb,a%nnodb(ielty),(a%ltopo(ielty)+1)*2,&
          posgw(1,1,ielty),weigw(1,ielty))
     do inodb=1,a%nnodb(ielty)   
        call shafga(&
             posgw(1,inodb,ielty),a%ndimb,a%ltopo(ielty),a%ngaub(ielty),&
             a%shagb(1,inodb,ielty))
     end do
  end do
  
  !Ltypb
  call a%Memor%alloc(a%mnodb,a%ltypb,'ltypb','dombou')
   a%ltypb = 0
   do ielty = 1,a%nelty
      a%ltypb(a%nnodb(ielty)) = ielty
   enddo

end subroutine cshder
