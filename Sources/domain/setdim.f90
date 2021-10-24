subroutine setdim(&
     nnode,ngaus,lquad,ndime,&
     llapl,ltopo,lrule,nnodb,ngaub,hnatu)

!-----------------------------------------------------------------------
!      
! This routines checks possible errors and sets some dimensions
!
! ltopo ... 0 bricks
!       ... 1 simplices
!       ... 2 Prisms
!
! llapl ... 0 Laplacian is not needed
!       ... 1 Laplacian is nedded
!
! lquad ... 0 Open quadrature
!       ... 1 Closed quadrature
!
!-----------------------------------------------------------------------
  use      def_parame
  use      typre
  implicit none
  integer(ip), intent(in)  :: nnode,ngaus,ndime,lquad
  integer(ip), intent(out) :: nnodb,ngaub,llapl,ltopo,lrule
  integer(ip)              :: degel
  real(rp)   , intent(out) :: hnatu
!
! Check errors.
!
  if(ndime==2) then
     if(  (nnode/= 3).and.(nnode/= 4).and.&
          (nnode/= 6).and.(nnode/= 9).and.&
          (nnode/=10).and.(nnode/=16).and.&
          (nnode/=15).and.(nnode/=25))    &
          call  runend('COUNTV: WRONG NUMBER OF NODES PER ELEMENTS')
  else if(ndime==3) then
     if(  (nnode/= 4).and.(nnode/= 8).and.&
          (nnode/=10).and.(nnode/=27).and.&
          (nnode/=20).and.(nnode/=64).and.&
          (nnode/=35).and.(nnode/=125).and.&
          (nnode/= 6))                    &
          call  runend('COUNTV: WRONG NUMBER OF NODES PER ELEMENTS')
  else
     call runend('COUNTV: WRONG NDIME')
  end if
!
! Define some logical flags and compute derived variables.
!
  ltopo = 1                                             ! Elements are simplices
  degel = 1                                             ! Elements are linear
  if(ndime==2) then
     if((nnode== 4).or.(nnode== 9).or.(nnode==16).or.(nnode==25)) ltopo = 0
     if((nnode== 6).or.(nnode== 9)) degel = 2
     if((nnode==10).or.(nnode==16)) degel = 3
     if((nnode==15).or.(nnode==25)) degel = 4
  else if(ndime==3) then
     if((nnode== 8).or.(nnode==27).or.(nnode==64).or.(nnode==125)) ltopo = 0
     if((nnode== 6))                                              ltopo = 2
     if((nnode==10).or.(nnode== 27)) degel = 2
     if((nnode==20).or.(nnode== 64)) degel = 3
     if((nnode==35).or.(nnode==125)) degel = 4
     if((nnode==6))                  degel = 0
  end if
  if(degel==1) then
     llapl = 0
  else
     llapl = 1
  end if
  if(ndime==2.and.nnode==4) llapl=1
  if(ndime==3.and.nnode==8) llapl=1  
!
! Number of integration rule
!
  if(ltopo==0) then
     if(lquad==0) lrule=1
     if(lquad==1) lrule=2
  else if(ltopo==1) then
     if(lquad==0) lrule=3
     if(lquad==1) lrule=4
  else if(ltopo==2) then
     if(lquad==0) lrule=5
     if(lquad==1) lrule=6
  end if
!
! Number of nodes and integration points per boundary element. For
! open integration rules, the number of boundary integration points
! is determined in order to keep the accuracy on the integration rule,
! whereas in the case of closed integration it is determined by matching
! the integration points as belonging to the element and to the boundary
!
  if(ndime==2) then
     nnodb=2                       
     if(ltopo==1) then                            !simplex
        if(lquad==0) then                               !open rule
           if(ngaus<=1) then       
              ngaub=1                 !exact for P1 
           else if(ngaus<=4) then  
              ngaub=2                 !exact for P3
           else if(ngaus<=7) then  
              ngaub=3                 !exact for P5
           else if(ngaus<=15) then 
              ngaub=4                 !exact for P7
           end if
        else if(lquad==1) then                          !closed rule
           if(ngaus<=4) then
              ngaub=2
           else if(ngaus<=7) then
              ngaub=3
           else if(ngaus<=15) then
              ngaub=4
           end if
        end if
        if(nnode== 6) nnodb=3
        if(nnode==10) nnodb=4
        if(nnode==15) nnodb=5		
     else                                     !quadrilateral
        ngaub=int(sqrt(real(ngaus,rp)))
        if(nnode== 9) nnodb=3
        if(nnode==16) nnodb=4
        if(nnode==25) nnodb=5		        
     end if
  else if(ndime==3) then
     if(ltopo==1) then                            !simplex
        if(lquad==0) then                               !open rule
           if(ngaus==1) then
              ngaub=1                 !exact for P1
           else if(ngaus<=5) then   
              ngaub=ngaus-1           !exact for P2/P3
           else if(ngaus==11) then
              ngaub=6                 !exact for P4
           else if(ngaus==14) then
              ngaub=7                 !exact for P5
           end if
        else if(lquad==1) then                          !closed rule
           if(ngaus<=5) then
              ngaub=3
           else if(ngaus<=11) then
              ngaub=6
           else if(ngaus==15) then
              ngaub=7
           else if(ngaus==20) then
              ngaub=10
           end if
        end if
        nnodb=3
        if(nnode==10) nnodb=6
        if(nnode==11) nnodb=6          
        if(nnode==15) nnodb=7
     else if(ltopo==0) then                            !brick
        nnodb=4
        ngaub=nint(float(ngaus)**(2.0d0/3.0d0))
        if(nnode==20) nnodb=8
        if(nnode==27) nnodb=9
     else                                              !prism
        nnodb=0    !WARNING?
        ngaub=0
     end if
  end if
!
! Natural element length
!
  if(ltopo==1) then
     hnatu=1.0_rp			!Simplices
  else if(ltopo==0) then
     hnatu=2.0_rp			!Bricks
  else  
     hnatu=1.0_rp			!Prism
  end if

end subroutine setdim
