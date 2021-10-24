subroutine countv(a)
!-----------------------------------------------------------------------
!****f* Domain/countv
! NAME
!    countv
! DESCRIPTION
!   Calculate dimensions:
!
!   NELTY         Number of element types
!   NELEM         Number of elements
!   NBOUN         Number of boundary elements
!   NNODE(NELTY)  Number of nodes per element
!   NPOIN         Number of nodal points
!   NDIME         Number of space dimensions
!   NSLAV         Number of slave points      
!   NGIVE         Number of give and take
!-----------------------------------------------------------------------
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh)            :: a
   integer(ip)              :: ielty,istat,jelem,ipara
   integer(ip)              :: icalc=0
   integer(ip), allocatable :: knode(:)
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL()
   integer(ip), pointer  :: nnwor => NULL()
   
   !Initializations.
   a%kfl_mnelt =  0     ! =1 a%nelty is given by user, =0 computed

   a%gnpoin     = -1    ! Obligatory
   a%gnelem     = -1    ! Obligatory
   a%ndime     = -1     ! Obligatory
   
   a%nboun     = -1     ! Optional
   !nslav     = -1      ! Optional
   a%gnskew     = -1    ! Optional
   a%nelty     = -1
   !nnods     = 0
   !nnodb     = 0
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   
   
   ! Read dimensions.
   call a%Listener%listen('countv')
   geomeloop : do while(words(1)/='GEOME')
      if(words(1)=='DIMEN') then
         dimeloop : do while(words(1)/='ENDDI')
            if(words(1)=='NODAL') then
               a%gnpoin=a%Listener%getint('NODAL',0,'*NUMBER OF NODAL POINTS')
            else if(words(1)=='ELEME') then
               a%gnelem=a%Listener%getint('ELEME',0,'*NUMBER OF ELEMENTS')
            else if(words(1)=='SPACE') then
               a%ndime=a%Listener%getint('SPACE',2,'*SPACE DIMENSION')
            else if(words(1)=='SKEWS') then
               a%gnskew=a%Listener%getint('SKEWS',0,'*NUMBER OF SKEW SYSTEMS')
            else if(words(1)=='NODES') then
               a%nelty=0
               do ipara=1,nnpar
                  if(int(param(ipara))>0) a%nelty=a%nelty+1
               end do
               
               call a%Memor%alloc(a%nelty,a%nnode,'NNODE','countv')
               call a%Memor%alloc(a%nelty,a%nface,'NFACE','countv')
               a%kfl_mnelt=1
               ielty=0
               do ipara=1,nnpar
                  if(int(param(ipara))>0) then
                     ielty=ielty+1
                     a%nnode(ielty)=int(param(ipara))
                  end if
               end do
               
            end if
            call a%Listener%listen('countv')
         end do dimeloop
         if(  a%gnpoin/=-1.and.a%gnelem/=-1.and.a%ndime/=-1.and.&
               a%nboun/=-1.and.a%gnskew/=-1.and.&
               a%nelty/=-1) exit geomeloop
      end if
      
      if(a%gnpoin==-1) call runend('COUNTV: WRONG NUMBER OF NODES')
      if(a%gnelem==-1) call runend('COUNTV: WRONG NUMBER OF ELEMENTS')
      if(a%nelty==-1) call runend('COUNTV: WRONG NUMBER OF ELEMENT TYPES')
      if(a%ndime==-1) call runend('COUNTV: WRONG NUMBER OF SPACE DIMENSIONS')
      !a%nboun=max(a%nboun,0)
      a%gnskew=max(a%gnskew,0)
      exit geomeloop
   end do geomeloop
   
   a%mnode=maxval(a%nnode)

end subroutine countv

