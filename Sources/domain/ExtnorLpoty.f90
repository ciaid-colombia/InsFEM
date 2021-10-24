subroutine ExtnorLpoty(a)
!-----------------------------------------------------------------------
! NAME
!    extnor
! DESCRIPTION
!    This routine computes exterior normals and the list of point type (LPOTY)
! OUTPUT
!    NBOPO                   : Number of points located on the boundary
!    EXNOR(NDIME,NDIME,NBOPO): Local basis at boundary nodes,
!    LPOTY(NPOIN)            : defined by   
!    LPOTY(IPOIN) = 0          If IPOIN is interior
!                 = IBOPO      Number of boundary point
!-----------------------------------------------------------------------
   use typre
   use Mod_Setext
   use Mod_Element
   use Mod_Mesh

   implicit none
   
   class(FemMesh) :: a
   
   integer(ip)           :: ielty,istat,lnodf(a%mnode)
   real(rp)              :: zerod
   
   real(rp), allocatable, target :: exwor(:,:)
   integer(ip)           :: ierr
   
   integer(ip)           :: islav,imast,ibopo,iboun,idime,ielem,igaus,inodb,ipoin,ispos,ispos2,linea,pnode,pblty,inode
   logical(lg)           :: slerr
   
   class(FiniteElement), pointer :: e => NULL()
   real(rp)              :: dvol,rpow
   real(rp)              :: auxcartd(a%ndime,a%mnode)
   
   integer(ip) :: auxkfl_alemov = 0
  

  
   call a%Timer%ExtnorLpoty%Tic
   
   a%vodom=0.0_rp            !Domain volume
   !FOR EXNOR and LPOTY 
   call a%Memor%alloc(a%npoin,a%lpoty,'lpoty','exnor')
   call a%Memor%alloc(a%ndime,a%npoin,exwor,'exwor','exnor')
  
   !For the first pass, a%displ has not been set yet
   if ((a%kfl_alemov == 1)  ) then
      if (.not. (associated(a%displ)) ) then
         auxkfl_alemov = 1
         a%kfl_alemov = 0
      elseif (size(a%displ,2) /= a%npoin) then
         !If adaptive this might not be the same size
         auxkfl_alemov = 1
         a%kfl_alemov = 0
      endif
   endif
  
   call a%ElementAlloc(e,a%Memor,'DefaultRule','ExtnorLpoty')
   elements : do ielem = 1,a%nelem
      call a%ElementLoad(ielem,e)
      
      call e%elmdel
      
      gausspoints : do igaus = 1,e%pgaus
         e%igaus = igaus
         
         call e%elmder
         dvol = e%weigp(e%igaus)*e%detjm
         a%vodom = a%vodom + dvol
         !exwor(:,e%lnods(1:e%pnode)) = exwor(:,e%lnods(1:e%pnode)) + e%cartd(:,1:e%pnode)*dvol
         auxcartd(:,1:e%pnode) = e%cartd(:,1:e%pnode)*dvol
         call a%AssemblyToArray(e,a%ndime,auxcartd,exwor)
      enddo gausspoints
   enddo elements
   !call a%EndAssemblyToArray(a%ndime,exwor)
   
   call a%ElementDealloc(e,a%Memor,'DefaultRule','ExtnorLpoty')
   
   !weighting of exwor with vmass
   !weighting is done to avoid too small or too big values of exwor (mesh size independent)
   ! in 2D a good scaling is obtained with rpow=0.5
   ! in 3D a good scaling is obtained with rpow=2.0/3.0
   if (a%ndime==2) then
      rpow = 0.5_rp
   else
      rpow = 2.0_rp/3.0_rp
   end if
   do ipoin=1,a%npoin
      exwor(:,ipoin) = exwor(:,ipoin)*(a%vmass(ipoin)**rpow)
   end do
   
   !Prescribed ExternalNormal
   do ipoin  =1,a%npoinLocal
      !If External Normal systems, set the normal as read 
      if (a%IsExnor(ipoin) .eqv. .true.) then
         exwor(1:a%ndime,ipoin) = a%ExternalNormal(ipoin)%a(1:a%ndime)
      endif
   enddo
   
   call a%ArrayCommunicator%GhostCommunicate(a%ndime,exwor)
   
   !Periodic Boundary Conditions
   if (a%kfl_perio == 1) call a%MasterToSlave(a%ndime,exwor)
   
   !Hanging nodes
   !if (a%kfl_HangingNodes .eqv. .true.) call a%InterpolateHangingValues(a%ndime,exwor)
   
   !Counting of the boundary points
   zerod=1.0e-5_rp
   a%nbopo=0
   do ipoin=1,a%npoin
      dimensions: do idime=1,a%ndime
         if(abs(exwor(idime,ipoin))>zerod) then
            a%nbopo=a%nbopo+1
            a%lpoty(ipoin)=1
            exit dimensions
         end if
      end do dimensions
   end do
   
   !Allocate a%Memory for the vector EXNOR
   call a%Memor%alloc(a%ndime,a%ndime,a%nbopo,a%exnor,'exnor','exnor')
   
   if(a%nbopo>0) then
      !Compute the vectors EXNOR and LPOTY
      call setext(a,zerod,exwor)
   end if

   !Dealloc a%Memory
   call a%Memor%dealloc(a%ndime,a%npoin,exwor,'exwor','exnor')
   
    !undo first pass correction
   if (auxkfl_alemov == 1) a%kfl_alemov = 1

   call a%Timer%ExtnorLpoty%Toc
   
end subroutine ExtnorLpoty


