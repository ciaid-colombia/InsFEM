module Mod_nsm_elmdir

contains

   subroutine nsm_elmdir(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      use Mod_php_elmdir
      use Mod_NavierStokes
      implicit none
      class(NavierStokesProblem) :: a
      class(FiniteElement)       :: e
      real(rp), intent(inout)    :: elmat(:,:,:,:), elrhs(:,:)
      integer(ip) :: aux,aux1
      integer(ip) :: suction_kfl_fixno(e%ndime,e%pnode)
   
      aux=a%ndofn
      aux1=a%ndofbcstart
      
      call nsm_pressdir(a,e,a%ndofn,elmat,elrhs)
      call nsm_rotdir(a,e,a%ndofn,elmat,elrhs)
      if (a%kfl_SuctionDirichletBC == 1) call nsm_suctiondir(a,e,suction_kfl_fixno,1)
      call php_elmdir(a,e,a%ndofn,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
      if (a%kfl_SuctionDirichletBC == 1) call nsm_suctiondir(a,e,suction_kfl_fixno,2)
      call nsm_signdir(a,e,a%ndofn,elmat,elrhs)
   end subroutine      
   
   subroutine nsm_suctiondir(a,e,suction_kfl_fixno,itask)
      use typre
      use Mod_Element
      use Mod_NavierStokes
      implicit none
      class(NavierStokesProblem) :: a
      class(FiniteElement)       :: e
      integer(ip) :: itask,suction_kfl_fixno(e%ndime,e%pnode)
      integer(ip) :: idime,inode,PointType
      
      !This subroutine modifies boundary conditions for free surface problems 
      !It prevents boundary conditions to be applied if there is suction (negative pressure)
      !It should avoid fluid stuck in the ceiling
      if (itask == 1) then
         suction_kfl_fixno = a%kfl_fixno(1:e%ndime,e%lnods(1:e%pnode))
         
         do inode = 1,e%pnode
            do idime = 1,e%ndime
               if (a%kfl_fixno(idime,e%lnods(inode)) == 5) then
                  if (a%press(e%lnods(inode),1) > 0) then
                     a%kfl_fixno(idime,e%lnods(inode)) = 1
                  endif
               endif
            enddo
         enddo
      elseif (itask == 2) then
          a%kfl_fixno(1:e%ndime,e%lnods(1:e%pnode)) = suction_kfl_fixno
      endif
   end subroutine
         
   subroutine nsm_pressdir(a,e,ndofn,elmat,elrhs)
      !This subroutine prescribes pressure at a point if the flow is confined
      use typre
      use Mod_Element
      use Mod_NavierStokes
      implicit none
      class(NavierStokesProblem) :: a
      class(FiniteElement)       :: e
      real(rp), intent(inout)    :: elmat(ndofn,e%mnode,ndofn,e%mnode), elrhs(ndofn,e%mnode)
      integer(ip), intent(in)    :: ndofn
      real(rp)    :: adiag
      integer(ip) :: inode,jnode,ipoin,idime,idofn,jdime
      
      !Prescribe one pressure degree of freedom if the flow is confined
      if(a%kfl_confi==1) then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            if (ipoin==a%nodpr) then
               adiag=elmat(ndofn,inode,ndofn,inode)
               elmat(ndofn,inode,1:ndofn,1:e%pnode)=0.0_rp       ! Row
               
               !Column only if kfl_DeleteDirichletColumns
               if (a%kfl_DeleteDirichletColumns) then
                  do jnode=1,e%pnode
                     do jdime=1,ndofn
                        elrhs(jdime,jnode)=elrhs(jdime,jnode)-elmat(jdime,jnode,ndofn,jnode)*a%prepr
                     end do
                  end do
                  elmat(1:ndofn,1:e%pnode,ndofn,inode) = 0.0_rp     ! Column
               endif
               elmat(ndofn,inode,ndofn,inode)=adiag
               elrhs(ndofn,inode)=adiag*a%prepr
            end if
         end do
      end if
   end subroutine nsm_pressdir

   subroutine nsm_rotdir(a,e,ndofn,elmat,elrhs)
      !------------------------------------------------------------------------
      !
      ! This routine modifies the element stiffness matrix of the Navier-
      ! Stokes equations to impose the correct boundary conditions for
      ! local boundary conditions (in skew-systems). In this case the rotation 
      ! of the nodal matrices if this system is the tangent is understood to be  
      !
      !  { NORMAL , TANGENT 1 , TANGENT 2 }
      !
      ! The tangent vectors have been computed according to a good conditioning
      ! of the calculations. Therefore, only the possibility of prescribing
      ! the normal component must be considered when this option is used.
      ! Also, if the flow is confined a pressure value is prescribed.
      !
      !------------------------------------------------------------------------
      use typre
      use Mod_Element
      use Mod_NavierStokes
      implicit none
      class(NavierStokesProblem) :: a
      class(FiniteElement)       :: e
      real(rp), intent(inout)    :: elmat(ndofn,e%mnode,ndofn,e%mnode), elrhs(ndofn,e%mnode)
      integer(ip), intent(in)    :: ndofn
      real(rp)          :: adiag
      integer(ip)       :: inode,jnode,ipoin,idime,idofn,jdime
      integer(ip)       :: iroty,ibopo
      integer(ip)       :: iffix_aux
      real(rp), pointer :: exnor(:,:) => NULL()
      
      !Rotate the nodal matrices, the element nodal velocities and RHS to 
      !prescribe boundary conditions in a skew-system, either the tangent
      !one or another prescribed by the user. 
      if(a%kfl_local==1) then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
            if(ibopo>0) then
               iroty=a%kfl_fixrs(ipoin)
               if(iroty==-1) then                                    ! Tangent system 
                  call nsm_rotmat(inode,e%pnode,e%mnode,e%ndime,ndofn,elmat,elrhs,exnor)
               else if(iroty>=1) then                                ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               else if(iroty==-2) then                              ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               end if
            end if
         end do
      end if
   end subroutine nsm_rotdir

   subroutine nsm_signdir(a,e,ndofn,elmat,elrhs)
      use typre
      use Mod_Element
      use Mod_NavierStokes
      implicit none
      class(NavierStokesProblem) :: a
      class(FiniteElement)       :: e
      real(rp), intent(inout)    :: elmat(ndofn,e%mnode,ndofn,e%mnode), elrhs(ndofn,e%mnode)
      integer(ip), intent(in)    :: ndofn
      real(rp)          :: adiag
      integer(ip)       :: inode,jnode,ipoin,idime,idofn,jdime
      integer(ip)       :: iroty,ibopo
      integer(ip)       :: iffix_aux
      real(rp), pointer :: exnor(:,:) => NULL()
      
      !Canvi de Signe 
      !Modifies the monolithic matrix so that one has G & -Gt (pos def as in faust) 
      !instead of G & Gt (simetric -Zephyr & Fantom)
      !To be left commented out 
      !It is related to some observations Oriol made where the convergence of the 
      !iterative solver was much better with pos def form
      elmat(ndofn,1:e%pnode,1:ndofn,1:e%pnode) = - elmat(ndofn,1:e%pnode,1:ndofn,1:e%pnode)
      elrhs(ndofn,1:e%pnode) = -elrhs(ndofn,1:e%pnode)
   end subroutine nsm_signdir

end module
