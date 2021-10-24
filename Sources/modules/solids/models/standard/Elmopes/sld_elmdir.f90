module Mod_sld_elmdir

contains

   subroutine sld_elmdir(a,e,elmat,elrhs)
      use typre
      use Mod_Element
      use Mod_php_elmdir
      use Mod_Solids
      implicit none
      class(SolidsProblem)      :: a
      class(FiniteElement)      :: e
      real(rp), intent(inout)   :: elmat(:,:,:,:), elrhs(:,:)
   
      call sld_rotdir(a,e,a%ndofn,elmat,elrhs)
      call php_elmdir(a,e,a%ndofn,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
   end subroutine      
         
   subroutine sld_rotdir(a,e,ndofn,elmat,elrhs)
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
      use Mod_Solids
      implicit none
      class(SolidsProblem)             :: a
      integer(ip)                            :: ndofn
      class(FiniteElement), intent(in)        :: e
      real(rp), intent(inout)                :: elmat(ndofn,e%mnode,ndofn,e%mnode), elrhs(ndofn,e%mnode)

      real(rp)                   :: adiag
      integer(ip)                :: inode,jnode,ipoin,idime,idofn,jdime
      integer(ip)                :: iroty,ibopo
      integer(ip)                :: iffix_aux
      real(rp), pointer          :: exnor(:,:) => NULL()
      
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
                  call sld_rotmat(inode,e%pnode,e%mnode,e%ndime,ndofn,&
                        &      elmat,elrhs,exnor)
               else if(iroty>=1) then                                ! Given system
                  call runend('sld_elmdir: skcos not ready')
               else if(iroty==-2) then                              ! Given system
                  call runend('sld_elmdir: skcos not ready')
               end if
            end if
         end do
      end if
end subroutine sld_rotdir

end module
