module Mod_nsm_ElasticBoundaryDir
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersElasticBoundaryDir, EndsteElasticBoundaryDir
   
   type, extends(PointerSetter) :: SPElasticBoundaryDir
contains
      procedure :: SpecificSet => SpecificSetElasticBoundaryDir
   end type
   type(SPElasticBoundaryDir) :: SetPointersElasticBoundaryDir
 
   real(rp)                 :: gplev(1)
   integer(ip)              :: elemStatus(1),ngauss_minus,ngauss_plus,ngaus_total
   real(rp), allocatable    :: weigp(:),xloc(:,:)
   
contains

   subroutine SpecificSetElasticBoundaryDir(d)
      implicit none
      class(SPElasticBoundaryDir) :: d
            
      if(a%kfl_ElasticBoundary == 1 )then
         call ConcatenateProcedures(ProcHook_PreDirichlet,nsm_ElasticBoundaryDir)
      end if
   end subroutine  
   
   subroutine nsm_rotrhs(a,e,elrhs)
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
      implicit none
      class(NavierStokesProblem) :: a
      class(FiniteElement)       :: e
      real(rp), intent(inout)    :: elrhs(e%ndime,e%mnode)
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
                  elrhs(1:e%ndime,inode) = matmul(transpose(exnor),elrhs(1:e%ndime,inode))      
               else if(iroty>=1) then                                ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               else if(iroty==-2) then                              ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               end if
            end if
         end do
      end if
   end subroutine nsm_rotrhs
   
   subroutine nsm_ElasticBoundaryDir
      implicit none
      real(rp), pointer :: exnor(:,:) => NULL()
      real(rp)    :: eldisp(e%ndime,e%pnode)
      integer(ip) :: idime,inode,ipoin,ibopo,iroty,jdime,jpoin,jnode
      real(rp)    :: myelvel(e%ndime,e%pnode)
      
      call e%gather(e%ndime,eldisp,a%EB_Displacement)
      myelvel = elvel(:,1:e%pnode,2)
      
      !Rotate displacements if necessary
      call nsm_rotrhs(a,e,eldisp)
      call nsm_rotrhs(a,e,myelvel)
      
      call e%elmlen
      
      !Rotate system
      if(a%kfl_local==1) then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
            if(ibopo>0) then
               iroty=a%kfl_fixrs(ipoin)
               if(iroty==-1) then                                    ! Tangent system 
                  call nsm_rotmat(inode,e%pnode,e%mnode,e%ndime,size(elmat,1),&
                        &      elmat,elrhs,exnor)
               else if(iroty>=1) then                                ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               else if(iroty==-2) then                              ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               end if
            end if
         end do
      endif
      
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         do idime = 1,e%ndime
            if (a%kfl_fixno(idime,ipoin) == 6) then
               elmat(idime,inode,:,:) = 0.0_rp
               elrhs(idime,inode) = 0.0_rp
               
               elmat(idime,inode,idime,inode) = a%EB_density*1+a%EB_Stiffness*a%dtime**2+a%EB_Damping
               elrhs(idime,inode) = myelvel(idime,inode)*a%EB_density -eldisp(idime,inode)*a%EB_Stiffness*a%dtime
               
               if (a%istep > a%EB_NDelaySteps) then
                  if (size(elmat,1)==e%ndime+1) then
                     elmat(idime,inode,e%ndime+1,inode) = -a%dtime
                  else
                     elrhs(idime,inode) = elrhs(idime,inode) + a%dtime*elpre(inode,1)
                  endif
               endif
            endif
         enddo
      enddo
      
      !Undo rotate system
      if(a%kfl_local==1) then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
            if(ibopo>0) then
               iroty=a%kfl_fixrs(ipoin)
               if(iroty==-1) then                                    ! Tangent system 
                  call nsm_rotmatInv(inode,e%pnode,e%mnode,e%ndime,size(elmat,1),&
                        &      elmat,elrhs,exnor)
               else if(iroty>=1) then                                ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               else if(iroty==-2) then                              ! Given system
                  call runend('nsm_elmdir: skcos not ready')
               end if
            end if
         end do
      endif
      
   end subroutine
             
   subroutine EndsteElasticBoundaryDir
      integer(ip) :: ielem, nelem, ndime, npoin, icomp, ipoin, nsteps,idime,ibopo,iroty
      real(rp)    :: prhs(1), dpdt, LHSDtinv,poinveloc(a%ndofn-1),poindisp(a%ndofn-1)
      real(rp), pointer :: exnor(:,:) => NULL()   

      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      do ipoin = 1,npoin
         poinveloc = a%veloc(:,ipoin,1)
         poindisp = a%EB_Displacement(:,ipoin)
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if(ibopo>0) then
            iroty=a%kfl_fixrs(ipoin)
            if(iroty==-1) then                                    ! Tangent system 
               poinveloc = matmul(transpose(exnor),poinveloc)
               poindisp = matmul(transpose(exnor),poindisp)      
            else if(iroty>=1) then                                ! Given system
               call runend('nsm_elmdir: skcos not ready')
            else if(iroty==-2) then                              ! Given system
               call runend('nsm_elmdir: skcos not ready')
            end if
         endif
      
         do idime = 1,ndime
            if (a%kfl_fixno(idime,ipoin) == 6)  poindisp(idime) = poindisp(idime)+poinveloc(idime)*a%dtime
         enddo
         
         if(ibopo>0) then
            iroty=a%kfl_fixrs(ipoin)
            if(iroty==-1) then                                    ! Tangent system 
               poindisp = matmul((exnor),poindisp)      
            else if(iroty>=1) then                                ! Given system
               call runend('nsm_elmdir: skcos not ready')
            else if(iroty==-2) then                              ! Given system
               call runend('nsm_elmdir: skcos not ready')
            end if
         endif
         a%EB_Displacement(:,ipoin) = poindisp
      enddo

   end subroutine
      
end module


