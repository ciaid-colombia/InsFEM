submodule(Mod_CauchyElement) CauchyElement_Irreducible

   implicit none

contains

   !---------------- Subroutines for RHS non NR ----------------
   
   module subroutine sld_TimeInt2extF(e,nsteps,Integrator,dtinv2,gpdisp,elext,densi)
      implicit none
      class(FiniteElement)   , intent(in) :: e
      type(TimeIntegratorDt2), intent(in) :: Integrator
      integer(ip) , intent(in)    :: nsteps
      real(rp)    , intent(in)    :: gpdisp(e%ndime,nsteps-1)
      real(rp)    , intent(in)    :: densi,dtinv2
      real(rp)    , intent(inout) :: elext(e%ndime)
      real(rp)                    :: gprhs(e%ndime)
      
      !Time integration
      call Integrator%GetRHS(e%ndime,gpdisp,gprhs)
      elext = elext + gprhs*densi*dtinv2
      
   end subroutine

   module subroutine sld_ComputeExternalForces(e,densi,grnor,gravi,traction,elext,dvolu,force_factor)
      implicit none
      class(FiniteElement) , intent(in) :: e
      real(rp), intent(in) :: densi,grnor,gravi(e%ndime)
      real(rp), intent(in) :: traction(e%ndime),dvolu,force_factor
      real(rp), intent(inout) :: elext(e%ndime)
      
      elext = elext + densi*grnor*gravi + traction
      elext = elext*force_factor

   end subroutine

   module subroutine sld_elmrhu(e,dvolu,elext,elrhs)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the rhs terms 
      !    (v, f) + (v, u_n/dt) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      real(rp),    intent(in)    :: elext(e%ndime)
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)
      integer(ip)                :: inode
      real(rp)                   :: aux(e%ndime)


      aux = elext(:)*dvolu 
      do inode=1,e%pnode
          elrhs(:,inode) = elrhs(:,inode) + e%shape(inode,e%igaus)*aux 
      end do

   end subroutine sld_elmrhu

   !---------------------------- Subroutines for LHS ------------------------------

   module subroutine sld_buildMassMatrix(e,dvol,ndime,densi,elmat,dtinv)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the mass matrix 
      !    Nj*Ni*rho*dv*dtinv_lhs
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip) , intent(in)    :: ndime
      real(rp)    , intent(in)    :: densi,dvol,dtinv
      real(rp)    , intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
      integer(ip)             :: idime,jdime,inode,jnode

      forall(inode=1:e%pnode,jnode=1:e%pnode,idime=1:e%ndime)

              elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) &
                  & + densi*dtinv*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*dvol

      end forall

   end subroutine sld_buildMassMatrix
   

   module subroutine sld_elmk(e,ndime,tn,c_elas,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental material matrix 
      !    B_t*D*B 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)    :: tn,ndime
      real(rp),    intent(in)    :: c_elas(tn,tn),dvol
      real(rp),    intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
      integer(ip)             :: idime,jdime,inode,jnode
      real(rp)                :: B_j(tn,ndime),B_i(ndime,tn)
      real(rp)                :: auxB_i(tn,ndime)
      real(rp)                :: mat_aux1(tn,ndime)
      real(rp)                :: mat_aux2(ndime,ndime)

      do inode=1,e%pnode

          call sld_calculateB(e,ndime,tn,inode,auxB_i)
          B_i   = transpose(auxB_i)

          do jnode=1,e%pnode

              call sld_calculateB(e,ndime,tn,jnode,B_j)

              mat_aux1=matmul(c_elas,B_j)
              mat_aux2=matmul(B_i,mat_aux1)*dvol

              elmat(1:ndime,inode,1:ndime,jnode) = elmat(1:ndime,inode,1:ndime,jnode) + mat_aux2(1:ndime,1:ndime)

          end do
      end do

   end subroutine sld_elmk

   module subroutine sld_elmgeo(e,ndime,tn,stress,dvol,elmat)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the elemental geometric matrix 
      !    B_t*K_geo*B  , see Belytschko for notation, used only in 
      ! non linear scenario
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in)    :: tn,ndime
      real(rp)   , intent(in)    :: stress(tn),dvol
      real(rp)   , intent(inout) :: elmat(ndime,e%pnode,ndime,e%pnode)
      integer(ip)                :: idime,jdime,inode,jnode
      real(rp)                   :: istress(ndime,ndime)
      real(rp)                   :: B_i(ndime),B_j(ndime)
      real(rp)                   :: vec_aux1(ndime)
      real(rp)                   :: res_aux

      call getStressTensor(tn,ndime,stress,istress)

      do inode=1,e%pnode

          call sld_calculateB_inditial(e,ndime,inode,B_i)

          do jnode=1,e%pnode

          call sld_calculateB_inditial(e,ndime,jnode,B_j)

              call matvec(ndime,ndime,istress,B_j,vec_aux1)
              res_aux=dot_product(B_i,vec_aux1)*dvol

              do idime=1,ndime
                  elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + res_aux
              end do

          end do
      end do

   end subroutine sld_elmgeo

   module subroutine getInternalForce(e,ndime,tn,stress,dvol,force)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the element internal force
      !    B_t*Sig  , see Belytschko for notation, used only in 
      ! non linear scenario
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement), intent(in) :: e
      integer(ip), intent(in) :: tn,ndime
      real(rp)   , intent(in) :: stress(tn),dvol
      real(rp)   , intent(inout) :: force(ndime,e%pnode)
      integer(ip)             :: idime,inode
      real(rp)                :: istress(ndime,ndime)
      real(rp)                :: B_i(ndime,tn)
      real(rp)                :: vec_aux1(ndime)

      !call getStressTensor(tn,ndime,stress,istress)

      !do inode=1,e%pnode

      !    call sld_calculateB_inditial(e,ndime,inode,B_i)

      !    call matvec(ndime,ndime,istress,B_i,vec_aux1)

      !    do idime=1,ndime
      !        force(idime,inode) = force(idime,inode) + vec_aux1(idime)*dvol
      !    end do

      !end do

      do inode=1,e%pnode

          call sld_calculateBt(e,ndime,tn,inode,B_i)

          vec_aux1 = matmul(B_i,stress)

          force(:,inode) = force(:,inode) + vec_aux1(:)*dvol

      end do

   end subroutine getInternalForce

end submodule CauchyElement_Irreducible
