module Mod_NSCompressibleExplicitElement
   use typre
   use Mod_Element
   implicit none

   real(rp), parameter :: zeroc = epsilon(0.0_rp)

contains

   subroutine nsc_ComputeMomentumExternalForces(e,grnor,gravi,elext)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: elext(e%ndime),grnor,gravi(e%ndime)
      
      elext = elext + grnor*gravi

   end subroutine

   subroutine nsc_ComputeEnergyExternalHeat(e,srce,elexh)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp)                   :: srce, elexh
      
      elexh = elexh + srce

   end subroutine

   subroutine nsc_res_unst(e,ndime,gpvar,dtinv,elres)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit unsteady residual 
      ! var - var_n / dt
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: gpvar(ndime,*)
      real(rp),    intent(in)    :: dtinv
      real(rp),    intent(inout) :: elres(ndime)

      real(rp)                   :: aux(ndime)
      integer(ip)                :: inode

      aux = (gpvar(:,1) - gpvar(:,2)) * dtinv
 
      elres = elres + aux 

   end subroutine nsc_res_unst

   subroutine nsc_res_dm(e,divmom,elcd)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit mass-momentum convective residual term
      ! div mom
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: divmom
      real(rp),    intent(inout) :: elcd(1)

      elcd(1) = elcd(1) + divmom

   end subroutine nsc_res_dm

   subroutine nsc_res_md(e,gpvel,mgden,elcm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit momentum-mass residual convective term  
      !  -(mom/rho)*mgden 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: mgden
      real(rp),    intent(inout) :: elcm(e%ndime)
  
      real(rp)                   :: tmp(e%ndime)
  
      tmp = mgden * gpvel

      elcm = elcm - tmp

   end subroutine nsc_res_md

   subroutine nsc_res_mdp(e,val,grden,elcm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit momentum-mass-pressure residual convective term  
      !  (gamma-1)(mom·mom/2 rho^2) *grad rho  
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: val
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(inout) :: elcm(e%ndime)
   
      elcm = elcm + val * grden

   end subroutine nsc_res_mdp

   subroutine nsc_res_mmd(e,gpvel,divmom,elcm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit momentum-momentum-div convective residual term
      ! (mom/rho)*div mom
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: divmom
      real(rp),    intent(inout) :: elcm(e%ndime)
  
      real(rp)                   :: divmom_mom(e%ndime)
  
      divmom_mom = divmom * gpvel
  
      elcm = elcm + divmom_mom

   end subroutine nsc_res_mmd

   subroutine nsc_res_mmg(e,mgmom,elcm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit momentum-momentum-grad convective residual term 
      !   mgmom 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: mgmom(e%ndime)
      real(rp),    intent(inout) :: elcm(e%ndime)
   
      elcm = elcm + mgmom

   end subroutine nsc_res_mmg

   subroutine nsc_res_mmp(e,val,gpvel,grmom,elcm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit momentum-momentum-pressure residual convective term  
      !  -(gamma-1)(mom/rho)·(grad mom)T 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: val
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elcm(e%ndime)
  
      real(rp)                   :: gradmom_mom(e%ndime)
  
      integer(ip)                :: idime
  
      gradmom_mom =  matmul(gpvel,grmom)

      elcm = elcm - val * gradmom_mom

   end subroutine nsc_res_mmp

   subroutine nsc_res_mep(e,val,grene,elcm)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit momentum-energy-pressure residual convective term  
      !  (gamma-1) * grad ene 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: val
      real(rp),    intent(in)    :: grene(e%ndime)
      real(rp),    intent(inout) :: elcm(e%ndime)
   
      elcm = elcm + val * grene

   end subroutine nsc_res_mep

   subroutine nsc_res_mf(e,gpden,elext,elrem)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit residual source term for momentum
    !   - rho f
    !
    !-----------------------------------------------------------------------
    implicit none
 
    class(FiniteElement) :: e
    real(rp),    intent(in)    :: gpden
    real(rp),    intent(in)    :: elext(e%ndime)
    real(rp),    intent(inout) :: elrem(e%ndime)
 
    
       elrem = elrem - gpden*elext

   end subroutine nsc_res_mf

   subroutine nsc_res_edvar(e,gpvar,mgden,elce)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit energy-mass(energy) residual convective term  
      !  -(var/rho)·mgden 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvar
      real(rp),    intent(in)    :: mgden
      real(rp),    intent(inout) :: elce(1)
  
      real(rp)                   :: tmp
  
      tmp = gpvar*mgden

      elce(1) = elce(1) - tmp

   end subroutine nsc_res_edvar

   subroutine nsc_res_edpp(e,val,mgden,elce)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit energy-mass-pressure residual convective term  
      ! (gamma-1)·(mom·mom/2 rho^2)·mgden 
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: val
      real(rp),    intent(in)    :: mgden
      real(rp),    intent(inout) :: elce(1)
  
      elce(1) = elce(1) + val*mgden

   end subroutine nsc_res_edpp

   subroutine nsc_res_emvar(e,gpvar,divmom,elce)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit energy-momentum-energy convective residual term
      ! (var/rho)·div mom
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvar
      real(rp),    intent(in)    :: divmom
      real(rp),    intent(inout) :: elce(1)
  
      real(rp)                   :: divmom_var
  
      divmom_var = divmom * gpvar
  
      elce(1) = elce(1) + divmom_var

   end subroutine nsc_res_emvar

   subroutine nsc_res_empp(e,val,gpvel,mgmom,elce)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit uGradmom convective residual term 
    !   -(gamma-1)·(mom/rho)·mgmom 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: val
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: mgmom(e%ndime)
    real(rp),    intent(inout) :: elce(1)

    real(rp)                   :: tmp

    tmp = dot_product(gpvel,mgmom) 

    elce(1) = elce(1) - val*tmp

   end subroutine nsc_res_empp

   subroutine nsc_res_eep(e,val,mgene,elce)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit energy-energy convective residual term 
    !   (gamma)·mgene 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: val
    real(rp),    intent(in)    :: mgene
    real(rp),    intent(inout) :: elce(1)

    elce(1) = elce(1) + val*mgene

   end subroutine nsc_res_eep

  subroutine nsc_res_ef(e,gpmom,elext,elree)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the gravity source residual term for energy
      !    - f·m
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpmom(e%ndime)
      real(rp),    intent(in)    :: elext(e%ndime)
      real(rp),    intent(inout) :: elree(1)

      elree(1) = elree(1) - dot_product(gpmom,elext)
  
   end subroutine nsc_res_ef
  
   subroutine nsc_res_es(e,gpden,elexh,elree)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the energy source residual term for energy
      !   -rho src
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: gpden
      real(rp),    intent(in)    :: elexh
      real(rp),    intent(inout) :: elree(1)

      elree(1) = elree(1) - gpden*elexh

   end subroutine nsc_res_es

   subroutine nsc_elmvst(e,gpden,gpmom,grden,mgden,grmom,divmom,novst)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit no-viscosity viscous stress tensor
      !    VT_ij = (mu/rho)*(grad + grad^T -(2/3)(div )I)mom
      !    - (mu/rho^2)*(mom_i grad_j + mom_j grad_i -(2/3)(mom_k div_k )I)rho
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e   
      real(rp),    intent(in)    :: gpden,mgden,divmom
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: gpmom(e%ndime)
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(inout) :: novst(e%ndime,e%ndime)

      integer(ip)                :: idime,jdime
      real(rp)                   :: tmp,gpaux(e%ndime)
      real(rp)                   :: divm,divd

      tmp = 1.0_rp/gpden
      gpaux = gpmom/gpden
      divm = (2.0_rp/3.0_rp)*tmp*divmom
      divd = (2.0_rp/3.0_rp)*tmp*mgden

      do idime = 1,e%ndime
         novst(idime,idime) = -divm+divd
         do jdime = 1,e%ndime
            novst(idime,jdime) = novst(idime,jdime) +&
             tmp*(grmom(idime,jdime) + grmom(jdime,idime)) -&
             tmp*(gpmom(idime)*grden(jdime) + gpmom(jdime)*grden(idime))
         enddo
      enddo

   end subroutine nsc_elmvst

   subroutine nsc_MolecularDiff(e,acvis,actco,visten,cndten)
      implicit none
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: acvis,actco
      real(rp),    intent(out)   :: visten(e%ndime,e%ndime),cndten(e%ndime,e%ndime)
      
      integer(ip)                :: idime

      do idime=1,e%ndime
         visten(idime,idime) = acvis
         cndten(idime,idime) = actco
      enddo
      
   end subroutine

   subroutine nsc_dvist_dh(e,gpvel,heden,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      ! -(mom_j/3rho^2) * Hess_jk rho 
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: heden(e%ndime,e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux(e%ndime)
   
      aux = matmul(gpvel,heden)

      dvist = dvist - aux 

   end subroutine nsc_dvist_dh

   subroutine nsc_dvist_dd(e,gpvel,grden,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      !  (2mom_k/rho^3) * grad rho · grad rho  
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvel(e%ndime),grden(e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux

      aux = dot_product(grden,grden)

      dvist = dvist + 2_rp*aux*gpvel

   end subroutine nsc_dvist_dd

   subroutine nsc_dvist_dl(e,gpvel,heden,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      ! -(mom_k/rho^2) * Lap rho 
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: gpvel(e%ndime),heden(e%ndime,e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux
      integer(ip)                :: idime

      aux= 0.0_rp
      do idime=1,e%ndime
         aux = aux + heden(idime,idime)
      enddo

      dvist = dvist - aux*gpvel

   end subroutine nsc_dvist_dl

   subroutine nsc_dvist_dm1(e,auxin,grden,grmom,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      !  -(1/3rho^2) * grad_j rho * grad_k m_j
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: auxin
      real(rp),    intent(in)    :: grden(e%ndime),grmom(e%ndime,e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux(e%ndime)

      aux = matmul(grden,grmom)

      dvist = dvist - auxin*aux

   end subroutine nsc_dvist_dm1

   subroutine nsc_dvist_dm2(e,auxin,grden,grmom,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      !  -(2/rho^2) * grad_j rho * grad_j m_k
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: auxin
      real(rp),    intent(in)    :: grden(e%ndime),grmom(e%ndime,e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux(e%ndime)

      aux = matmul(grmom,grden)

      dvist = dvist - auxin*aux

   end subroutine nsc_dvist_dm2

   subroutine nsc_dvist_mh(e,hemom,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      ! (1/3rho) * Hess_jk m_j 
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: hemom(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux
      integer(ip)                :: idime,jdime

      do idime=1,e%ndime
         aux= 0.0_rp
         do jdime=1,e%ndime
            aux = aux + hemom(jdime,idime,jdime)
         enddo
         dvist(idime) = dvist(idime) + aux
      enddo

   end subroutine nsc_dvist_mh

   subroutine nsc_dvist_ml(e,hemom,dvist)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the  
      ! (1/rho)*Lap m_k 
      ! term in the derivative of the viscous tensor
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: hemom(e%ndime,e%ndime,e%ndime)
      real(rp),    intent(inout) :: dvist(e%ndime)
   
      real(rp)  :: aux
      integer(ip)                :: idime,jdime

      do idime=1,e%ndime
         aux= 0.0_rp
         do jdime=1,e%ndime
            aux = aux + hemom(jdime,jdime,idime)
         enddo
         dvist(idime) = dvist(idime) + aux
      enddo

   end subroutine nsc_dvist_ml

   subroutine nsc_res_diff_ed(e,gpvel,elvst,grden,elree)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit visc dissipation residual term 
    !   (mom·elvst/rho^2)·grad rho 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: gpvel(e%ndime),grden(e%ndime)
    real(rp),    intent(in)    :: elvst(e%ndime,e%ndime)
    real(rp),    intent(inout) :: elree(1)

    real(rp)                   :: tmp(e%ndime),aux
     
    tmp = matmul(gpvel,elvst)
    
    aux = dot_product(tmp,grden) 

    elree(1) = elree(1) + aux

   end subroutine nsc_res_diff_ed

   subroutine nsc_res_diff_em(e,elvst,grmom,elree)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit visc dissipation residual term 
    !  -(elvst/rho):grad mom 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: elvst(e%ndime,e%ndime)
    real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
    real(rp),    intent(inout) :: elree(1)

    real(rp)                   :: aux
    integer(ip)                :: idime,jdime

    aux = 0.0_rp
    do idime=1,e%ndime
         do jdime = 1,e%ndime
            aux = aux + elvst(idime,jdime)*grmom(idime,jdime)
         enddo
    end do

    elree(1) = elree(1) - aux

   end subroutine nsc_res_diff_em

   subroutine nsc_res_diff_ev(e,gpvel,visten,dvist,elree)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit visc dissipation residual term 
    !  - (mom/rho)·dvist 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: visten(e%ndime,e%ndime),dvist(e%ndime)
    real(rp),    intent(inout) :: elree(1)

    real(rp)                   :: aux,rdvist(e%ndime)
     
    rdvist = matmul(visten,dvist)
    aux = dot_product(gpvel,rdvist) 

    elree(1) = elree(1) - aux

   end subroutine nsc_res_diff_ev

   subroutine nsc_dhflx_edd(e,aux,grden,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    ! 1/(rho*cv) *(3 mom·mom/rho^3 - 2s/rho^2) * grad_j rho * grad_j rho
    ! term in the derivative of the heat flux
    ! 
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: aux
    real(rp),    intent(in)    :: grden(e%ndime)
    real(rp),    intent(inout) :: dhflx
   
    integer(ip)                :: idime,jdime
    real(rp)  :: aux_prod

    aux_prod = dot_product(grden,grden)

    dhflx = dhflx + aux*aux_prod

   end subroutine nsc_dhflx_edd

   subroutine nsc_dhflx_edl(e,aux,heden,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    ! 1/(rho*cv) *(s/rho - mom·mom/rho^2) * Lap rho
    ! term in the derivative of the heat flux
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: aux
    real(rp),    intent(in)    :: heden(e%ndime,e%ndime)
    real(rp),    intent(inout) :: dhflx
   
    integer(ip)                :: idime
    real(rp)  :: aux_sum

    aux_sum = 0.0_rp
    do idime=1,e%ndime
       aux_sum = aux_sum + heden(idime,idime)
    end do
    dhflx = dhflx + aux*aux_sum

   end subroutine nsc_dhflx_edl

   subroutine nsc_dhflx_ede(e,aux,grden,grene,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    !  (2/(rho^2*cv)) * (grad_j rho * grad_j s )
    ! term in the derivative of the heat flux
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: aux
    real(rp),    intent(in)    :: grden(e%ndime),grene(e%ndime)
    real(rp),    intent(inout) :: dhflx
   
    real(rp)  :: aux_prod

    aux_prod = dot_product(grden,grene)

    dhflx = dhflx + aux*aux_prod 

   end subroutine nsc_dhflx_ede

   subroutine nsc_dhflx_edm(e,aux,gpvel,grden,grmom,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    !  -(4*mom_l/(rho^3*cv))*(grad_j rho * grad_j m_l)
    ! term in the derivative of the heat flux
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: aux
    real(rp),    intent(in)    :: gpvel(e%ndime),grden(e%ndime)
    real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
    real(rp),    intent(inout) :: dhflx
   
    real(rp)  :: tmp(e%ndime),aux_prod

    tmp = matmul(gpvel,grmom)
    aux_prod = dot_product(tmp,grden)

    dhflx = dhflx - aux*aux_prod

   end subroutine nsc_dhflx_edm

   subroutine nsc_dhflx_emm(e,aux,grmom,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    !  (1/(rho^2*cv))*(grad_j m_l * grad_j m_l)
    ! term in the derivative of the heat flux
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: aux
    real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
    real(rp),    intent(inout) :: dhflx
  
    integer(ip)                :: idime
    real(rp)  :: tmp 

    tmp = 0.0_rp
    do idime=1,e%ndime
       tmp = tmp + dot_product(grmom(idime,:),grmom(idime,:))
    end do

    dhflx = dhflx + aux*tmp

   end subroutine nsc_dhflx_emm

   subroutine nsc_dhflx_emh(e,aux,gpvel,hemom,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    ! mom_l/(rho^2*cv) * Lap mom_l
    ! term in the derivative of the heat flux
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    real(rp),    intent(in)    :: aux
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: hemom(e%ndime,e%ndime,e%ndime)
    real(rp),    intent(inout) :: dhflx
   
    integer(ip)                :: idime,jdime
    real(rp)  :: tmp,aux_sum(e%ndime)


    aux_sum = 0.0_rp
    do idime=1,e%ndime
       do jdime = 1,e%ndime
           aux_sum(idime) = aux_sum(idime) + hemom(jdime,jdime,idime)
       enddo
    end do

    tmp = dot_product(aux_sum,gpvel)

    dhflx = dhflx + aux*tmp

   end subroutine nsc_dhflx_emh

   subroutine nsc_dhflx_eeh(e,actcn,heene,dhflx)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the  
    ! -1/(rho*cv) * Lap ene
    ! term in the derivative of the heat flux
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: actcn
    real(rp),    intent(in)    :: heene(e%ndime,e%ndime)
    real(rp),    intent(inout) :: dhflx
   
    integer(ip)                :: idime,jdime
    real(rp)  :: aux_sum

    aux_sum = 0.0_rp
    do idime=1,e%ndime
           aux_sum = aux_sum + heene(idime,idime)
    end do

    dhflx = dhflx - actcn*aux_sum

   end subroutine nsc_dhflx_eeh

   subroutine nsc_elmc(e,ndime,dvolu,elc,elr)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit convective terms
    !  - (v, vector)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    integer(ip), intent(in)    :: ndime
    real(rp),    intent(in)    :: elc(ndime)
    real(rp),    intent(inout) :: elr(ndime,e%pnode)

    integer(ip)                :: idime

    do idime=1,ndime
    elr(idime,1:e%pnode) = elr(idime,1:e%pnode) - e%shape(1:e%pnode,e%igaus)*elc(idime)*dvolu
    end do

   end subroutine nsc_elmc

   subroutine nsc_elmc1(e,dvolu,elc,elr)
    !-----------------------------------------------------------------------
    !
    ! This routine computes
    !  - (v, scalar)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: elc
    real(rp),    intent(inout) :: elr(1,e%pnode)


    elr(1,1:e%pnode) = elr(1,1:e%pnode) - e%shape(1:e%pnode,e%igaus)*elc*dvolu

   end subroutine nsc_elmc1

   subroutine nsc_elmrmn_mf(e,dvolu,gpden,elext,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the source rhs terms for momentum
    !   + (v, rho f)  
    !
    !-----------------------------------------------------------------------
    implicit none
 
    class(FiniteElement) :: e
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: gpden
    real(rp),    intent(in)    :: elext(e%ndime)
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)
 
    integer(ip)                :: idime
    
    do idime=1,e%ndime
       elrmn(idime,1:e%pnode) = elrmn(idime,1:e%pnode) + &
                e%shape(1:e%pnode,e%igaus)*dvolu*gpden*elext(idime)
    end do

   end subroutine nsc_elmrmn_mf
   
   subroutine nsc_elmreg_ef(e,dvolu,gpmom,elext,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the gravity source term for energy
      !    (g, rho f.u)
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: gpmom(e%ndime)
      real(rp),    intent(in)    :: elext(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

      real(rp)                   :: tmp

      tmp = dot_product(gpmom,elext)
      elreg = elreg + e%shape(:,e%igaus)*dvolu*tmp

   end subroutine nsc_elmreg_ef

   subroutine nsc_elmreg_es(e,dvolu,gpden,elexh,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the energy source term for energy
      !    (g, rho src) 
      !
      !-----------------------------------------------------------------------
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: gpden
      real(rp),    intent(in)    :: elexh
      real(rp),    intent(inout) :: elreg(e%pnode)

      elreg = elreg + e%shape(:,e%igaus)*dvolu*gpden*elexh

   end subroutine nsc_elmreg_es

   subroutine nsc_elmgrt_i(e,actcn,gpden,gpvel,gpene,grden,grmom,grene,grtem)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit temperature gradient ideal state law
      !    grad tem = - (1/rho*cv) * (grad ene - (ene/rho)*grad rho
      !    - (1/rho)*(mom_j grad_i mom_j + (mom_k mom_k/rho^2)*grad rho)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e   
      real(rp),    intent(in)    :: actcn,gpden,gpene
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: grene(e%ndime)
      real(rp),    intent(inout) :: grtem(e%ndime)

      real(rp)                   :: tmp,gpaux(e%ndime)
      real(rp)                   :: aux,aux1(e%ndime)

      tmp = dot_product(gpvel,gpvel)
      aux = tmp - (gpene/gpden)
      gpaux = matmul(gpvel,grmom)

      aux1 =  (aux*grden) - gpaux + grene
      
      grtem = actcn * aux1

   end subroutine nsc_elmgrt_i

   subroutine nsc_elmr_dif(e,ndime,dvolu,grvar,elr)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit weak form of the diffusive term
      !    -(grad v, grad u)
      !
      !-----------------------------------------------------------------------

      implicit none
      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvolu
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: grvar(ndime,e%ndime)
      real(rp),    intent(inout) :: elr(ndime,e%pnode)

      real(rp)                   :: aux(e%pnode)
      integer(ip)                :: idime

      do idime=1,ndime
      aux = matmul(grvar(idime,:),e%cartd)
      elr(idime,:) = elr(idime,:) - dvolu*aux(:)
      end do

   end subroutine nsc_elmr_dif

   subroutine nsc_elmrdd_dm(e,dvolu,testcm,elrem,elrdd)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the explicit stabilization term for mass
      !  - (grad_i w, tau2*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolu
      real(rp),    intent(in)    :: testcm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrdd(1,e%pnode)

      integer(ip)                :: idime
      real(rp)                   :: aux(e%pnode)
 
      aux = 0.0_rp      

      do idime=1,e%ndime
         aux(:) = aux(:) + testcm(idime,:) * elrem(idime)
      end do
      
      elrdd(1,:) = elrdd(1,:) - aux(:)*dvolu

   end subroutine nsc_elmrdd_dm

   subroutine nsc_elmrmn_md(e,dvolu,gpvel,testfd,elred,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization term for momentum
    !  -(-((mom/rho)·grad v)·(mom/rho)tau1*R(rho)) 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: testfd(e%pnode)
    real(rp),    intent(in)    :: elred
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    integer(ip)                :: idime

    do idime=1,e%ndime
      elrmn(idime,:) = elrmn(idime,:) + testfd(:)*gpvel(idime)*elred*dvolu
    end do

   end subroutine nsc_elmrmn_md

   subroutine nsc_elmrmn_mdp(e,dvolu,val,testcd,elred,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization term for momentum
    !  - ((gamma-1)(mom·mom)/2 rho^2)(div v)tau1*R(rho))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu,val
    real(rp),    intent(in)    :: testcd(e%ndime,e%pnode)
    real(rp),    intent(in)    :: elred
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    integer(ip)                :: idime

    do idime=1,e%ndime
      elrmn(idime,:) = elrmn(idime,:) - val*testcd(idime,:)*elred*dvolu
    end do

   end subroutine nsc_elmrmn_mdp

   subroutine nsc_elmrmn_mmd(e,dvolu,gpvel,testcm,elrem,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for momentum
    !  - ((tau2*R(mom_k)·grad_k v_i)mom_i/rho) 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: testcm(e%ndime,e%pnode)
    real(rp),    intent(in)    :: elrem(e%ndime)
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    real(rp)                   :: aux(e%pnode)
    integer(ip)                :: idime

    aux = matmul(elrem,testcm)

    do idime=1,e%ndime
       elrmn(idime,:) = elrmn(idime,:) - aux(:)*gpvel(idime)*dvolu
    end do

   end subroutine nsc_elmrmn_mmd

   subroutine nsc_elmrmn_mmg(e,dvolu,testfm,elrem,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for momentum
    !  - (((mom/rho)·grad v_i) tau2*R(mom_i)) 
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: testfm(e%pnode)
    real(rp),    intent(in)    :: elrem(e%ndime)
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    integer(ip)                :: idime

    do idime=1,e%ndime
       elrmn(idime,:) = elrmn(idime,:) - testfm(:)*elrem(idime)*dvolu
    end do

   end subroutine nsc_elmrmn_mmg

   subroutine nsc_elmrmn_mmp(e,dvolu,val,gpvel,testcm,elrem,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for momentum
    !  - (-(gamma-1)(mom/rho)·tau2*R(mom)(div v))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu,val
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: testcm(e%ndime,e%pnode)
    real(rp),    intent(in)    :: elrem(e%ndime)
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    real(rp)                   :: aux
    integer(ip)                :: idime

    aux = val * dot_product(gpvel,elrem)
    do idime=1,e%ndime
       elrmn(idime,:) = elrmn(idime,:) + aux*testcm(idime,:)*dvolu
    end do

   end subroutine nsc_elmrmn_mmp

   subroutine nsc_elmrmn_mep(e,dvolu,val,testce,elree,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for momentum
    !  - ((gamma-1)(div v) tau3*R(s))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu,val
    real(rp),    intent(in)    :: testce(e%ndime,e%pnode)
    real(rp),    intent(in)    :: elree
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    integer(ip)                :: idime

    do idime=1,e%ndime
      elrmn(idime,:) = elrmn(idime,:) - val*testce(idime,:)*elree*dvolu
    end do

   end subroutine nsc_elmrmn_mep

   subroutine nsc_elmreg_edvar(e,dvolu,gpvar,testfd,elred,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for energy
    !  - ((-var/rho)((mom/rho)·grad g)tau1*R(rho))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e

    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: gpvar
    real(rp),    intent(in)    :: testfd(e%pnode)
    real(rp),    intent(in)    :: elred
    real(rp),    intent(inout) :: elreg(e%pnode)

    real(rp)                   :: aux1(e%pnode)

    aux1 = gpvar*elred*testfd

    elreg = elreg + aux1*dvolu

   end subroutine nsc_elmreg_edvar

   subroutine nsc_elmreg_edpp(e,dvolu,val,testfd,elred,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for energy
    !  - ((gamma-1)(mom·mom)/2 rho^2)((mom/rho)·grad g)tau1*R(rho))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e

    real(rp),    intent(in)    :: dvolu,val
    real(rp),    intent(in)    :: testfd(e%pnode)
    real(rp),    intent(in)    :: elred
    real(rp),    intent(inout) :: elreg(e%pnode)

    real(rp)                   :: aux1(e%pnode)

    aux1 = val*elred*testfd

    elreg = elreg - aux1*dvolu

   end subroutine nsc_elmreg_edpp

   subroutine nsc_elmreg_emvar(e,dvolu,gpvar,testcm,elrem,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for energy
    !  - ((tau2*R(mom_i)·grad_i g)var/rho)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: gpvar
    real(rp),    intent(in)    :: testcm(e%ndime,e%pnode)
    real(rp),    intent(in)    :: elrem(e%ndime)
    real(rp),    intent(inout) :: elreg(e%pnode)

    real(rp)                   :: aux(e%pnode)

    aux = matmul(elrem,testcm)

    elreg = elreg - aux*gpvar*dvolu

   end subroutine nsc_elmreg_emvar

   subroutine nsc_elmreg_empp(e,dvolu,val,gpvel,testfm,elrem,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for energy
    !  - (-(gamma-1)((mom/rho^2)·grad g)tau2*R(mom_i)·mom_i)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu,val
    real(rp),    intent(in)    :: gpvel(e%ndime)
    real(rp),    intent(in)    :: testfm(e%pnode)
    real(rp),    intent(in)    :: elrem(e%ndime)
    real(rp),    intent(inout) :: elreg(e%pnode)

    real(rp)                   :: aux

    aux = dot_product(elrem,gpvel)

    elreg = elreg + val*testfm*aux*dvolu

   end subroutine nsc_elmreg_empp

   subroutine nsc_elmreg_eep(e,dvolu,val,testfe,elree,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the explicit stabilization terms for energy
    !  - (gamma((mom/rho)·grad g)tau3*R(s))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: val
    real(rp),    intent(in)    :: testfe(e%pnode)
    real(rp),    intent(in)    :: elree
    real(rp),    intent(inout) :: elreg(e%pnode)

    elreg = elreg - val*testfe*elree*dvolu

   end subroutine nsc_elmreg_eep

   subroutine nsc_elmrmn_mdf(e,dvolstab,testd,elext,elred,elrmn)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the source stabilization terms for momentum
    !  - (elext testd, tau1*R(rho))
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: dvolstab,elred
    real(rp),    intent(in)    :: testd(e%pnode)
    real(rp),    intent(in)    :: elext(e%ndime)
    real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

    real(rp)                   :: aux(e%pnode)
    integer(ip)                :: idime

    do idime=1,e%ndime
      aux(:) = elext(idime)*testd(:)*elred
      elrmn(idime,:) = elrmn(idime,:) - aux(:)*dvolstab
    end do

   end subroutine nsc_elmrmn_mdf

   subroutine nsc_elmreg_edf(e,dvolstab,testg,elexh,elred,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the source estabilization terms for energy
    !   - (elexh testg, tau1*R(rho))  
    !
    !-----------------------------------------------------------------------
    implicit none
 
    class(FiniteElement) :: e
    real(rp),    intent(in)    :: dvolstab,elexh,elred
    real(rp),    intent(in)    :: testg(e%pnode)
    real(rp),    intent(inout) :: elreg(e%pnode)
 
    elreg = elreg - elexh*testg*elred*dvolstab

   end subroutine nsc_elmreg_edf

   subroutine nsc_elmreg_emf(e,dvolstab,testg,elext,elrem,elreg)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the source estabilization terms for energy
    !   - (elext testg, tau2*R(mom))  
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e

    real(rp),    intent(in)    :: dvolstab
    real(rp),    intent(in)    :: testg(e%pnode)
    real(rp),    intent(in)    :: elext(e%ndime), elrem(e%ndime)
    real(rp),    intent(inout) :: elreg(e%pnode)

    real(rp)                   :: aux
    integer(ip)                :: idime

    aux = dot_product(elext,elrem)

    elreg = elreg - aux*testg*dvolstab

   end subroutine nsc_elmreg_emf

   subroutine nsc_hessian(e,ndime,elvar,hevar)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the Hessian of a variable
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    integer(ip), intent(in)    :: ndime
    real(rp),    intent(in)    :: elvar(ndime,e%pnode)
    real(rp),    intent(inout) :: hevar(e%ndime,e%ndime,ndime)

    integer(ip) :: idime,itens
    real(rp)    :: tmphe(e%ntens)

    do idime=1,ndime
         do itens = 1,e%ntens
            tmphe(itens) = dot_product(e%hessi(itens,1:e%pnode),elvar(idime,1:e%pnode))
         enddo
         if(e%ndime.eq.2) then
           hevar(1,1,idime)=tmphe(1)
           hevar(1,2,idime)=tmphe(3)
           hevar(2,1,idime)=tmphe(3)
           hevar(2,2,idime)=tmphe(2)
         else
           hevar(1,1,idime)=tmphe(1)
           hevar(1,2,idime)=tmphe(4)
           hevar(2,1,idime)=tmphe(4)
           hevar(2,2,idime)=tmphe(2)
           hevar(1,3,idime)=tmphe(5)
           hevar(3,1,idime)=tmphe(5)
           hevar(2,3,idime)=tmphe(6)
           hevar(3,2,idime)=tmphe(6)
           hevar(3,3,idime)=tmphe(3)
         end if
    end do

   end subroutine nsc_hessian

   subroutine nsc_elmrm_nmdd1(e,dvolstab,acvis,gpden,mgden,testm,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (grad_i N_i (2*mu*mom_j/rho^3)*grad_j rho, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden,mgden
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

      real(rp)                   :: aux
 
      aux = (2.0_rp*acvis*mgden*elred)/(gpden*gpden) 

      elrmn = elrmn - aux*testm*dvolstab

   end subroutine nsc_elmrm_nmdd1

   subroutine nsc_elmrm_nmd1(e,dvolstab,acvis,gpden,divmom,testm,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-grad_i N_i (mu/rho^2)*div mom, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden,divmom
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

      real(rp)                   :: aux
 
      aux = (acvis*divmom*elred)/(gpden*gpden) 

      elrmn = elrmn + aux*testm*dvolstab

   end subroutine nsc_elmrm_nmd1

   subroutine nsc_elmrm_nmdd2(e,dvolstab,acvis,gpvel,grden,testm,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (grad_j N_i (2*mu*mom_i/rho^3)*grad_j rho, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime),grden(e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)
      integer(ip)                :: idime

      aux = 2.0_rp*acvis*elred
      auxprod = matmul(grden,testm)

      do idime=1,e%ndime
        elrmn(idime,:) = elrmn(idime,:) - aux*auxprod(:)*gpvel(idime)*dvolstab
      end do

   end subroutine nsc_elmrm_nmdd2

   subroutine nsc_elmrm_nmd2(e,dvolstab,acvis,gpden,grmom,testm,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-grad_j N_i (mu/rho^2)*grad_j mom_i, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%ndime,e%pnode)

      aux = (acvis*elred)/(gpden*gpden) 
      auxprod = matmul(grmom,testm)

      elrmn = elrmn + aux*auxprod*dvolstab

   end subroutine nsc_elmrm_nmd2

   subroutine nsc_elmrm_nmdd3(e,dvolstab,acvis,gpvel,grden,testm,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-grad_j N_i (4*mu*mom_j/3*rho^3)*grad_i rho, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime),grden(e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)
      integer(ip)                :: idime

      aux = 4.0_rp*acvis*elred/3.0_rp 
      auxprod = matmul(gpvel,testm)
 
      do idime=1,e%ndime
        elrmn(idime,:) = elrmn(idime,:) +aux*auxprod(:)*grden(idime)*dvolstab
      end do

   end subroutine nsc_elmrm_nmdd3

   subroutine nsc_elmrm_nmd3(e,dvolstab,acvis,gpden,grmom,testm,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (grad_j N_i (2*mu/3*rho^2)*grad_i mom_j, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%ndime,e%pnode)

      aux = (2.0_rp*acvis*elred)/(3.0_rp*gpden*gpden)
      auxprod = matmul(transpose(grmom),testm)

      elrmn = elrmn - aux*auxprod*dvolstab

   end subroutine nsc_elmrm_nmd3

   subroutine nsc_elmrm_hnmd(e,dvolstab,acvis,gpvel,testh,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-hess_ij N_i (mu/3*rho^2)*mom_j, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: testh(e%ndime,e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)
      integer(ip)                :: idime

      aux = acvis*elred/3.0_rp

      do idime=1,e%ndime
         auxprod = matmul(gpvel,testh(idime,:,:))
         elrmn(idime,:) = elrmn(idime,:) + aux*auxprod(:)*dvolstab
      end do

   end subroutine nsc_elmrm_hnmd

   subroutine nsc_elmrm_lnmd(e,dvolstab,acvis,gpvel,testh,elred,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-hess_jj N_i (mu/rho^2)*mom_i, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: testh(e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux
      integer(ip)                :: idime

      aux = acvis*elred

      do idime=1,e%ndime
         elrmn(idime,:) = elrmn(idime,:) + aux*testh(:)*gpvel(idime)*dvolstab
      end do

   end subroutine nsc_elmrm_lnmd

   subroutine nsc_elmrm_ndm1(e,dvolstab,acvis,gpden,grden,testm,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-grad_i N_i (mu/rho^2)*grad_j rho, tau2*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod

      auxprod = dot_product(grden,elrem)
      aux = (acvis*auxprod)/(gpden*gpden)

      elrmn = elrmn + aux*testm*dvolstab

   end subroutine nsc_elmrm_ndm1

   subroutine nsc_elmrm_ndm2(e,dvolstab,acvis,gpden,grden,testm,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (-grad_j N_i (mu/rho^2)*grad_j rho, tau2*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)
      integer(ip)                :: idime

      aux = acvis/(gpden*gpden)
      auxprod = matmul(grden,testm)

      do idime=1,e%ndime
      elrmn(idime,:) = elrmn(idime,:) + aux*auxprod(:)*elrem(idime)*dvolstab
      end do

   end subroutine nsc_elmrm_ndm2

   subroutine nsc_elmrm_ndm3(e,dvolstab,acvis,gpden,grden,testm,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (grad_j N_i (2*mu/3*rho^2)*grad_i rho, tau2*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: testm(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)
      integer(ip)                :: idime

      aux = (2.0_rp*acvis)/(3.0_rp*gpden*gpden)
      auxprod = matmul(elrem,testm)

      do idime=1,e%ndime
      elrmn(idime,:) = elrmn(idime,:) - aux*auxprod(:)*grden(idime)*dvolstab
      end do

   end subroutine nsc_elmrm_ndm3

   subroutine nsc_elmrm_hnm(e,dvolstab,acvis,gpden,testh,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (hess_ij N_i (mu/3*rho), tau2*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: testh(e%ndime,e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)
      integer(ip)                :: idime

      aux = acvis/(3.0_rp*gpden)

      do idime=1,e%ndime
         auxprod = matmul(elrem,testh(idime,:,:))
         elrmn(idime,:) = elrmn(idime,:) - aux*auxprod(:)*dvolstab
      end do

   end subroutine nsc_elmrm_hnm

   subroutine nsc_elmrm_lnm(e,dvolstab,acvis,gpden,testh,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for momentum
      !  - (hess_jj N_i (mu/rho), tau2*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: testh(e%pnode)
      real(rp),    intent(in)    :: elrem(e%pnode)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux
      integer(ip)                :: idime

      aux = acvis/gpden

      do idime=1,e%ndime
         elrmn(idime,:) = elrmn(idime,:) - aux*testh(:)*elrem(idime)*dvolstab
      end do

   end subroutine nsc_elmrm_lnm

   subroutine nsc_elmre_gmmdd1(e,dvolstab,acdiff,aux_k,grden,teste,elred,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (grad_i G (diff*scalar/rho^2)*grad_i rho, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acdiff,aux_k
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = acdiff*aux_k*elred
      auxprod = matmul(grden,teste)

      elreg = elreg - aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmmdd1

   subroutine nsc_elmre_gmmd1(e,dvolstab,acdiff,gpvel,grmom,teste,elred,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (diff*mom_j/rho^3)*grad_i mom_j, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acdiff
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%ndime)
      real(rp)                   :: auxprodtest(e%pnode)

      aux = acdiff*elred
      auxprod = matmul(gpvel,grmom)
      auxprodtest = matmul(auxprod,teste)

      elreg = elreg + aux*auxprodtest*dvolstab

   end subroutine nsc_elmre_gmmd1

   subroutine nsc_elmre_gmmd2(e,dvolstab,acvis,gpvel,scalar,teste,elred,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (grad_i G (mu*mom_i/rho^2)*scalar, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: scalar,gpvel(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = acvis*scalar*elred
      auxprod = matmul(gpvel,teste)

      elreg = elreg - aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmmd2

   subroutine nsc_elmre_gmmd3(e,dvolstab,acvis,gpden,mgmom,teste,elred,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (mu*mom_j/3*rho^3)*grad_j mom_i, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: mgmom(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = (acvis*elred)/(3.0_rp*gpden*gpden)
      auxprod = matmul(mgmom,teste)

      elreg = elreg + aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmmd3

   subroutine nsc_elmre_ged(e,dvolstab,lambcv,gpden,grene,teste,elred,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (lambda/cv*rho^2)*grad_i ene, tau1*R(rho))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,lambcv,gpden
      real(rp),    intent(in)    :: grene(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elred
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = (lambcv*elred)/(gpden*gpden)
      auxprod = matmul(grene,teste)

      elreg = elreg + aux*auxprod*dvolstab

   end subroutine nsc_elmre_ged

   subroutine nsc_elmre_gmdm1(e,dvolstab,acvis,gpvel,grden,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (2*mu*mom_j/rho^3)*grad_i rho, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1,auxprod(e%pnode)

      aux = 2.0_rp*acvis
      aux1 = dot_product(gpvel,elrem)
      auxprod = matmul(grden,teste)

      elreg = elreg + aux*aux1*auxprod*dvolstab

   end subroutine nsc_elmre_gmdm1

   subroutine nsc_elmre_gmm1(e,dvolstab,acvis,gpden,grmom,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (grad_i G (mu/rho^2)*grad_i mom_j, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1(e%ndime),auxprod(e%pnode)

      aux = acvis/(gpden*gpden)
      aux1 = matmul(elrem,grmom)
      auxprod = matmul(aux1,teste)

      elreg = elreg - aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmm1

   subroutine nsc_elmre_gmdm2(e,dvolstab,acvis,gpvel,grden,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (2*mu*mom_i/rho^3)*grad_j rho, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1,auxprod(e%pnode)

      aux = 2.0_rp*acvis
      aux1 = dot_product(grden,elrem)
      auxprod = matmul(gpvel,teste)

      elreg = elreg + aux*aux1*auxprod*dvolstab

   end subroutine nsc_elmre_gmdm2

   subroutine nsc_elmre_gmm2(e,dvolstab,acvis,gpden,grmom,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (grad_i G (mu/rho^2)*grad_j mom_i, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1(e%ndime),auxprod(e%pnode)

      aux = acvis/(gpden*gpden)
      aux1 = matmul(grmom,elrem)
      auxprod = matmul(aux1,teste)

      elreg = elreg - aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmm2

   subroutine nsc_elmre_gmdm3(e,dvolstab,acvis,gpden,mgden,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (grad_i G (4*mu*mom_j/3*rho^3)*grad_j rho, tau1*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: mgden
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = (4.0_rp*acvis*mgden)/(3.0_rp*gpden*gpden*gpden)
      auxprod = matmul(elrem,teste)

      elreg = elreg - aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmdm3

   subroutine nsc_elmre_gmm3(e,dvolstab,acvis,gpden,divmom,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (2*mu/3*rho^2)*grad_j mom_j, tau1*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis,gpden
      real(rp),    intent(in)    :: divmom
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = (2.0_rp*acvis*divmom)/(3.0_rp*gpden*gpden)
      auxprod = matmul(elrem,teste)

      elreg = elreg + aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmm3

   subroutine nsc_elmre_gmdm4(e,dvolstab,lambcv,gpvel,grden,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (grad_i G (2*lambda*mom_j/cv*rho^3)*grad_i rho, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,lambcv
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1,auxprod(e%pnode)

      aux = 2.0_rp*lambcv
      aux1 = dot_product(gpvel,elrem)
      auxprod = matmul(grden,teste)

      elreg = elreg - aux*aux1*auxprod*dvolstab

   end subroutine nsc_elmre_gmdm4

   subroutine nsc_elmre_gmm4(e,dvolstab,lambcv,gpden,grmom,teste,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (lambda/cv*rho^2)*grad_i mom_j, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,lambcv,gpden
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1(e%ndime),auxprod(e%pnode)

      aux = lambcv/(gpden*gpden)
      aux1 = matmul(elrem,grmom)
      auxprod = matmul(aux1,teste)

      elreg = elreg + aux*auxprod*dvolstab

   end subroutine nsc_elmre_gmm4

   subroutine nsc_elmre_lgmm1(e,dvolstab,acvis,gpvel,testl,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (hess_jj G (mu/rho^2)*mom_i, tau1*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: testl(e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: auxprod

      auxprod = dot_product(elrem,gpvel)

      elreg = elreg - acvis*auxprod*testl*dvolstab

   end subroutine nsc_elmre_lgmm1

   subroutine nsc_elmre_hgmm1(e,dvolstab,acvis,gpvel,testh,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (hess_ij G (mu/rho^2)*mom_j, tau1*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: testh(e%ndime,e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1(e%ndime,e%pnode)
      real(rp)                   :: auxprod(e%pnode)
      integer(ip)                :: idime

      do idime=1,e%ndime
         aux1(idime,:) = matmul(gpvel,testh(idime,:,:))
      end do
      auxprod = matmul(elrem,aux1) 

      elreg = elreg - acvis*auxprod*dvolstab

   end subroutine nsc_elmre_hgmm1

   subroutine nsc_elmre_hgmm2(e,dvolstab,acvis,gpvel,testh,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-hess_ij G (2*mu/3*rho^2)*mom_i, tau1*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,acvis
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: testh(e%ndime,e%ndime,e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,aux1(e%ndime,e%pnode)
      real(rp)                   :: auxprod(e%pnode)
      integer(ip)                :: idime

      aux = 2.0_rp*acvis/3.0_rp
      do idime=1,e%ndime
         aux1(idime,:) = matmul(elrem,testh(idime,:,:))
      end do
      auxprod = matmul(gpvel,aux1) 

      elreg = elreg + aux*auxprod*dvolstab

   end subroutine nsc_elmre_hgmm2

   subroutine nsc_elmre_lgmm2(e,dvolstab,lambcv,gpvel,testl,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-hess_jj G (lambda/cv*rho^2)*mom_i, tau1*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,lambcv
      real(rp),    intent(in)    :: gpvel(e%ndime)
      real(rp),    intent(in)    :: testl(e%pnode)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: auxprod

      auxprod = dot_product(elrem,gpvel)

      elreg = elreg + lambcv*auxprod*testl*dvolstab

   end subroutine nsc_elmre_lgmm2

   subroutine nsc_elmre_lgme(e,dvolstab,lambcv,gpden,testl,elree,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (hess_jj G (lambda/cv*rho), tau1*R(s))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,lambcv,gpden
      real(rp),    intent(in)    :: testl(e%pnode)
      real(rp),    intent(in)    :: elree
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux

      aux = (lambcv*elree)/gpden

      elreg = elreg - aux*testl*dvolstab

   end subroutine nsc_elmre_lgme

   subroutine nsc_elmre_gde(e,dvolstab,lambcv,gpden,grden,teste,elree,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (-grad_i G (lambda/cv*rho)*grad_i rho, tau1*R(e))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,lambcv,gpden
      real(rp),    intent(in)    :: grden(e%ndime)
      real(rp),    intent(in)    :: teste(e%ndime,e%pnode)
      real(rp),    intent(in)    :: elree
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux,auxprod(e%pnode)

      aux = (lambcv*elree)/gpden
      auxprod = matmul(grden,teste)

      elreg = elreg + aux*auxprod*dvolstab

   end subroutine nsc_elmre_gde

   subroutine nsc_testmomentum(e,gpmom,testm)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the projection into the Finite element space
    !  (v, var)
    !
    !-----------------------------------------------------------------------
    implicit none
    class(FiniteElement) :: e
    
    real(rp),    intent(in)    :: gpmom(e%ndime)
    real(rp),    intent(inout) :: testm(e%ndime,e%pnode)

    integer(ip)                :: idime

    do idime=1,e%ndime
       testm(idime,:) = e%shape(:,e%igaus)*gpmom(idime)
    end do

   end subroutine nsc_testmomentum

   subroutine nsc_elmre_gsca(e,dvolstab,aux,scalar,elres,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (G aux*scalar, tau*Res)
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,aux
      real(rp),    intent(in)    :: scalar
      real(rp),    intent(in)    :: elres
      real(rp),    intent(inout) :: elreg(e%pnode)

      elreg = elreg - aux*scalar*elres*e%shape(:,e%igaus)*dvolstab

   end subroutine nsc_elmre_gsca

   subroutine nsc_elmre_ggr(e,dvolstab,aux,grene,elrem,elreg)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - ( G aux*grad_i ene, tau1*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,aux
      real(rp),    intent(in)    :: grene(e%ndime)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elreg(e%pnode)

 
      real(rp)                   :: aux1

      aux1 = dot_product(grene,elrem)
      
      elreg = elreg - aux*aux1*e%shape(:,e%igaus)*dvolstab

   end subroutine nsc_elmre_ggr

   subroutine nsc_elmrm_ngm1(e,dvolstab,aux,grmom,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (N_i aux*grad_j mom_i, tau2*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,aux
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux1(e%ndime)
      integer(ip)                :: idime

      aux1 = matmul(grmom,elrem)

      do idime=1,e%ndime
         elrmn(idime,:) = elrmn(idime,:) &
		          - aux*aux1(idime)*e%shape(:,e%igaus)*dvolstab
      end do

   end subroutine nsc_elmrm_ngm1

   subroutine nsc_elmrm_ngm2(e,dvolstab,aux,divmom,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (N_i aux*grad_j mom_j, tau2*R(mom_i))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,aux
      real(rp),    intent(in)    :: divmom
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

      integer(ip)                :: idime

      do idime=1,e%ndime
         elrmn(idime,1:e%pnode) = elrmn(idime,1:e%pnode) &
		          - aux*divmom*elrem(idime)*e%shape(1:e%pnode,e%igaus)*dvolstab
      end do

   end subroutine nsc_elmrm_ngm2

   subroutine nsc_elmrm_ngm3(e,dvolstab,aux,grmom,elrem,elrmn)
      !-----------------------------------------------------------------------
      !
      ! This routine computes the diffusion stabilization term for energy
      !  - (N_i aux*grad_i mom_j, tau2*R(mom_j))
      !
      !-----------------------------------------------------------------------
      implicit none
      class(FiniteElement) :: e
      
      real(rp),    intent(in)    :: dvolstab,aux
      real(rp),    intent(in)    :: grmom(e%ndime,e%ndime)
      real(rp),    intent(in)    :: elrem(e%ndime)
      real(rp),    intent(inout) :: elrmn(e%ndime,e%pnode)

 
      real(rp)                   :: aux1(e%ndime)
      integer(ip)                :: idime

      aux1 = matmul(elrem,grmom)

      do idime=1,e%ndime
         elrmn(idime,:) = elrmn(idime,:) &
		          - aux*aux1(idime)*e%shape(:,e%igaus)*dvolstab
      end do

   end subroutine nsc_elmrm_ngm3

end module
