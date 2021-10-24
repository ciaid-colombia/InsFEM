module Mod_NSCompressibleSubroutines
contains

   subroutine nsc_ComputeEnergy(ndime,accvh,tempe,veloc,energ)
      use typre
      use def_parame
      implicit none
      
      integer(ip), intent(in)    :: ndime
      real(rp),    intent(in)    :: accvh,tempe
      real(rp),    intent(in)    :: veloc(ndime)
      real(rp),    intent(out)   :: energ

      real(rp)                   :: tmpvel(ndime),tmp

      tmpvel = veloc
      tmp = dot_product(tmpvel,veloc) / 2.0_rp
      energ = accvh*tempe + tmp
      
   end subroutine   

   subroutine nsc_ComputeAdvectionVelocity(e,pcden,pcmom,pcvno)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: pcden,pcmom(e%ndime)
      real(rp),    intent(out)   :: pcvno

      real(rp)                   :: tmp(e%ndime)
     
      tmp = pcmom/pcden 
      call vecnor(tmp,e%ndime,pcvno,2)
      
   end subroutine     

   subroutine nsc_ComputeTemperature(a,e,pcden,pcmom,pcene,pctem)
      use typre
      use Mod_NSCompressible
      use Mod_Element
      implicit none
      class(NSCompressibleProblem) :: a
      class(FiniteElement)        :: e
      real(rp),    intent(in)    :: pcden,pcmom(e%ndime),pcene
      real(rp),    intent(out)   :: pctem

      real(rp)                   :: acvis,actco,accph,accvh
      real(rp)                   :: invcvh,invpcd,aux_k

      !Physical Parameters
      call a%GetPhysicalParameters(acvis,actco,accph,accvh)

      invcvh = 1.0_rp / accvh 
      invpcd = 1.0_rp / pcden
      aux_k = dot_product(pcmom,pcmom) / 2.0_rp
     
      pctem = invpcd * invcvh * (pcene - (aux_k*invpcd))
      
   end subroutine     

   subroutine nsc_ComputeSoundSpeed(accph,accvh,pctem,pcspd)
      use typre
      implicit none
      real(rp),    intent(in)    :: accph,accvh,pctem
      real(rp),    intent(out)   :: pcspd

      real(rp)                   :: cgamma

      cgamma = accph/accvh
      pcspd = sqrt(accph * (cgamma-1.0_rp) * pctem)
      
   end subroutine     

end module
