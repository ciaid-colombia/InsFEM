subroutine rom_basisFilter(a)
   use typre
   use def_parame
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   real(rp)    :: basis_energy
   integer(ip) :: ibasis

   call a%Memor%alloc(a%nconv,a%bEnergy,'bEnergy','rom_Memall')
   a%bEnergyT = 0.0_rp
   a%bEnergy = 0.0_rp
   basis_energy = 0.0_rp

   do ibasis=2,a%nconv
      a%bEnergyT = a%bEnergyT + a%sigma(ibasis)
      a%bEnergy(ibasis) = a%bEnergy(ibasis-1) + a%sigma(ibasis)
   end do

   if (a%kfl_basis_number) then
       a%ndofr = a%basis_number
       basis_energy = a%bEnergy(a%basis_number)/a%bEnergyT
       a%basis_energy = a%bEnergy(a%basis_number)/a%bEnergyT
   else
      ibasis = 1
      do while (basis_energy < a%basis_energy .and. ibasis < a%nconv)
         ibasis = ibasis+1
         basis_energy = a%bEnergy(ibasis)/a%bEnergyT
      end do
      a%ndofr  = ibasis
      a%basis_number = ibasis
   end if
   
   call a%EigenSystem%SetBasisSize(a%ndofr)
   if (a%MPIrank == a%MPIroot) then
      write(a%lun_outro,*) "Retained energy:         ", (basis_energy)
      write(a%lun_outro,*) "Basis size:              ", (a%ndofr)
   end if

end subroutine
