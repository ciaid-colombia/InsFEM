!> \brief Computes a Diagonal Lumped Boundary Mass matrix for the Current Mesh. <BR>
!> The Matrix is stored inverted (1/mass), so that no divisions are required when applying it
!! @param bvmass Boundary Inverted Mass matrix
subroutine nsc_CalculateBoundaryVmass(a)
   use typre
   use def_parame
   use Mod_Memor
   use Mod_Mesh
   use Mod_Element
   use Mod_NSCompressible

   implicit none
   class(NSCompressibleProblem) :: a
   

   integer(ip)             :: iboun,nboun,npoin
   integer(ip)             :: ipoin,ndime,igaub
   integer(ip)             :: calculate_statistics

   class(FiniteElement), pointer :: e => NULL()
   real(rp)    :: auxshape(a%Mesh%mnode)
   real(rp)    :: dsurf

   !Initializations
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNboun(nboun)   

   !FOR BVMASS
   call a%Memor%alloc(npoin,a%bvmass,'BVMASS','nsc_BoundaryVmass')
   a%bvmass=0.0_rp

   !We force a closed integration rule for computing vmass
   !Loop over elements and compute vmass
   call a%Mesh%ElementAlloc(e,a%Memor,'ForceClosedRule','nsc_BoundaryVmass')

   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)
      
      call e%elmdel
      
      !Gauss Point Loop
      gauss_points : do igaub=1,e%pgaub
        e%igaub = igaub
        !Calculate exterior Normal
        call e%bounor
        
        auxshape = 0.0_rp
        dsurf=e%weigb(e%igaub)*e%eucta
        auxshape(e%lboel(1:e%pnodb)) = dsurf*e%shapb(1:e%pnodb,e%igaub)
        call a%Mesh%AssemblyToArray(e,1_ip,auxshape,a%bvmass)
      enddo   gauss_points

   end do boundaries 

   call a%Mesh%ElementDealloc(e,a%Memor,'ForceClosedRule','nsc_BoundaryVmass')

   !Communicate between subdomains
   call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,a%bvmass)
   
   !Store the inverse of vmass
   do ipoin=1,npoin   
      if (a%bvmass(ipoin) /= 0.0_rp) then
         a%bvmass(ipoin)=1.0_rp/a%bvmass(ipoin)
         if (a%bvmass(ipoin) < 0.0_rp) call runend('Negative bvmass value')
      endif
   end do

end subroutine nsc_CalculateBoundaryVmass
