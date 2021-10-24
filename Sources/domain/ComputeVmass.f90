!> \brief Computes a Diagonal Lumped Mass matrix for the Current Mesh. <BR>
!> The Matrix is stored inverted (1/mass), so that no divisions are required when applying it
!! @param vmass  Inverted Mass matrix
subroutine ComputeVmass(a)
   use typre
   use def_parame
   use Mod_int2str
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh), target :: a

   integer(ip)           :: ielty,istat,lnodf(a%mnode),prule
   integer(ip)           :: ierr  !< my house
   integer(ip)           :: islav,imast,ibopo,iboun,idime,ielem,igaus,inodb,ipoin,ispos,ispos2,linea,pnode,pblty,inode
   logical(lg)           :: slerr
   real(rp)    :: posgc(a%ndime,a%mnode,a%nelty)       
   real(rp)    :: shapc(a%mnode,a%mnode,a%nelty)
   real(rp),target    :: deric(a%ndime,a%mnode,a%mnode,a%nelty)
   real(rp)    :: heslc(a%ntens,a%mnode,a%mnode,a%nelty)
   real(rp),target    :: weigc(a%mnode,a%nelty)
   real(rp)    :: auxshape(a%mnode)
   
   class(FiniteElement), pointer :: e => NULL()  
   real(rp)    :: dvol
   
   integer(ip) :: auxkfl_alemov = 0, initialPoint

   call a%Timer%ComputeVmass%Tic
   
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
   
   !FOR VMASS
   call a%Memor%alloc(a%npoin,a%vmass,'VMASS','exnor')
   !We force a closed integration rule for computing vmass
   !Loop over elements and compute vmass
   call a%ElementAlloc(e,a%Memor,'ForceClosedRule','vmass')
   
   elements : do ielem = 1,a%nelem
      !Load Element
      call a%ElementLoad(ielem,e)    
      call e%elmdel
      
      !Gauss Point Loop
      gauss_points : do igaus=1,e%pgaus
        e%igaus = igaus
        call e%elmder
        dvol = e%weigp(e%igaus)*e%detjm
        !a%vmass(e%lnods(1:e%pnode)) = a%vmass(e%lnods(1:e%pnode))+dvol*e%shape(1:e%pnode,e%igaus)
        auxshape(1:e%pnode) = dvol*e%shape(1:e%pnode,e%igaus)
        call a%AssemblyToArray(e,1_ip,auxshape,a%vmass)
      enddo   gauss_points
   enddo elements 
   !call a%EndAssemblyToArray(1_ip,a%vmass)
   
   call a%ElementDealloc(e,a%Memor,'ForceClosedRule','vmass')

   
   !Communicate between subdomains
   call a%ArrayCommunicator%GhostCommunicate(1_ip,a%Vmass)
   
   !Periodic Boundary conditions
   if (a%kfl_perio == 1) call a%MasterToSlave(1_ip,a%Vmass)
   
   !Store the inverse of vmass
   do ipoin=1,a%npoin   
      if (a%vmass(ipoin) /= 0.0_rp) then
         a%vmass(ipoin)=1.0_rp/a%vmass(ipoin)
         if (a%vmass(ipoin) <= 0.0_rp) then
            call a%Global2Initial(ipoin,initialPoint)
            call runend('Negative vmass value, initial point: '//adjustl(trim(int2str(initialPoint))))
         endif
      endif
   end do
   
   !undo first pass correction
   if (auxkfl_alemov == 1) a%kfl_alemov = 1
  
   call a%Timer%ComputeVmass%Toc
  
end subroutine
