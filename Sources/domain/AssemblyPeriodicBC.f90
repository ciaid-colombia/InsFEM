subroutine AssemblyPeriodicBC(a,ndofn,LinearSystem,Memor)
   use typre
   use Mod_ParallelSystemInterface
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   class(ParallelSystemInterface) :: LinearSystem
   type(MemoryMan) :: Memor
   real(rp) :: elmat(ndofn,2,ndofn,2),elrhs(ndofn,2)
   real(rp), allocatable :: MatDiag(:,:)
   type(FiniteElement) :: e
   integer(ip) :: idofn,imaster,islave,ipoin
   integer(ip) :: iwa(ndofn*100)
   real(rp)    :: rwa(ndofn*100)
   real(rp)    :: diag(ndofn)
   
   elmat = 0.0_rp
   elrhs = 0.0_rp

   !Manual allocation of the element e
   call Memor%palloc(2_ip,e%lnods,'e%lnods','AssemblyPeriodicBC')
   e%pnode = 2
   e%mnode = 2
   
   do ipoin = 1,a%nslave
      islave = a%SlaveList(ipoin)
      imaster = a%MasterSlave(islave)

      e%lnods(1) = islave
      e%lnods(2) = imaster
      
      do idofn = 1,ndofn
         elmat(idofn,1,idofn,1) = 1
      enddo
      
      e%lnods(1) = islave
      e%lnods(2) = imaster
      
      !AssemblyPeriodicBC
      call LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo   

   call Memor%pdealloc(2,e%lnods,'e%lnods','AssemblyPeriodicBC')

end subroutine

subroutine AssemblyPeriodicBCToZero(a,ndofn,LinearSystem,Memor)
   use typre
   use Mod_ParallelSystemInterface
   use Mod_Memor
   use Mod_Element
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip) :: ndofn
   class(ParallelSystemInterface) :: LinearSystem
   type(MemoryMan) :: Memor
   
   real(rp) :: elmat(ndofn,2,ndofn,2),elrhs(ndofn,2)
   real(rp), allocatable :: MatDiag(:,:)
   
   type(FiniteElement) :: e
   integer(ip) :: idofn,imaster,islave,ipoin
   
   integer(ip) :: iwa(ndofn*100)
   real(rp)    :: rwa(ndofn*100)
   
   real(rp)    :: diag(ndofn)
   
   elmat = 0.0_rp
   elrhs = 0.0_rp

   !Manual allocation of the element e
   call Memor%palloc(2_ip,e%lnods,'e%lnods','AssemblyPeriodicBC')
   e%pnode = 2
   e%mnode = 2
   
   do ipoin = 1,a%nslave
      islave = a%SlaveList(ipoin)
      imaster = a%MasterSlave(islave)

      e%lnods(1) = islave
      e%lnods(2) = imaster
      
      do idofn = 1,ndofn
         elmat(idofn,1,idofn,1) = 0
      enddo
      
      e%lnods(1) = islave
      e%lnods(2) = imaster
      
      !AssemblyPeriodicBC
      call LinearSystem%Assembly(e,elmat,elrhs)
      
   enddo   

   call Memor%pdealloc(2,e%lnods,'e%lnods','AssemblyPeriodicBC')

end subroutine
