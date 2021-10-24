module Mod_plcd_ExacsoError
   use typre
   use Mod_plcd_BaseElmope
   use Mod_plcdExacso
   use MPI
   use Mod_plcd_StrainGenerator
   implicit none
   private
   public SetPointersExacsoError
   
   type, extends(PointerSetter) :: SPExacsoError
contains
      procedure :: SpecificSet => SpecificSetExacsoError
   end type
   
   type(SPExacsoError) :: SetPointersExacsoError
   
   integer(ip), allocatable :: kfl_IsSet
   
   type(plcdExacso) :: Exacso
   real(rp), allocatable :: L2error(:), NodalError(:,:), EnergyNormError(:),elexactdisp(:,:), elexactpre(:), PErrorL2(:),elpre(:,:)
   
   real(rp) :: gpcod(3)
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetExacsoError(d)
      class(SPExacsoError) :: d
         
      !Element size required by EMDS
      if (a%kfl_exacs > 0) then
         call ConcatenateProcedures(ProcHook%Initializations,Initializations)
         call ConcatenateProcedures(ProcHook%Finalizations,Finalizations)
         call ConcatenateProcedures(ProcHook%PreGauss,Gathers)
         call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmats)
      endif
   end subroutine   
   
   
   subroutine Initializations
      implicit none
      integer(ip) :: ndime,npoin
      
      call a%Mesh%GetNelem(nelem)
      call a%Memor%alloc(nelem,L2error,'error','plcd_exacso')
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Memor%alloc(ndime,npoin,NodalError,'error','plcd_exacso')
      
      call a%Memor%alloc(nelem,EnergyNormError,'EnergyNormError','plcd_exacso')
      
      call a%Memor%alloc(e%ndime,e%mnode,1,eldisp,'eldisp','plcd_exacso')
      
      if (a%UseUPFormulation) then
         call a%Memor%alloc(nelem,PErrorL2,'PErrorL2','plcd_exacso')
         call a%Memor%alloc(1,e%mnode,elpre,'elpre','plcd_exacso')
      endif
   
   end subroutine
   
   subroutine Gathers
      implicit none
      integer(ip) :: npoinLocal,nghost
      
      call e%gather(e%ndime,eldisp,a%displacement)
      call e%elmdcg
      call e%elmlen
      
      if (a%UseUPFormulation) then
         call e%gather(1,elpre,a%pressure)
      endif
   end subroutine
   
   subroutine InGaussElmats
      implicit none
      real(rp) :: exvel(3), exveg(e%ndime,e%ndime)
      real(rp) :: errnorm2, gpdisp(3)
      
      class(PLCDMaterial), pointer :: Material
      integer(ip) :: vsize,ndime,inode,ipoin
      real(rp) ::  ErrorStress(6), ErrorStrain(6),EnergyNorm2, gDispError(e%ndime,e%ndime),gDisp(e%ndime,e%ndime), expre, expregr(3), gppre(1),G, Cnorm
      
      real(rp), pointer :: coord(:)
      real(rp) :: grapre(3), tau
   
   
      call e%interpg(e%ndime,e%elcod,gpcod)
      call Exacso%ComputeSolution(e%ndime,gpcod,a)
      call Exacso%GetDisplacement(e%ndime,exvel,exveg)
      
      call ElementMatData%GetMaterialPointer(Material)
      call Material%CT%GetCNorm(Cnorm)
      
      !-------------------------------------------------------------------------
      !For L2error
      
      call e%interpg(e%ndime,eldisp,gpdisp)
      !L2 error
      call e%elmlen
      errnorm2 = dot_product(gpdisp(1:e%ndime)-exvel(1:e%ndime),gpdisp(1:e%ndime)-exvel(1:e%ndime))*Cnorm/(e%hleng(1)**2)
      L2error(ielem) = L2error(ielem) + e%detjm*e%weigp(e%igaus)*errnorm2
      
!       write(*,*) L2error(ielem), e%detjm, e%weigp(e%igaus), errnorm2
      
      if (a%UseUPFormulation) then
         call ElementMatData%GetSecantCNorm(e%igaus,G)
         call Exacso%GetPressure(e%ndime,expre,expregr)
         call e%interpg(1,elpre,gppre)
         errnorm2 = (gppre(1)-expre)**2/G
         PErrorL2(ielem) = PErrorL2(ielem) + e%detjm*e%weigp(e%igaus)*errnorm2
      endif
      !-------------------------------------------------------------------------
      
      
      !-------------------------------------------------------------------------
      !For energy error
      call GetVoigtsize(e%ndime,vsize)
      !EnergyNormError
      !Equation 4.19 from Badia Baiges 2013
      EnergyNorm2 = 0.0_rp
      !ComputeDisplacementGradients
      call e%gradient(e%ndime,eldisp(:,:,1),gDisp)
      gDispError(1:e%ndime,1:e%ndime) = gDisp(1:e%ndime,1:e%ndime) - exveg(1:e%ndime,1:e%ndime)      
   
      call ElementMatData%GetMaterialPointer(Material)
      call Material%CT%GetStrain(gDispError,ErrorStrain)            
      call ElementMatData%GetConstitutiveTensorPointer(e%igaus,C)
      
      call Material%CT%GetCNorm(Cnorm)
      
      ErrorStress(1:vsize) =  matmul(C,ErrorStrain(1:vsize))
      
      !call Material%CT%ComputeStressNorm(Stress,StressNorm)
      
      EnergyNorm2  =EnergyNorm2  + dot_product(ErrorStress(1:vsize),ErrorStrain(1:vsize))
      
      if (a%useUpFormulation) then
         call ElementMatData%GetSecantCNorm(e%igaus,G)
         !energynorm
         EnergyNorm2 = EnergyNorm2 + ((gppre(1)-expre)**2)/(2*G)
         
         call e%gradient(1,elpre,grapre)
         call e%elmlen
         
         tau = a%staco(1)*e%hleng(1)**2/(2*G)
         
         EnergyNorm2  =EnergyNorm2  + tau*dot_product(grapre(1:e%ndime)-expregr(1:e%ndime),grapre(1:e%ndime)-expregr(1:e%ndime))
      endif
      
      EnergyNormError(ielem) = EnergyNormError(ielem) + e%detjm*e%weigp(e%igaus)*EnergyNorm2
      
      !NodalError
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         
         call a%Mesh%GetPointCoord(ipoin,coord)
         
         call Exacso%ComputeSolution(e%ndime,coord,a)
         call Exacso%GetDisplacement(e%ndime,exvel,exveg)
         
         NodalError(:,ipoin) = eldisp(:,inode,1)- exvel(1:e%ndime)
      enddo
      
   end subroutine
   
   subroutine Finalizations
      implicit none
      
      real(rp) :: TotalL2error,TotalL2Error0, TotalEnergyError, TotalEnergyError0, TotalPErrorL2, TotalPErrorL20
      integer(ip) :: ndime,npoin,ierr,npoinLocal
      real(rp) :: weightfactor
      
      
      call a%Mesh%GetNpoinLocal(npoinLocal)

      TotalL2error0 = 0.0_rp
      TotalEnergyError0 = 0.0_rp
      TotalPErrorL20 = 0.0_rp
      do ielem = 1,nelem
         call a%Mesh%ElementLoad(ielem,e)
         weightfactor = 1.0-real(count(e%lnods(1:e%pnode)>npoinLocal))/real(e%pnode)
      
      
         TotalL2error0 = TotalL2error0 + L2error(ielem)*weightfactor
!          write(*,*) L2error(ielem)*weightfactor
         TotalEnergyError0 = TotalEnergyError0 + EnergyNormError(ielem)*weightfactor
         if (a%UseUPFormulation) TotalPErrorL20 = TotalPErrorL20 + PErrorL2(ielem)*weightfactor
      enddo
      call  MPI_REDUCE(TotalL2Error0, TotalL2Error, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE(TotalEnergyError0, TotalEnergyError, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      if (a%UseUPFormulation)  call  MPI_REDUCE(TotalPErrorL20, TotalPErrorL2, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      TotalL2Error = sqrt(TotalL2Error)
      TotalEnergyError = sqrt(TotalEnergyError)
      TotalPErrorL2 = sqrt(TotalPErrorL2)
      
      
      if (a%MPIrank == a%MPIroot) write(*,*) 'Total L2 error: ', TotalL2error
      if (a%MPIrank == a%MPIroot) write(*,*) 'Total EnergyNorm error: ', TotalEnergyError
      
      call a%FilePostpr%postpr(NodalError,'NodalError',a%istep,a%ctime,a%Mesh)
      
      call a%FilePostpr%postgp(L2error,'L2error',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postgp(EnergyNormError,'EnergyNormError',a%istep,a%ctime,a%Mesh)
      
      if (a%UseUPFormulation) then
         if (a%MPIrank == a%MPIroot) write(*,*) 'Total P L2 error: ', TotalPerrorL2
         call a%FilePostpr%postgp(PErrorL2,'PErrorL2',a%istep,a%ctime,a%Mesh)
         
         !Deallocates
         call a%Memor%dealloc(nelem,PErrorL2,'PErrorL2','plcd_exacso')
         call a%Memor%dealloc(1,e%mnode,elpre,'elpre','plcd_exacso')
      endif
      
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      call a%Memor%dealloc(ndime,npoin,NodalError,'error','plcd_exacso')
      
      call a%Memor%dealloc(nelem,L2error,'error','plcd_exacso')
      call a%Memor%dealloc(nelem,EnergyNormError,'EnergyNormError','plcd_exacso')
      call a%Memor%dealloc(e%ndime,e%mnode,1,eldisp,'eldisp','plcd_exacso')
   end subroutine
   
   
   
end module
