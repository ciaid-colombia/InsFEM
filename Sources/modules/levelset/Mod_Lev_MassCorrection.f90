module Mod_Lev_MassCorrection
   use MPI   
   use Mod_LevelSet
   use Mod_lev_BaseElmope
   use Mod_PointerSetter
   implicit none
   
   real(rp) :: MeanLevsetGradient, TotalFreeSurfaceArea, TotalVolume
   
   real(rp), allocatable :: PointLevsetGradient(:), PointPdSurf(:)
   
   
   
   type, extends(PointerSetter) :: SPMassCorrection
contains
      procedure :: SpecificSet => SpecificSetMassCorrection
   end type
   
   type(SPMassCorrection) :: SetPointersMassCorrection
   
   
   
   
   
   
   
   
contains

   subroutine SpecificSetMassCorrection(d)
      implicit none
      class(SPMassCorrection) :: d
      
      !a%kfl_MassCorrection = 1
      if (a%kfl_MassCorrection == 1) then
         call ConcatenateProcedures(ProcHook%Initializations,InitializationsMassCorrection)
         call ConcatenateProcedures(ProcHook%PostGaussElmats,PostGaussMassCorrection)
         call ConcatenateProcedures(ProcHook%Finalizations,FinalizationsMassCorrection)
      endif
   
   end subroutine



   subroutine InitializationsMassCorrection
      implicit none
      
      integer(ip) :: npoin
      
      TotalVolume = 0.0_rp
      MeanLevsetGradient = 0.0_rp
      TotalFreeSurfaceArea = 0.0_rp
      
      call a%Mesh%Getnpoin(npoin)
      
      call a%Memor%alloc(npoin,PointLevsetGradient,'PointLevsetGradient','MassCorrection')
      call a%Memor%alloc(npoin,PointPdSurf,'PointPdSurf','MassCorrection')
      
   
   end subroutine
   
   
   subroutine PostGaussMassCorrection
      implicit none
      
      real(rp) :: grlev(3), UnitOrthogonalVector(3)
      
      integer(ip) :: elemStatus, ngaus_total, ngauss_plus, ngauss_minus
      real(rp)    :: weigp(40)
      real(rp)    :: pdsurf, NormalLevsetGradient
      
      integer(ip) :: inode,ipoin,npoinLocal
      real(rp)    :: weightfactor
      
      
      
      !Count the number of nodes which are mine
      call a%Mesh%GetNpoinLocal(npoinLocal)
      weightfactor = 0.0_rp
      do inode = 1,e%pnode
         ipoin = e%lnods(inode)
         if (ipoin <= npoinLocal) then
            weightfactor = weightfactor + 1.0_rp
         endif
      enddo
      weightfactor = weightfactor/real(e%pnode,rp)
         
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
         
      if (elemStatus == 1) then
         !TotalVolume
         TotalVolume = TotalVolume + weightfactor*sum(e%weigp(1:e%pgaus))*e%detjm
      
      elseif(elemStatus ==0)then    
         
         ngaus_total=0
   
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         call a%CutMesh%GetWeigpCut(ielem,e%ndime,weigp)
         
         !TotalVolume
         TotalVolume = TotalVolume + weightfactor*sum(weigp(ngauss_minus+1:ngauss_minus+ngauss_plus))*e%detjm
        
         call a%CutMesh%GetExternalNormal(ielem,e,UnitOrthogonalVector)
         call a%CutMesh%GetSurfaceIntersection(ielem,e,pdsurf)
         
         call e%gradient(1_ip,ellev(:,1),grlev)
          
         NormalLevsetGradient = dot_product(grlev(1:e%ndime),UnitOrthogonalVector(1:e%ndime))
         
         !PointLevsetGradient
         PointLevsetGradient(e%lnods(1:e%pnode)) = PointLevsetGradient(e%lnods(1:e%pnode)) + NormalLevsetGradient*pdsurf
         PointPdSurf(e%lnods(1:e%pnode)) = PointPdSurf(e%lnods(1:e%pnode)) + pdsurf
         
         !MeanLevsetGradient 
         MeanLevsetGradient = MeanLevsetGradient + weightfactor*NormalLevsetGradient*pdsurf
         TotalFreeSurfaceArea = TotalFreeSurfaceArea + weightfactor*pdsurf
     
      end if
   end subroutine
   
   
   subroutine FinalizationsMassCorrection
      implicit none
      
      real(rp) :: CorrectionVolume, CorrectionLevset, PointCorrection,CorrectionDistance
      integer(ip) :: npoin,ipoin
      real(rp) :: auxTotalSurface, auxTotalVolume
      integer(ip) :: ierr

      !Communications, volume and surface
      auxTotalVolume = TotalVolume
      call  MPI_AllREDUCE( auxTotalVolume, TotalVolume, 1, MPI_REAL8, MPI_SUM,a%MPIcomm, ierr )
      
      auxTotalSurface = TotalFreeSurfaceArea
      call  MPI_AllREDUCE( auxTotalSurface, TotalFreeSurfaceArea, 1, MPI_REAL8, MPI_SUM,a%MPIcomm, ierr )
      
      !Ghost communcations
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,PointLevsetGradient)
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,PointPdSurf)
      
      
      !The first time we only count the volume
      if (a%ipass == 0) then
         a%ReferenceVolume = TotalVolume
         a%ipass = 1
         
      else
         CorrectionVolume = a%ReferenceVolume - TotalVolume
         
         CorrectionDistance = CorrectionVolume/TotalFreeSurfaceArea
         
         !call a%FilePostpr%postpr(a%level(:,1),'levelprecorrection',a%istep,a%ctime,a%Mesh)
         
         call a%Mesh%GetNPoin(npoin)
         
         do ipoin = 1,npoin
            if (PointPdSurf(ipoin) /= 0.0_rp) then
               PointCorrection = CorrectionDistance*PointLevsetGradient(ipoin)/PointPdSurf(ipoin)
               a%level(ipoin,1) = a%level(ipoin,1) - PointCorrection
            endif
         enddo
         
         write(*,*) 'CorrectionDistance: ', CorrectionDistance
         !call a%FilePostpr%postpr(a%level(:,1),'levelpostcorrection',a%istep,a%ctime,a%Mesh)
         
         
      endif
      
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%dealloc(npoin,PointLevsetGradient,'PointLevsetGradient','MassCorrection')
      call a%Memor%dealloc(npoin,PointPdSurf,'PointPdSurf','MassCorrection')
         
   
   end subroutine




end module

