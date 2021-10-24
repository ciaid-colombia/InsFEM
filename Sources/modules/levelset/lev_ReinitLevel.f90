subroutine lev_ReinitLevel(a)
!DESCRIPTION
!   This get the type of element: cut=1 and non cut=0
!-----------------------------------------------------------------------
   use typre
   use def_parame
   use MPI
   use Mod_Memor   
   use Mod_Mesh
   use Mod_Element  
   use Mod_CutMesh
   use Mod_LevelSet
   implicit none
   class(LevelSetProblem) :: a  
   
   interface
      subroutine lev_LayerOut(a, NelemWithAtLeastOneauxPointDone,ElementDone,PointDone,ElementListByLayers)
         use typre
         use Mod_LevelSet
         use Mod_Element
         implicit none
         class(LevelSetProblem) :: a
         integer(ip) :: NelemWithAtLeastOneauxPointDone
         logical :: ElementDone(*),PointDone(*)
         integer(ip) :: ElementListByLayers(*)
      end subroutine
      
      subroutine lev_reinitElmope(LevelProblem,auxlevel)
         use typre
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem), target :: LevelProblem
         real(rp) :: auxlevel(:)
      end subroutine
   
   end interface
   
   class(FiniteElement), pointer :: e => NULL()
   
   integer(ip) :: icomp,ipoin,npoin,ielem,nelem,nnode,inode,idime,ndime,jnode
   real(rp), allocatable :: ellev(:,:)  
   real(rp), allocatable :: elmat(:,:), elrhs(:)

   
   !auxiliar variables
   integer(ip)              :: auxlelemcut,elemi,elemk,elemj
   integer(ip)              :: nNodeDone,nNodeToDo,nelemList,auxPointList,nPointDone
   real(rp)                 :: distance,mindist,ratio,hipo,hipo2,auxhipo
   integer(ip)              :: jpoin,kpoin,knode,nodej,poini,nodek,jelem,nipoin,igaus,poinj
   integer(ip)              :: lpoin,hpoin,pgaus,kelem,hnode
   
   real(rp), allocatable    :: grlev(:)  
   real(rp), pointer        :: points(:,:) => NULL()    
   real(rp)                 :: iratio(4),pdsurf,grnorm,weigp(4)
   
   real(rp), pointer        :: xloc(:,:) => NULL()
   
   logical, allocatable     :: ElementDone(:),PointDone(:),PointListaux(:)
   
   integer(ip), allocatable :: PointList(:)
   
   !Assembly Variables
   real(rp), allocatable    :: rhs(:),lhs(:)   
   real(rp), allocatable    :: auxlevel(:)
   integer(ip)              :: gnode,nodeg,fnode,nodef,jdime,idir,fpoin,gpoin,poinStatus,ncutelem   

   integer(ip) :: WeAreDone, IamDone, StartElementLayer, auxAdd, NelemWithAtLeastOneauxPointDone
   logical, allocatable :: auxpointDone(:)
   integer(ip) :: ierr
  
   interface
      subroutine lev_zigzag(a,e,mnode,StartElementLayer,auxPointList,PointDone,PointListaux,PointList,auxlevel)
      !DESCRIPTION
      !   This get the type of element: cut=1 and non cut=0
      !-----------------------------------------------------------------------
         use typre
         use Mod_Mesh
         use Mod_Memor   
         use Mod_LevelSet
         use Mod_Element
         use Mod_CutMesh
         use def_parame
         implicit none
         class(LevelSetProblem) :: a   
         class(FiniteElement) :: e
         
         integer(ip), intent(in)    :: mnode  
         integer(ip)                :: nelem,ndime,npoin
         real(rp), intent(inout)    :: auxlevel(*)
         logical, intent(inout)     :: PointDone(*),PointListaux(*)
         integer(ip), intent(inout) :: PointList(*),auxPointList
         integer(ip)                :: StartElementLayer
      end subroutine
      
      subroutine lev_distance(a,e,mnode,StartElementLayer,auxPointList,PointDone,PointListaux,PointList,auxlevel)
      !DESCRIPTION
      !   This get the type of element: cut=1 and non cut=0
      !-----------------------------------------------------------------------
         use typre
         use Mod_Mesh
         use Mod_Memor   
         use Mod_LevelSet
         use Mod_Element
         use Mod_CutMesh
         use def_parame
         implicit none
         class(LevelSetProblem) :: a   
         class(FiniteElement) :: e
         
         integer(ip), intent(in)    :: mnode  
         integer(ip)                :: nelem,ndime,npoin
         real(rp), intent(inout)    :: auxlevel(*)
         logical, intent(inout)     :: PointDone(*),PointListaux(*)
         integer(ip), intent(inout) :: PointList(*),auxPointList
         integer(ip) :: StartElementLayer
      end subroutine      
   
   
   end interface     
   
   real(rp), pointer :: coord(:,:) => NULL()
   
   call a%Mesh%GetNpoin(npoin)   
   call a%Mesh%GetNelem(nelem)  
   call a%Mesh%GetNdime(ndime) 
   call a%Mesh%GetCoord(coord)
   
   call a%CutMesh%GetNCutElements(ncutelem)   
   
   !allocate
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_ReinitLevel')
   
   call a%Memor%alloc(e%mnode,1,ellev,'ellev','lev_ReinitLevel')
   call a%Memor%alloc(nelem,ElementDone,'ElementDone','lev_ReinitLevel')
   call a%Memor%alloc(npoin,PointDone,'PointDone','lev_ReinitLevel')
 
   !call a%Memor%alloc(ndime,4,xloc,'xloc','lev_ReinitLevel') 
   call a%Memor%alloc(ndime,grlev,'grlev','lev_ReinitLevel')  
   !elmat and elrhs
   call a%Memor%alloc(e%mnode,e%mnode,elmat,'elmat','lev_ReinitLevel')
   call a%Memor%alloc(e%mnode,elrhs,'elrhs','lev_ReinitLevel')
   !Assembly arrays
   call a%Memor%alloc(npoin,rhs,'rhs','lev_ReinitLevel')
   call a%Memor%alloc(npoin,lhs,'lhs','lev_ReinitLevel')
   call a%Memor%alloc(npoin,PointList,'PointList','lev_ReinitLevel')
   !auxiliar list
   call a%Memor%alloc(npoin,PointListaux,'PointListaux','lev_ReinitLevel')   
   call a%Memor%alloc(npoin,auxlevel,'auxlevel','lev_ReinitLevel')
   
   auxlevel = a%level(:,3)

   
   !initialization
   PointDone= .false.
   ElementDone = .false.   
   PointListaux = .false.
   auxPointList = 0
  
   elemi=1   
   !Interface Elements
   do elemi = 1,ncutelem   
      !List of Cut elements
      ielem = a%ElementListByLayers(elemi)       
         
      call a%Mesh%ElementLoad(ielem,e)  
      call e%gather(1,ellev(:,1),a%level(:,3))   
      
      call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
      call a%CutMesh%GetSurfaceIntersection(ielem,e,pdsurf)
      call a%CutMesh%GetNinters(ielem,nipoin)
      
        
      !Compute shape functions at the intersection points
      do kpoin = 1, nipoin         
         weigp(kpoin) = 1.0_rp/nipoin
      end do      
         
      !The rutine gives the needed shape functions associated 
      call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)
         
      !Levelset Gradient 
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg      
      
      call e%gradient(1,ellev(:,1),grlev)      
      call vecnor(grlev,ndime,grnorm,2)
      if(grnorm==0.0_rp) grnorm=1.0_rp    

      do inode= 1,e%pnode
         !Initialization 
         elmat = 0.0_rp
         elrhs = 0.0_rp      
      
         ipoin=e%lnods(inode)          
         
         !distance function
         call a%CutMesh%GetPointType(ipoin,poinStatus)
         if(poinStatus == 1)then
            
            if(PointListaux(ipoin) .eqv. .false.)then
               PointListaux(ipoin) = .true.
               !Auxiliar List
               auxPointList = auxPointList + 1
               PointList(auxPointList) = ipoin            
            end if
            
            do igaus = 1,e%pgaus
               e%igaus = igaus
                 
               elmat(inode,inode) = elmat(inode,inode) + weigp(e%igaus)*pdsurf
               elrhs(inode) = elrhs(inode) + weigp(e%igaus)*pdsurf*a%level(ipoin,3)/grnorm
            enddo 
            
            !Assembly
            lhs(ipoin) = lhs(ipoin) + elmat(inode,inode) 
            rhs(ipoin) = rhs(ipoin) + elrhs(inode)  
            
         end if
      end do      
   end do
    
   !Solve the system 
   poini=1
   do poini = 1,auxPointList
      ipoin = PointList(poini)      
      PointDone(ipoin)=.true.         
      auxlevel(ipoin) = rhs(ipoin)/lhs(ipoin)   
   end do
  
   nPointDone = auxPointList  
     
   !Mass matrix approach  
   elemk=1  
   do elemk = 1,ncutelem   
      !List of Cut elements
      kelem = a%ElementListByLayers(elemk)      
         
      call a%Mesh%ElementLoad(kelem,e)  
      call e%gather(1,ellev(:,1),a%level(:,3))      
      
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg  
       
      call a%CutMesh%GetLocalIntersectionPoints(kelem,xloc)
      call a%CutMesh%GetSurfaceIntersection(kelem,e,pdsurf)
      call a%CutMesh%GetNinters(kelem,nipoin)
      
      !Compute shape functions at the intersection points
      do kpoin = 1, nipoin         
         weigp(kpoin) = 1.0_rp/nipoin
      end do      
         
      !The rutine give the needed shape functions associated 
      call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)  
      
      do inode= 1,e%pnode
      
         !Initialization 
         elmat = 0.0_rp
         elrhs = 0.0_rp
         
         call a%CutMesh%GetPointType(e%lnods(inode),poinStatus)         
         if(poinStatus == -1)then
         
            ipoin=e%lnods(inode)          
                  
            if(PointDone(ipoin) .eqv. .false.)then            
             
               if(PointListaux(ipoin) .eqv. .false.)then
                  PointListaux(ipoin) = .true.
                  !Auxiliar List
                  auxPointList = auxPointList + 1
                  PointList(auxPointList) = ipoin            
               end if            
               
               do igaus = 1,e%pgaus
                  e%igaus = igaus                    
 
                  elmat(inode,inode) = elmat(inode,inode) &
                    + e%shape(inode,e%igaus)*e%shape(inode,e%igaus)*weigp(e%igaus)*pdsurf

               enddo
               
               do jnode = 1,e%pnode
                  jpoin = e%lnods(jnode)                  

                  if(PointDone(jpoin) .eqv. .true.)then
               
                     do igaus = 1,e%pgaus
                        e%igaus = igaus
                        elrhs(inode) = elrhs(inode) &
                           - e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*weigp(e%igaus)*pdsurf*auxlevel(jpoin)
                     end do
                  endif
               enddo   

               !Assembly
               lhs(ipoin) = lhs(ipoin) + elmat(inode,inode)                
               rhs(ipoin) = rhs(ipoin) + elrhs(inode)
            endif               
         end if
      end do   
   end do    
   
   !Solve the system
   poinj=1
   do poinj = nPointDone+1,auxPointList
      jpoin = PointList(poinj)      
      PointDone(jpoin)=.true.         
      auxlevel(jpoin) = rhs(jpoin)/lhs(jpoin)
   end do
   
   nPointDone = auxPointList   
   
   
   if (a%kfl_ReinitLevel == 1 .or. a%kfl_ReinitLevel == 2) then
      !For Zig-Zag and distance
      WeareDone = 0
      IamDone = 0
      
      if (ncutelem > 0) then
         StartElementLayer = ncutelem
      else
         StartElementLayer = -1
      endif
      
      call a%Memor%alloc(npoin,auxPointDone,'auxPointDone','lev_ReinitLevel')
      
      do while (WeAreDone == 0) 
         
         if (StartElementLayer > -1) then
            IAmDone = 1
         endif
         
         if (StartElementLayer > -1) then
         
            !remaining layers of elements   
            !---------------------------------------------------------------------------------------------
            if(a%kfl_ReinitLevel==2)then
               call lev_distance(a,e,e%mnode,StartElementLayer,auxPointList,PointDone,PointListaux,PointList,auxlevel)
            
            elseif(a%kfl_ReinitLevel==1)then
            
               call lev_zigzag(a,e,e%mnode,StartElementLayer,auxPointList,PointDone,PointListaux,PointList,auxlevel)
            endif
         endif
         
         call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,PointDone)
         call a%Mesh%ArrayCommunicator%GhostCommunicate(1_ip,auxlevel)
         
         auxPointDone = PointDone
         
         if (StartElementLayer == -1) then
            kelem = 0
            call a%Mesh%GetNelem(nelem)
            do ielem = 1,nelem
               ElementDone(ielem) = .false.
               call a%Mesh%ElementLoad(ielem,e)  

               auxAdd = 0
               do inode = 1,e%pnode
                  ipoin = e%lnods(inode)
                  if (PointDone(ipoin) .eqv. .true.) then
                     auxAdd = auxAdd+1
                  endif
               enddo
               ! in linear elements we add elements when only one node is not calculated
               if (auxAdd>=ndime) then
                  kelem = kelem + 1
                  a%ElementListByLayers(kelem) = ielem
                  ElementDone(ielem) = .true.
                  auxPointDone(e%lnods(1:e%pnode)) = .true.
               endif
               
            enddo
            NelemWithAtLeastOneauxPointDone = kelem
            
            
            call lev_LayerOut(a, NelemWithAtLeastOneauxPointDone,ElementDone,auxPointDone,a%ElementListByLayers)
            
            
            if (all(ElementDone)) then
               StartElementLayer = 0
            endif
         endif
         
         call  MPI_AllREDUCE(IAmDone, WeareDone, 1, MPI_INTEGER4, MPI_MIN,a%MPIcomm, ierr )
         
      enddo
      
      call a%Memor%dealloc(npoin,auxPointDone,'auxPointDone','lev_ReinitLevel')
   
   elseif (a%kfl_ReinitLevel == 3) then
      
      !Solve a poisson problem
      call lev_reinitElmope(a,auxlevel)
   
   endif
   
   
   
   
   
   
   
   
   
   
   
   
   !-------------------------------------------------------------------------------------------------------------
   !End of the reinitialization
   a%level(:,3) = auxlevel  
   
   !allocate
   call a%Memor%dealloc(e%mnode,1,ellev,'ellev','lev_ReinitLevel')
   call a%Memor%dealloc(nelem,ElementDone,'ElementDone','lev_ReinitLevel')
   call a%Memor%dealloc(npoin,PointDone,'PointDone','lev_ReinitLevel') 

   !call a%Memor%dealloc(ndime,4,xloc,'xloc','lev_ReinitLevel')   
   call a%Memor%dealloc(ndime,grlev,'grlev','lev_ReinitLevel')
   !elmat and elrhs
   call a%Memor%dealloc(e%mnode,e%mnode,elmat,'elmat','lev_ReinitLevel')
   call a%Memor%dealloc(e%mnode,elrhs,'elrhs','lev_ReinitLevel')   
   !Assembly arrays
   call a%Memor%dealloc(npoin,rhs,'rhs','lev_ReinitLevel')
   call a%Memor%dealloc(npoin,lhs,'lhs','lev_ReinitLevel')
   call a%Memor%dealloc(npoin,PointList,'PointList','lev_ReinitLevel')   
   !auxiliar list
   call a%Memor%dealloc(npoin,PointListaux,'PointListaux','lev_ReinitLevel')   
   call a%Memor%dealloc(npoin,auxlevel,'auxlevel','lev_ReinitLevel')

   !dealloc element
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','lev_ReinitLevel')     
   
end subroutine
