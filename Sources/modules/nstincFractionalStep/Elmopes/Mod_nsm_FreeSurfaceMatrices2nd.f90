module Mod_nsm_FreeSurfaceMatrices2nd
   use Mod_nsm_BaseElmope
   use Mod_CutMesh
   use typre
   implicit none
   private
   public :: SetPointersFreeSurfaceMatrices2nd
   
   !External Procedures
   procedure() :: NULLSUB   
   

   !Coupling with levelset
   real(rp)                 :: gplev(1)
   integer(ip)              :: elemStatus(1),ngauss_minus,ngauss_plus,ngaus_total
   real(rp), allocatable    :: weigp(:)
   real(rp), pointer        :: xloc(:,:) => NULL()
   integer(ip), allocatable :: kfl_IsSet 
   !To define the pressure boundary condition over free surface
   !Nodal coordinates
   real(rp), pointer       :: coord(:) => NULL()   

   
contains

   !For setting the pointers
   subroutine SetPointersFreeSurfaceMatrices2nd(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            !Pointers are already set
            kfl_IsSet = 1
         
            if (a%kfl_colev == 1) then 
               !Level Set Coupling
               call ConcatenateProcedures(ProcHook_Initializations,AllocLev)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocLev)            

               if(a%kfl_fsurf==1)then
                  call ConcatenateProcedures(ProcHook_PreDirichlet,FreeSurfMats)
                  
                  if(a%kfl_fsurfLapla==1)then
                     call PrependProcedure(ProcHook_PreDirichlet,ElmatsToLapla)                     
                  end if
                  
                  !Always after FreeSurfMatsTF
                  call ConcatenateProcedures(ProcHook_PreDirichlet,PressureZeroFurf_2nd)                  
               end if      
               
            endif
            
            
         endif
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   subroutine AllocLev
      implicit none
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2
      
      call a%Memor%alloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')

      
   end subroutine   
   
   
   subroutine DeallocLev
      implicit none      
      integer(ip) :: auxdim,auxdim1
      
      auxdim=(e%ndime+1)
      auxdim1=(e%ndime-1)*(e%ndime-1)+2    
      
      call a%Memor%dealloc(auxdim*auxdim1,weigp,'weigp','nsm_elmope')
     
      
   end subroutine
   

   
   subroutine FreeSurfMats
      implicit none 
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==-1)then
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%CutMesh%GetPointType(ipoin,poinStatus)            
            if(poinStatus==-1)then
            
               elmat(:,inode,:,1:e%pnode) = 0.0_rp
               elrhs(:,inode)=0.0_rp

            
            end if         
         end do        
      end if
   end subroutine 
   
   
   subroutine ElmatsToLapla
      implicit none
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin,idime 
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      
      if(elemStatus==-1)then
         
         elmat(:,:,:,:) = 0.0_rp
         elrhs(:,:)=0.0_rp     
         wrmat1=0.0_rp                  
         
         ! Viscosity terms : we only consider mu*(grad v, grad u)         
         call elmvis(e,dvolt0,acvis,wrmat1)     
         

         elmat(1,1:e%pnode,1,1:e%pnode) = elmat(1,1:e%pnode,1,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)   
        
         do inode = 1,e%pnode
            !Rhs
            elrhs(1,inode) = elrhs(1,inode) + acvis*dot_product(e%cartd(:,inode),acden*a%grnor*a%gravi(1:e%ndime))*dvolt0 
         enddo      
           
         
      end if      
      
   end subroutine
   
   subroutine PressureZeroFurf_2nd
      implicit none
      integer(ip) :: elemStatus,poinStatus,JpoinStatus,ndone,ntodo
      integer(ip) :: nipoin,ipoin,jpoin,jnode
      real(rp)    :: pdsurf,toleralce
      !real(rp)    :: boundaryVector(3),hipo,tmp
      !real(rp)    :: UnitBoundaryVector(3),,BoundaryCenter(3),checkvector(3)
      !3d case
      !real(rp)    :: boundaryVector2(3),hipo2,localx3(3),hipo3,projCheck
      !integer(ip) :: ialone,idime,jdime,inode !node alone
      real(rp)    :: Residuo(3),alfa,gppre(1),beta,grpreN
      integer(ip) :: idime,inode
      real(rp) :: tmp, UnitOrthogonalVector(3)

      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      toleralce=1e-6
      
      if(elemStatus==0)then  
      
         !To integrate over the surface or line intersection
         call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
         call a%CutMesh%GetSurfaceIntersection(ielem,e,pdsurf)
         call a%CutMesh%GetNinters(ielem,nipoin)
         call a%CutMesh%GetExternalNormal(ielem,e,UnitOrthogonalVector)
         
               
         !Compute shape functions at the intersection points
         do jpoin = 1, nipoin         
            weigp(jpoin) = 1.0_rp/nipoin
         end do      
            
         !The rutine give the needed shape functions associated 
         call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)    
         
         !-----------------------------------------------------------------------------------------------
!             !Approximate imposition of boundary conditions in immersed boundary methods
            
            !Ramon-Joan Method         
         if(pdsurf>toleralce)then

            
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               !if(poinStatus==-1)then  
                  !elmat(:,inode,:,1:e%pnode) = 0.0_rp
                  !elrhs(:,inode) = 0.0_rp
                  do jnode = 1,e%pnode
                     do igaus = 1,e%pgaus
                        e%igaus = igaus                    
                        tmp=1.0_rp/(a%MatProp(1)%densi*LHSdtinv)
      
                        elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
                        + 100*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*(tmp/chale(1))*weigp(e%igaus)*pdsurf

                     enddo  
                  enddo                
               
               !end if         
            end do           
            
            !terms from the integrations by parts
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               if(poinStatus==1)then
                  do jnode = 1,e%pnode
                     do idime=1,e%ndime                  
                        do igaus = 1,e%pgaus
                           e%igaus = igaus
                           tmp=1.0_rp/(a%MatProp(1)%densi*LHSdtinv)
                           
                           elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) &
                              - e%shape(inode,e%igaus)*e%cartd(idime,jnode)*UnitOrthogonalVector(idime)*(tmp)*(weigp(e%igaus)*pdsurf)     
                        end do
                     enddo                  
                  enddo                
               
               end if         
            end do  
            
            !terms from the integrations by parts
            do inode=1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,poinStatus)            
               if(poinStatus==1)then
                  do idime=1,e%ndime                  
                     do igaus = 1,e%pgaus
                        e%igaus = igaus
                        tmp=1.0_rp/(a%MatProp(1)%densi*LHSdtinv)                

                           
                        elrhs(1,inode) = elrhs(1,inode) &
                           - e%shape(inode,e%igaus)*grpre(1,idime)*tmp*UnitOrthogonalVector(idime)*(weigp(e%igaus)*pdsurf)                             
                              
                              
                           
                     end do
                  enddo                  
               
               end if         
            end do           
         
         elseif(pdsurf<=toleralce)then

         end if
      
      end if
   
   end subroutine   
   
end module
