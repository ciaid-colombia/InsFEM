   subroutine DeallocStages(a)
      use typre
      use Mod_PLCD
      implicit none
      class(PLCDProblem) :: a
      integer(ip) :: ndime,npoin,istage,nsubstages,nboun
      
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         call a%Mesh%GetNboun(nboun)
                       
         do istage = 1,a%NumberOfStages
            call a%Memor%dealloc(ndime,npoin,a%Stages(istage)%NodalForces,'NodalForces','plcd_memall')
            call a%Memor%dealloc(ndime,npoin,a%Stages(istage)%kfl_fixno,'kfl_fixno','plcd_memall')
            call a%Memor%dealloc(ndime,npoin,a%Stages(istage)%bvess,'bvess','plcd_memall')
            call a%Memor%dealloc(nboun,a%Stages(istage)%kfl_fixbo,'kfl_fixbo','plcd_memall')
            call a%Memor%dealloc(nboun,a%Stages(istage)%bvnat,'bvnat','plcd_memall')
            call a%Memor%dealloc(nboun,a%Stages(istage)%kfl_funbo,'kfl_funbo','plcd_memall')
            
            nsubstages = a%Stages(istage)%NumberOfSubStages
            deallocate(a%Stages(istage)%Substages)
            call a%Memor%deallocObj(0,'SubStages','plcd_memall',1*nsubstages)
            !Will only work as long as no allocatables in substage are present
         enddo 
      deallocate(a%Stages)
      call a%Memor%deallocObj(0,'Stages','plcd_memall',1*a%NumberOfStages)
   
   end subroutine

   subroutine SetStages(a)
      use typre
      use Mod_PLCD
      use Mod_plcd_StrainGenerator
      implicit none
      class(PLCDProblem), target :: a      
      integer(ip) :: vsize,ndime,npoin,istage,ipoin,ibopo
      real(rp) :: coord(3), GD(3,3), time
      real(rp), pointer :: coordpointer(:)
      real(rp), pointer :: exnor(:,:)
         
         call a%Mesh%GetNdime(ndime)
         call a%Mesh%GetNpoin(npoin)
         
         call GetVoigtsize(ndime,vsize)
         a%NumberOfStages = vsize
      
         allocate(a%Stages(a%NumberOfStages))
         call a%Memor%allocObj(0,'Stages','plcd_memall',1*a%NumberOfStages)
         time = 0.0_rp
                 
         do istage = 1,a%NumberOfStages
         
            a%Stages(istage)%NumberOfSubstages = 1
            allocate(a%Stages(istage)%Substages(a%Stages(istage)%NumberOfSubstages))
            call a%Memor%allocObj(0,'SubStages','plcd_memall',1*a%Stages(istage)%NumberOfSubstages)
            !Setup stages and substages
            a%Stages(istage)%Substages(1)%TimeInterval = 1
            a%Stages(istage)%Substages(1)%IniTime = time
            time = time + a%Stages(istage)%Substages(1)%TimeInterval
            a%Stages(istage)%Substages(1)%EndTime = time
 
            call a%Memor%alloc(ndime,npoin,a%Stages(istage)%NodalForces,'NodalForces','plcd_memall')
            call a%Memor%alloc(ndime,npoin,a%Stages(istage)%kfl_fixno,'kfl_fixno','plcd_memall')
            call a%Memor%alloc(ndime,npoin,a%Stages(istage)%bvess,'bvess','plcd_memall')
            
            !Default is no boundary condition
            a%Stages(istage)%kfl_fixno = -1
            
            call GetF(istage,ndime,GD)
         
            do ipoin = 1, npoin

               call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
               if(ibopo>0) then
               
                  coord = 0
                  call a%Mesh%GetPointCoord(ipoin,coordpointer)
                  coord(1:ndime) = coordpointer
                  a%Stages(istage)%kfl_fixno(1:a%ndofbc,ipoin) = 1
                  a%Stages(istage)%bvess(1:a%ndofbc,ipoin) = matmul(GD(1:ndime,1:ndime),coord(1:ndime))
               endif
            enddo
         
         enddo
         
         a%cs => a%Stages(1)
         a%css => a%Stages(1)%Substages(1)
    
   end subroutine
   
   Subroutine GetF(istage,ndime,GD)
   use typre
   implicit none
   
   integer(ip) :: istage,ndime
   real(rp) :: GD(3,3)
   
   GD = 0
     
   Select case(istage)
   
      case(1) ! x direction
         GD(1,1) = 1
      
      case(2) ! y direction
         GD(2,2) = 1
      
      case(3) ! Z or XY direction
         Select case(ndime)
            
            case(3) ! 3D dimension 
               GD(3,3) = 1
            
            case(2) ! 2D dimension
               GD(1,2) = 0.5
               GD(2,1) = 0.5
            
            end select
      case(4) ! XY direction
      
         GD(1,2) = 0.5
         GD(2,1) = 0.5
      
      case(5) ! XZ direction
         
         GD(1,3) = 0.5
         GD(3,1) = 0.5
      
      case(6) ! YZ direction
      
         GD(2,3) = 0.5
         GD(3,2) = 0.5
      
   end select
     
   end subroutine