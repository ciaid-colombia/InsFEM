module Mod_RVEpropDriver
    use typre
    use MPI
    use Mod_BroadCastBuffer
    use Mod_Listen
    use Mod_caseVariables
    use Mod_PhysicalProblemDriver
    use Mod_PhysicalProblem
    use Mod_DistributedContainer 
    use Mod_DC_rp
    use Mod_InChannel
    use Mod_MasterVariables
    use Mod_PLCD
    use Mod_iofile
   implicit none

   private

   public :: RVEpropDriver, RVEpropDriver_Const
   
   type(PLCDProblem), pointer :: PLCD => NULL()

   type, extends(DriverInterface) :: RVEpropDriver 
      
      type(PLCDProblem) :: PLCD
      
      integer(ip) :: propfile

!       !Couplings 
!       type(InChannelGroupInfo) :: InNstincInfo
     

contains

      procedure :: Lev1Reapro             => RVEprop_Lev1Reapro
      procedure :: Lev1ReaproMPI          => RVEprop_Lev1ReaproMPI
      procedure :: SetOutChannels         => RVEprop_SetOutChannels
      procedure :: SetInChannels          => RVEprop_SetInChannels
      procedure :: UpdateChannels         => RVEprop_UpdateChannels
      procedure :: Turnon                 => RVEprop_Turnon
      procedure :: GetTimeStep            => RVEprop_GetTimeStep
      procedure :: SetTimeStep            => RVEprop_SetTimeStep
      procedure :: GetRefinementCriteria  => RVEprop_GetRefinementCriteria
      procedure :: PreRefine              => RVEprop_PreRefine
      procedure :: Refine                 => RVEprop_Refine
      procedure :: Rebalance              => RVEprop_Rebalance
      procedure :: Begste                 => RVEprop_Begste
      procedure :: MeshProjections        => RVEprop_MeshProjections
      procedure :: MeshAdvections         => RVEprop_MeshAdvections
      procedure :: Doiter                 => RVEprop_Doiter
      procedure :: Convergence            => RVEprop_Convergence
      procedure :: Endste                 => RVEprop_Endste
      procedure :: Output                 => RVEprop_Output
      procedure :: Turnof                 => RVEprop_Turnof
      procedure :: WriteTimes             => RVEprop_WriteTimes
      procedure :: GetPhysicalProblem     => RVEprop_GetPhysicalProblem
      
   end type RVEpropDriver

   interface RVEpropDriver_Const
      procedure constructor
   end interface 

contains

   function constructor()
       class(RVEpropDriver), pointer :: constructor

       allocate(constructor)
   end function constructor


   subroutine RVEprop_Lev1Reapro(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb


   end subroutine RVEprop_Lev1Reapro
   
   subroutine RVEprop_Lev1ReaproMPI(a,c,bb)
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      type(BroadCastBuffer) :: bb
      
   end subroutine RVEprop_Lev1ReaproMPI
   
   subroutine RVEprop_SetOutChannels(a,task)
      class(RVEpropDriver) :: a
      character(*) :: task

   end subroutine RVEprop_SetOutChannels
   
   subroutine RVEprop_SetInChannels(a)
      class(RVEpropDriver) :: a

   end subroutine RVEprop_SetInChannels
      
   subroutine RVEprop_UpdateChannels(a)
      implicit none
      class(RVEpropDriver) :: a
      class(DistributedContainer), pointer :: myDC => NULL()
      real(rp), pointer :: trac(:,:) => NULL()

   end subroutine RVEprop_UpdateChannels
   
   subroutine RVEprop_Turnon(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      
      call point_RVEprop(a)
      call a%PLCD%SetComputeMean(.true.)
      
      if (a%MPIrank == a%MPIroot) call iofile(zero,a%propfile,'MechProp.txt','PropertiesRVE')
      
      !PhysicalProblem operations
      call physical_Turnon(a,c,a%PLCD)
      
      Call a%PLCD%DeallocStages
      Call a%PLCD%SetStages

   end subroutine RVEprop_Turnon
   
   subroutine point_RVEprop(a)
      implicit none
      class(RVEpropDriver), target :: a
      
      PLCD => a%PLCD
   end subroutine point_RVEprop
   
   
   subroutine RVEprop_GetTimeStep(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetTimeStep(a,c,a%PLCD)
   end subroutine RVEprop_GetTimeStep
   
   subroutine RVEprop_SetTimeStep(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_SetTimeStep(a,c,a%PLCD)
   end subroutine RVEprop_SetTimeStep
   
   subroutine RVEprop_GetRefinementCriteria(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_GetRefinementCriteria(a,c,a%PLCD)
   end subroutine RVEprop_GetRefinementCriteria
   
   subroutine RVEprop_PreRefine(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_PreRefine(a,c,a%PLCD)
   end subroutine RVEprop_PreRefine
   
   subroutine RVEprop_Refine(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Refine(a,c,a%PLCD)
   end subroutine RVEprop_Refine
   
   subroutine RVEprop_Rebalance(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Rebalance(a,c,a%PLCD)
   end subroutine RVEprop_Rebalance
   
   subroutine RVEprop_Begste(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
 
      !PhysicalProblem operations
      call physical_Begste(a,c,a%PLCD)
   end subroutine RVEprop_Begste
   
   subroutine RVEprop_MeshProjections(a,c,itask)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshProjections(a,c,a%PLCD,itask)
   end subroutine RVEprop_MeshProjections
   
   subroutine RVEprop_MeshAdvections(a,c,itask)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      integer(ip)         :: itask
  
      !PhysicalProblem operations
      call physical_MeshAdvections(a,c,a%PLCD,itask)
   end subroutine RVEprop_MeshAdvections
   
   subroutine RVEprop_Doiter(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      real(rp), pointer :: StressMean(:) => NULL()
  
      !PhysicalProblem operations
      call physical_Doiter(a,c,a%PLCD)
      
   end subroutine RVEprop_Doiter
   
   subroutine RVEprop_Convergence(a,c,glres)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      real(rp)            :: glres
  
      !PhysicalProblem operations
      call physical_Convergence(a,c,a%PLCD,glres)
   end subroutine RVEprop_Convergence
   
   subroutine RVEprop_Endste(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Endste(a,c,a%PLCD)
   end subroutine RVEprop_Endste

   subroutine RVEprop_Output(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
  
      !PhysicalProblem operations
      call physical_Output(a,c,a%PLCD)
   end subroutine RVEprop_Output
   
   subroutine RVEprop_Turnof(a,c)
      implicit none
      class(RVEpropDriver) :: a
      type(caseVariables) :: c
      real(rp), pointer :: ConstTensor(:,:) => NULL()
      integer(ip) :: i,j
      
      call a%PLCD%GetMeanConst(ConstTensor)
      
      if (a%MPIrank == a%MPIroot) call PrintRVEprop(a,ConstTensor)
      !PhysicalProblem operations
      call physical_Turnof(a,c,a%PLCD)
      
      if (a%MPIrank == a%MPIroot) call iofile(two,a%propfile,'MechProp.txt','PropertiesRVE')
      
   end subroutine RVEprop_Turnof
   
   subroutine RVEprop_WriteTimes(a)
      implicit none
      class(RVEpropDriver) :: a
      
      !PhysicalProblem operations
      call physical_WriteTimes(a,a%PLCD)
   end subroutine RVEprop_WriteTimes
   
   subroutine RVEprop_GetPhysicalProblem(a,PhysicalPro)
      implicit none
      class(RVEpropDriver), target :: a
      class(PhysicalProblem), pointer :: PhysicalPro
      
      PhysicalPro => a%PLCD
   end subroutine

   subroutine PrintRVEprop(a,C)
      implicit none
      class(RVEpropDriver) :: a
      real(rp), pointer :: C(:,:)
      real(rp) :: maxV, InvC(6,6)
      integer(ip) :: i,j,nn,ndime
      character(len=70) :: fmt,fnn,fmt1,fmt2,fmt3,fmt4,fmt5
      
      Real(rp) :: Exx,Eyy,Ezz,Gxy,Gxz,Gyz,Pxy,Pyx,Pxz,Pzx,Pyz,Pzy
      
      ndime = 6
      !a%PLCD%Mesh%GetNdime(ndime)
      
      maxV = maxval(C(1:ndime,1:ndime))
      nn = 6
      
      do i = 1,20
        if((maxV/10**i)<1) then
         nn = i + 4
         exit 
        end if
      end do
      write(fnn,'(I2)') nn
      
      write(a%propfile,*)'Estimated Mechanical Properties of the RVE'
      write(a%propfile,*)
      write(a%propfile,*)'The properties have the same metric unid of the constituint RVE materials'
      write(a%propfile,*)

      if (ndime==6) then
      
         Call InvertMatrix(C(1:ndime,1:ndime),InvC(1:ndime,1:ndime))

         Exx = 1/InvC(1,1)   ;   Eyy = 1/InvC(2,2) ;   Ezz = 1/InvC(3,3)
         Gxy = 1/InvC(4,4)   ;   Gxz = 1/InvC(5,5) ;   Gyz = 1/InvC(6,6)
         Pxy = -InvC(1,2)*Eyy    ;   Pyx = -InvC(2,1)*Exx
         Pxz = -InvC(1,3)*Ezz    ;   Pzx = -InvC(3,1)*Exx
         Pyz = -InvC(2,3)*Ezz    ;   Pzy = -InvC(3,2)*Eyy
      
         fmt1 = "(' Exx = ',F"//trim(fnn)//".2,2X,'Eyy = ',F"//trim(fnn)//".2,2X,'Ezz = ',F"//trim(fnn)//".2)"
         fmt2 = "(' Gxy = ',F"//trim(fnn)//".2,2X,'Gxz = ',F"//trim(fnn)//".2,2X,'Gyz = ',F"//trim(fnn)//".2)"
         fmt3 = "(' Pxy = ',F7.4,2X,'Pyx = ',F7.4)"
         fmt4 = "(' Pxz = ',F7.4,2X,'Pzx = ',F7.4)"
         fmt5 = "(' Pyz = ',F7.4,2X,'Pzy = ',F7.4)"
         
         write(a%propfile,fmt1) Exx,Eyy,Ezz
         write(a%propfile,fmt2) Gxy,Gxz,Gyz
         write(a%propfile,fmt3) Pxy,Pyx
         write(a%propfile,fmt4) Pxz,Pzx
         Write(a%propfile,fmt5) Pyz,Pzy
      
      else if(ndime==3) then
      
         Call InvertMatrix(C(1:ndime,1:ndime),InvC(1:ndime,1:ndime))

         Exx = 1/InvC(1,1)   ;   Eyy = 1/InvC(2,2)   ;   Gxy = 1/InvC(3,3)
         Pxy = -InvC(1,2)*Eyy    ;   Pyx = -InvC(2,1)*Exx
      
         fmt1 = "(' Exx = ',F"//trim(fnn)//".2,2X,'Eyy = ',F"//trim(fnn)//".2,2X,'Gxy = ',F"//trim(fnn)//".2)"
         fmt2 = "(' Pxy = ',F7.4,2X,'Pyx = ',F7.4)"
         
         write(a%propfile,fmt1) Exx,Eyy,Gxy
         write(a%propfile,fmt2) Pxy,Pyx
     
      end if

      write(a%propfile,*)
      if (ndime==6) then
        write(a%propfile,*) 'Constitutive Tensor of the RVE (x,y,z,xy,xz,yz)'
        fmt = "(6(F"//trim(fnn)//".2,2X))"
      else if(ndime==3) then
        write(a%propfile,*) 'Constitutive Tensor of the RVE (x,y,xy)'
        fmt = "(3(F"//trim(fnn)//".2,2X))"
      end if
      write(a%propfile,*)
      
      do i = 1,ndime
         write(a%propfile,fmt) (C(i,j), j=1,ndime)
  
      end do
   
   end subroutine PrintRVEprop
   
  Subroutine InvertMatrix(Matrix,IMatrix)
    Implicit none

    Real(rp), Intent(in), Dimension(:,:) :: Matrix
    Real(rp), Intent(out), Dimension(:,:) :: IMatrix
    Real(rp), Allocatable :: AuxA(:,:)
    Real(rp), Parameter :: zero = 1.0e-30
    Real(rp) :: d, f
    Integer(ip) :: n, i, j

    n = size(Matrix,1)

    Allocate(AuxA(n,n))

    AuxA = Matrix   ;   IMatrix = 0.0d0

    do i=1,n
      if(abs(AuxA(i,i)) < zero) Call RunEnd('A zero value on the matrix diagonal, it is not possible get the invert')
      IMatrix(i,i) = 1.0d0
    
    end do

    Do i = 1,n
      d = AuxA(i,i)
      
      Do j = 1,n
        AuxA(i,j) = AuxA(i,j)/d
        IMatrix(i,j) = IMatrix(i,j)/d

      End do

      Do j = i+1,n
        f = AuxA(j,i)
        AuxA(j,:) = AuxA(j,:) - AuxA(i,:)*f
        IMatrix(j,:) = IMatrix(j,:) - IMatrix(i,:)*f

      End do

      Do j = i-1,1,-1
        f = AuxA(j,i)
        AuxA(j,:) = AuxA(j,:) - AuxA(i,:)*f
        IMatrix(j,:) = IMatrix(j,:) - IMatrix(i,:)*f

      End do

    End do

    deallocate(AuxA)

  end subroutine
   
end module 
