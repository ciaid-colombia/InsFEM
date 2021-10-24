module Mod_plcd_StrainGenerator 
   use typre
   use Mod_MPIObject
   use Mod_Listen
   implicit none
   private
   public AllocateGetStrain,GetVoigtsize, AllocateAddSphericalComponent, AllocateGetStressTensor, AllocateGetStrainTensor, AllocateGetStressVector, AllocateGetStrainVector
   
contains
   
   subroutine GetStrain3D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)

      
      vector(1) = tensor(1,1)
      vector(2) = tensor(2,2)
      vector(3) = tensor(3,3)
      vector(4) = tensor(1,2)+tensor(2,1)
      vector(5) = tensor(1,3)+tensor(3,1)
      vector(6) = tensor(2,3)+tensor(3,2)
   end subroutine
   
   subroutine GetStrain2D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
      
      
      vector(1) = tensor(1,1)
      vector(2) = tensor(2,2)
      vector(3) = tensor(1,2)+tensor(2,1)
   end subroutine
   
   subroutine AddSphericalComponent3D(value,vector)
      use typre
      implicit none
      real(rp) :: value
      real(rp) :: vector(6)

      vector(1:3) = vector(1:3) + value
   end subroutine
   
   subroutine AddSphericalComponent2D(value,vector)
      use typre
      implicit none
      real(rp) :: value
      real(rp) :: vector(3)

      vector(1:2) = vector(1:2) + value
   end subroutine
   
   subroutine GetAlmansiStrain3D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)

      vector(1) = tensor(1,1)-(0.5_rp)*tensor(1,1)*tensor(1,1)-(0.5_rp)*tensor(2,1)*tensor(2,1)-(0.5_rp)*tensor(3,1)*tensor(3,1)
      vector(2) = tensor(2,2)-(0.5_rp)*tensor(1,2)*tensor(1,2)-(0.5_rp)*tensor(2,2)*tensor(2,2)-(0.5_rp)*tensor(3,2)*tensor(3,2)
      vector(3) = tensor(3,3)-(0.5_rp)*tensor(1,3)*tensor(1,3)-(0.5_rp)*tensor(2,3)*tensor(2,3)-(0.5_rp)*tensor(3,3)*tensor(3,3)
      vector(4) = tensor(1,2)+tensor(2,1)-tensor(1,1)*tensor(1,2)-tensor(2,1)*tensor(2,2)-tensor(3,1)*tensor(3,2)
      vector(5) = tensor(1,3)+tensor(3,1)-tensor(1,1)*tensor(1,3)-tensor(2,1)*tensor(2,3)-tensor(3,1)*tensor(3,3)
      vector(6) = tensor(2,3)+tensor(3,2)-tensor(1,2)*tensor(1,3)-tensor(2,2)*tensor(2,3)-tensor(3,2)*tensor(3,3)
   end subroutine
   
   subroutine GetAlmansiStrain2D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
      
      vector(1) = tensor(1,1)-(0.5_rp)*tensor(1,1)*tensor(1,1)-(0.5_rp)*tensor(2,1)*tensor(2,1)
      vector(2) = tensor(2,2)-(0.5_rp)*tensor(1,2)*tensor(1,2)-(0.5_rp)*tensor(2,2)*tensor(2,2)
      vector(3) = tensor(1,2)+tensor(2,1)-tensor(1,1)*tensor(1,2)-tensor(2,1)*tensor(2,2)
   end subroutine
   
   subroutine GetGreenLagrangeStrain3D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)

      vector(1) = tensor(1,1)+(0.5_rp)*tensor(1,1)*tensor(1,1)+(0.5_rp)*tensor(2,1)*tensor(2,1)+(0.5_rp)*tensor(3,1)*tensor(3,1)
      vector(2) = tensor(2,2)+(0.5_rp)*tensor(1,2)*tensor(1,2)+(0.5_rp)*tensor(2,2)*tensor(2,2)+(0.5_rp)*tensor(3,2)*tensor(3,2)
      vector(3) = tensor(3,3)+(0.5_rp)*tensor(1,3)*tensor(1,3)+(0.5_rp)*tensor(2,3)*tensor(2,3)+(0.5_rp)*tensor(3,3)*tensor(3,3)
      vector(4) = tensor(1,2)+tensor(2,1)+tensor(1,1)*tensor(1,2)+tensor(2,1)*tensor(2,2)+tensor(3,1)*tensor(3,2)
      vector(5) = tensor(1,3)+tensor(3,1)+tensor(1,1)*tensor(1,3)+tensor(2,1)*tensor(2,3)+tensor(3,1)*tensor(3,3)
      vector(6) = tensor(2,3)+tensor(3,2)+tensor(1,2)*tensor(1,3)+tensor(2,2)*tensor(2,3)+tensor(3,2)*tensor(3,3)
   end subroutine
   
   subroutine GetGreenLagrangeStrain2D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
      
      vector(1) = tensor(1,1)+(0.5_rp)*tensor(1,1)*tensor(1,1)+(0.5_rp)*tensor(2,1)*tensor(2,1)
      vector(2) = tensor(2,2)+(0.5_rp)*tensor(1,2)*tensor(1,2)+(0.5_rp)*tensor(2,2)*tensor(2,2)
      vector(3) = tensor(1,2)+tensor(2,1)+tensor(1,1)*tensor(1,2)+tensor(2,1)*tensor(2,2)
   end subroutine
   
   subroutine GetStressTensor2D(vector,tensor)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
      
      tensor =  reshape([vector(1),vector(3),vector(3),vector(2)],[2,2])
   end subroutine
   
   subroutine GetStressTensor3D(vector,tensor)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)
      
      tensor =  reshape([ vector(1),vector(4),vector(5),vector(4),vector(2),vector(6),vector(5),vector(6),vector(3)],[3,3])
   end subroutine
   
   subroutine GetStrainTensor2D(vector,tensor)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
      
      tensor =  reshape([vector(1),0.5_rp*vector(3),0.5_rp*vector(3),vector(2)],[2,2])
   end subroutine
   
   subroutine GetStrainTensor3D(vector,tensor)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)
      
      tensor =  reshape([vector(1),0.5_rp*vector(4),0.5_rp*vector(5),0.5_rp*vector(4),vector(2),0.5_rp*vector(6),&
      0.5_rp*vector(5),0.5_rp*vector(6),vector(3)],[3,3])
   end subroutine
   
   subroutine GetStressVector2D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3,1)
      
      vector =  reshape([tensor(1,1),tensor(2,2),tensor(1,2)],[3,1])
   end subroutine
   
   subroutine GetStressVector3D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6,1)
      
      vector =  reshape([tensor(1,1),tensor(2,2),tensor(3,3),tensor(1,2),tensor(1,3),tensor(2,3)],[6,1])
   end subroutine
   
   subroutine GetStrainVector2D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3,1)
      
      vector =  reshape([tensor(1,1), tensor(2,2), 2.0_rp*tensor(1,2)],[3,1])
   end subroutine
   
   subroutine GetStrainVector3D(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6,1)
      
      vector =  reshape([tensor(1,1),tensor(2,2),tensor(3,3),2.0_rp*tensor(1,2),2.0_rp*tensor(1,3),2.0_rp*tensor(2,2)], &
      [6,1])
   end subroutine
   
   subroutine GetVoigtsize(ndime,vsize)
      integer(ip) :: ndime,vsize
      
      if (ndime == 2) then
         vsize = 3
      elseif (ndime == 3) then
         vsize = 6
      endif 
   end subroutine
   
   subroutine AllocateGetStrain(ndime,GetStrain,LargeStrainsflag)
      integer(ip) :: ndime
      procedure(), pointer :: GetStrain
      integer(ip) :: LargeStrainsflag
      if (LargeStrainsflag == 0) then
         if (ndime == 2) then
            GetStrain => GetStrain2D
         elseif (ndime == 3) then
            GetStrain => GetStrain3D
         endif
      elseif (LargeStrainsflag == 1) then
         if (ndime == 2) then
            GetStrain => GetAlmansiStrain2D
         elseif (ndime == 3) then
            GetStrain => GetAlmansiStrain3D
         endif
      elseif (LargeStrainsflag == 2) then
         if (ndime == 2) then
            GetStrain => GetGreenLagrangeStrain2D
         elseif (ndime == 3) then
            GetStrain => GetGreenLagrangeStrain3D
         endif
      endif
   end subroutine
   
   subroutine AllocateAddSphericalComponent(ndime,AddSphericalComponent)
      integer(ip) :: ndime
      procedure(), pointer :: AddSphericalComponent
      if (ndime == 2) then
         AddSphericalComponent => AddSphericalComponent2D
      elseif (ndime == 3) then
         AddSphericalComponent => AddSphericalComponent3D
      endif
   end subroutine
   
   subroutine AllocateGetStressTensor(ndime,GetStressTensor)
      integer(ip) :: ndime
      procedure(), pointer :: GetStressTensor
      if (ndime == 2) then
         GetStressTensor => GetStressTensor2D
      elseif (ndime == 3) then
         GetStressTensor => GetStressTensor3D
      endif
   end subroutine
   
   subroutine AllocateGetStrainTensor(ndime,GetStrainTensor)
      integer(ip) :: ndime
      procedure(), pointer :: GetStrainTensor
      if (ndime == 2) then
         GetStrainTensor => GetStrainTensor2D
      elseif (ndime == 3) then
         GetStrainTensor => GetStrainTensor3D
      endif
   end subroutine
   
   subroutine AllocateGetStressVector(ndime,GetStressVector)
      integer(ip) :: ndime
      procedure(), pointer :: GetStressVector
      if (ndime == 2) then
         GetStressVector => GetStressVector2D
      elseif (ndime == 3) then
         GetStressVector => GetStressVector3D
      endif
   end subroutine
   
   subroutine AllocateGetStrainVector(ndime,GetStrainVector)
      integer(ip) :: ndime
      procedure(), pointer :: GetStrainVector
      if (ndime == 2) then
         GetStrainVector => GetStrainVector2D
      elseif (ndime == 3) then
         GetStrainVector => GetStrainVector3D
      endif
   end subroutine
   
end module
