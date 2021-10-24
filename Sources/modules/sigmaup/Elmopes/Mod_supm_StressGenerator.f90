module Mod_supm_StressGenerator 
   use typre
   use Mod_MPIObject
   use Mod_Listen
   use Mod_supm_MatrixVector
   use Mod_Element
   implicit none
   
contains
   
!------------------------------ 
! scalar per tensor ht
!------------------------------

   subroutine GetStress2D(tensor,vector)
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
       
      vector(1) = tensor(1,1)
      vector(2) = tensor(2,2)
      vector(3) = tensor(1,2)+tensor(2,1)
   end subroutine
   
   
     subroutine GetStress2D_sym(tensor,vector)
      implicit none
      real(rp) :: tensor(2,2)
      real(rp) :: vector(3)
       
      vector(1) = tensor(1,1)
      vector(2) = tensor(2,2)
      vector(3) = tensor(1,2)
   end subroutine
   
   subroutine AddSphericalComponent2D(value,vector)
      use typre
      implicit none
      real(rp) :: value
      real(rp) :: vector(3)

      vector(1:2) = vector(1:2) + value
   end subroutine
 
   subroutine GetStress3D(tensor,vector)
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)

       
      vector(1) = tensor(1,1)
      vector(2) = tensor(2,2)
      vector(3) = tensor(3,3)
      vector(4) = tensor(2,3)+tensor(3,2)
      vector(5) = tensor(1,3)+tensor(3,1)
      vector(6) = tensor(1,2)+tensor(2,1)
   end subroutine
   
   
    subroutine GetStress3D_sym(tensor,vector)
      implicit none
      real(rp) :: tensor(3,3)
      real(rp) :: vector(6)

       
      vector(1) = tensor(1,1)
      vector(2) = tensor(2,2)
      vector(3) = tensor(3,3)
      vector(4) = tensor(2,3)
      vector(5) = tensor(1,3)
      vector(6) = tensor(1,2)
   end subroutine
   
   subroutine AddSphericalComponent3D(value,vector)
      implicit none
      real(rp) :: value
      real(rp) :: vector(6)

      vector(1:3) = vector(1:3) + value
   end subroutine
   
   
  subroutine GetStress(ndime,tensor,vector)
      integer(ip) :: ndime
      real(rp)    :: vector(:), tensor(:,:)
      
      if (ndime == 2) call GetStress2D(tensor,vector)
      if (ndime == 3) call GetStress3D(tensor,vector)
  end subroutine
  
  
    subroutine GetStress_sym(ndime,tensor,vector)
      integer(ip) :: ndime
      real(rp)    :: vector(:), tensor(:,:)
      
      if (ndime == 2) call GetStress2D_sym(tensor,vector)
      if (ndime == 3) call GetStress3D_sym(tensor,vector)
      
  end subroutine
    
   
!------------------------------ 
! two tensors product st
!------------------------------
   subroutine GetStress2D_tt(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(2,2), tensor2(3,3)
      
      tensor2(1,1)=tensor1(1,1)
      tensor2(1,2)=0
      tensor2(1,3)=tensor1(1,2)
      tensor2(2,1)=0
      tensor2(2,2)=tensor1(2,2)
      tensor2(2,3)=tensor1(2,1)
      tensor2(3,1)=tensor1(2,1)
      tensor2(3,2)=tensor1(1,2)
      tensor2(3,3)=tensor1(1,1)+tensor1(2,2)      
    end subroutine  
   
   
   subroutine AddSphericalComponent2D_tt(value,tensor)
      implicit none
      real(rp) :: value
      real(rp) :: tensor(3,3)

      tensor(1,1) = tensor(1,1) + value
      tensor(2,2) = tensor(2,2) + value
      tensor(3,3) = tensor(3,3) + 2.0_rp*value
   end subroutine
   
   
   subroutine GetAdvStress2D_tt(auxPTT,vector,tensor)
      implicit none
      real(rp) :: tensor(3,3), vector(3), traza
      integer :: auxPTT
      traza=vector(1)+vector(2)
      
      tensor(1,1)= (1_ip + auxPTT)*vector(1)+(1_ip - auxPTT)*traza
      tensor(1,2)= (1_ip - auxPTT)*vector(1)
      tensor(1,3)= auxPTT*2.0_rp*vector(3)
      tensor(2,1)= (1_ip - auxPTT)*vector(2)
      tensor(2,2)= (1_ip + auxPTT)*vector(2)+(1_ip - auxPTT)*traza
      tensor(2,3)= auxPTT*2.0_rp*vector(3)
      tensor(3,1)= 2.0_rp*vector(3)
      tensor(3,2)= 2.0_rp*vector(3)
      tensor(3,3)= (1_ip + auxPTT)*(vector(1)+vector(2))+(1_ip - auxPTT)*traza
      
   end subroutine
   
   
   subroutine GetAdvStress2D_tt_B(auxPTT,vector,tensor)
      implicit none
      real(rp) :: tensor(3,3), vector(3), traza
      integer :: auxPTT
      traza=vector(1)+vector(2)
      
      tensor(1,1)= -auxPTT*vector(1)-(1_ip - auxPTT)*traza
      tensor(1,2)= 0
      tensor(1,3)= -auxPTT*vector(3)
      tensor(2,1)= 0
      tensor(2,2)= -auxPTT*vector(2)-(1_ip - auxPTT)*traza
      tensor(2,3)= -auxPTT*vector(3)
      tensor(3,1)= -auxPTT*vector(3)
      tensor(3,2)= -auxPTT*vector(3)
      tensor(3,3)= 2.0_rp*(-auxPTT*(vector(1)+vector(2))-(1_ip - auxPTT)*traza)
      
   end subroutine
    
    subroutine GetStress3D_tt(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(3,3), tensor2(6,6)
      
      tensor2(1,1)=tensor1(1,1)
      tensor2(1,2)=0
      tensor2(1,3)=0
      tensor2(1,4)=0
      tensor2(1,5)=tensor1(1,3)
      tensor2(1,6)=tensor1(1,2)
      
      tensor2(2,1)=0
      tensor2(2,2)=tensor1(2,2)
      tensor2(2,3)=0
      tensor2(2,4)=tensor1(2,3)
      tensor2(2,5)=0
      tensor2(2,6)=tensor1(2,1)
      
      tensor2(3,1)=0
      tensor2(3,2)=0
      tensor2(3,3)=tensor1(3,3)
      tensor2(3,4)=tensor1(3,2)
      tensor2(3,5)=tensor1(3,1)
      tensor2(3,6)=0
      
      tensor2(4,1)=0
      tensor2(4,2)=tensor1(3,2)
      tensor2(4,3)=tensor1(2,3)
      tensor2(4,4)=tensor1(3,3)+tensor1(2,2)
      tensor2(4,5)=tensor1(2,1)
      tensor2(4,6)=tensor1(3,1) 
      
      
      tensor2(5,1)=tensor1(3,1)
      tensor2(5,2)=0
      tensor2(5,3)=tensor1(1,3)
      tensor2(5,4)=tensor1(1,2)
      tensor2(5,5)=tensor1(1,1)+tensor1(3,3)
      tensor2(5,6)=tensor1(3,2) 
      
      
      tensor2(6,1)=tensor1(2,1)
      tensor2(6,2)=tensor1(1,2)
      tensor2(6,3)=0
      tensor2(6,4)=tensor1(1,3)
      tensor2(6,5)=tensor1(2,3)
      tensor2(6,6)=tensor1(1,1)+tensor1(2,2) 
  
    end subroutine  
    
       
   subroutine AddSphericalComponent3D_tt(value,tensor)
      implicit none
      real(rp) :: value
      real(rp) :: tensor(6,6)

      tensor(1,1) = tensor(1,1) + value
      tensor(2,2) = tensor(2,2) + value
      tensor(3,3) = tensor(3,3) + value
      tensor(4,4) = tensor(4,4) + 2.0_rp*value
      tensor(5,5) = tensor(5,5) + 2.0_rp*value
      tensor(6,6) = tensor(6,6) + 2.0_rp*value
   end subroutine
   
   subroutine GetAdvStress3D_tt(auxPTT,vector,tensor)                               
      implicit none
      real(rp) :: tensor(6,6), vector(6), traza
      integer :: auxPTT
      !if auxPTT =1 is OFF, else auxPTT=0 is ON
      traza=vector(1)+vector(2)+vector(3)
      
      tensor(1,1)=(1_ip + auxPTT)*vector(1)+(1_ip - auxPTT)*traza
      tensor(1,2)=(1 - auxPTT)*vector(1)
      tensor(1,3)=(1 - auxPTT)*vector(1)
      tensor(1,4)=0
      tensor(1,5)=auxPTT*2.0_rp*vector(5)
      tensor(1,6)=auxPTT*2.0_rp*vector(6)
      
      tensor(2,1)=(1 - auxPTT)*vector(2)
      tensor(2,2)=(1_ip + auxPTT)*vector(2)+(1_ip - auxPTT)*traza
      tensor(2,3)=(1 - auxPTT)*vector(2)
      tensor(2,4)=auxPTT*2.0_rp*vector(4)
      tensor(2,5)=0
      tensor(2,6)=auxPTT*2.0_rp*vector(6)
      
      tensor(3,1)=(1 - auxPTT)*vector(3)
      tensor(3,2)=(1 - auxPTT)*vector(3)
      tensor(3,3)=(1_ip + auxPTT)*vector(3)+(1_ip - auxPTT)*traza
      tensor(3,4)=auxPTT*2.0_rp*vector(4)
      tensor(3,5)=auxPTT*2.0_rp*vector(5)
      tensor(3,6)=0
      
      tensor(4,1)=(1 - auxPTT)*2.0_rp*vector(4)
      tensor(4,2)=2.0_rp*vector(4)
      tensor(4,3)=2.0_rp*vector(4)
      tensor(4,4)=(1_ip + auxPTT)*(vector(2)+vector(3))+(1_ip - auxPTT)*traza
      tensor(4,5)=auxPTT*2.0_rp*vector(6)
      tensor(4,6)=auxPTT*2.0_rp*vector(5)

      tensor(5,1)=2.0_rp*vector(5)
      tensor(5,2)=(1 - auxPTT)*2.0_rp*vector(5)
      tensor(5,3)=2.0_rp*vector(5)
      tensor(5,4)=auxPTT*2.0_rp*vector(6)
      tensor(5,5)=(1_ip + auxPTT)*(vector(1)+vector(3))+(1_ip - auxPTT)*traza
      tensor(5,6)=auxPTT*2.0_rp*vector(4)
 
      tensor(6,1)=2.0_rp*vector(6)
      tensor(6,2)=2.0_rp*vector(6)
      tensor(6,3)=(1 - auxPTT)*2.0_rp*vector(6)
      tensor(6,4)=auxPTT*2.0_rp*vector(5)
      tensor(6,5)=auxPTT*2.0_rp*vector(4)
      tensor(6,6)=(1_ip + auxPTT)*(vector(1)+vector(2))+(1_ip - auxPTT)*traza

   end subroutine
   
   
   subroutine GetAdvStress3D_tt_B(auxPTT,vector,tensor)
      implicit none
      real(rp) :: tensor(6,6), vector(6), traza
      integer :: auxPTT
      !if auxPTT =1 is OFF, else auxPTT=0 is ON
      traza=vector(1)+vector(2)+vector(3)
      
      tensor(1,1)=-auxPTT*vector(1)-(1_ip - auxPTT)*traza
      tensor(1,2)=0
      tensor(1,3)=0
      tensor(1,4)=0
      tensor(1,5)=-auxPTT*vector(5)
      tensor(1,6)=-auxPTT*vector(6)
      
      tensor(2,1)=0
      tensor(2,2)=-auxPTT*vector(2)-(1_ip - auxPTT)*traza
      tensor(2,3)=0
      tensor(2,4)=-auxPTT*vector(4)
      tensor(2,5)=0
      tensor(2,6)=-auxPTT*vector(6)
      
      tensor(3,1)=0
      tensor(3,2)=0
      tensor(3,3)=-auxPTT*vector(3)-(1_ip - auxPTT)*traza
      tensor(3,4)=-auxPTT*vector(4)
      tensor(3,5)=-auxPTT*vector(5)
      tensor(3,6)=0
      
      tensor(4,1)=0
      tensor(4,2)=-auxPTT*vector(4)
      tensor(4,3)=-auxPTT*vector(4)
      tensor(4,4)=-auxPTT*(vector(2)+vector(3))-(1_ip - auxPTT)*traza
      tensor(4,5)=-auxPTT*vector(6)
      tensor(4,6)=-auxPTT*vector(5)

      tensor(5,1)=-auxPTT*vector(5)
      tensor(5,2)=0
      tensor(5,3)=-auxPTT*vector(5)
      tensor(5,4)=-auxPTT*vector(6)
      tensor(5,5)=-auxPTT*(vector(1)+vector(3))-(1_ip - auxPTT)*traza
      tensor(5,6)=-auxPTT*vector(4)
 
      tensor(6,1)=-auxPTT*vector(6)
      tensor(6,2)=-auxPTT*vector(6)
      tensor(6,3)=0
      tensor(6,4)=-auxPTT*vector(5)
      tensor(6,5)=-auxPTT*vector(4)
      tensor(6,6)=-auxPTT*(vector(1)+vector(2))-(1_ip - auxPTT)*traza

   end subroutine
   
   
     subroutine GetStress_tt(ndime,tensor1,tensor2)
      integer(ip) :: ndime
      real(rp)    :: tensor1(ndime,ndime), tensor2(:,:)
      
      if (ndime == 2) call GetStress2D_tt(tensor1,tensor2)
      if (ndime == 3) call GetStress3D_tt(tensor1,tensor2)
   end subroutine
   
    subroutine AddSphericalComponent_tt(ndime,value,tensor1)
      integer(ip) :: ndime
      real(rp)    :: tensor1(:,:), value
      
      if (ndime == 2) call AddSphericalComponent2D_tt(value,tensor1)
      if (ndime == 3) call AddSphericalComponent3D_tt(value,tensor1)
   end subroutine
   
    subroutine GetAdvStress_tt(ndime,auxPTT,vector,tensor)
      integer(ip) :: ndime, auxPTT
      real(rp)    :: vector(:), tensor(:,:)
      
      if (ndime == 2) call GetAdvStress2D_tt(auxPTT,vector,tensor)
      if (ndime == 3) call GetAdvStress3D_tt(auxPTT,vector,tensor)
   end subroutine
   
   
    subroutine GetAdvStress_tt_B(ndime,auxPTT,vector,tensor)
      integer(ip) :: ndime, auxPTT
      real(rp)    :: vector(:), tensor(:,:)
      
      if (ndime == 2) call GetAdvStress2D_tt_B(auxPTT,vector,tensor)
      if (ndime == 3) call GetAdvStress3D_tt_B(auxPTT,vector,tensor)
   end subroutine
   
!------------------------------ 
! vector per tensor ut
!------------------------------   
   
   subroutine GetStress2D_ut(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(2,2), tensor2(3,2)
      
      tensor2(1,1)=tensor1(1,1)
      tensor2(1,2)=0
      tensor2(2,1)=0
      tensor2(2,2)=tensor1(2,2)
      tensor2(3,1)=tensor1(2,2)
      tensor2(3,2)=tensor1(1,1)    
   end subroutine  
   
    subroutine GetStress3D_ut(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(3,3), tensor2(6,3)
      
      tensor2(1,1)=tensor1(1,1)
      tensor2(1,2)=0
      tensor2(1,3)=0
      
      tensor2(2,1)=0
      tensor2(2,2)=tensor1(2,2)
      tensor2(2,3)=0
      
      tensor2(3,1)=0
      tensor2(3,2)=0 
      tensor2(3,3)=tensor1(3,3) 
      
      tensor2(4,1)=0
      tensor2(4,2)=tensor1(3,3)  
      tensor2(4,3)=tensor1(2,2)
      
      tensor2(5,1)=tensor1(3,3)
      tensor2(5,2)=0 
      tensor2(5,3)=tensor1(1,1)
      
      tensor2(6,1)=tensor1(2,2)
      tensor2(6,2)=tensor1(1,1)  
      tensor2(6,3)=0
   end subroutine  
   
   
   subroutine AddSigmaGradVel2D_ut(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(3,2), tensor2(3,2)
      
      tensor2(1,1)= tensor1(1,1) + tensor1(3,2) + tensor2(1,1)
      tensor2(1,2)= tensor2(1,2)
      tensor2(2,1)= tensor2(2,1)
      tensor2(2,2)= tensor1(3,1) + tensor1(2,2) + tensor2(2,2)
      tensor2(3,1)= tensor1(2,2) + 2.0_rp*tensor2(3,1)
      tensor2(3,2)= tensor1(1,1) + 2.0_rp*tensor2(3,2) 
   end subroutine  
   
   subroutine AddSigmaGradVel3D_ut(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(6,3), tensor2(6,3)
      
      tensor2(1,1)=(tensor1(1,1)+tensor1(6,2)+tensor1(5,3)) + tensor2(1,1)
      tensor2(1,2)=tensor2(1,2)
      tensor2(1,3)=tensor2(1,3)
      
      tensor2(2,1)=tensor2(2,1)
      tensor2(2,2)=(tensor1(6,1)+tensor1(2,2)+tensor1(4,3)) + tensor2(2,2)
      tensor2(2,3)=tensor2(2,3)
      
      tensor2(3,1)=tensor2(3,1)  
      tensor2(3,2)=tensor2(3,2)
      tensor2(3,3)=(tensor1(5,1)+tensor1(4,2)+tensor1(3,3)) + tensor2(3,3)
      
      tensor2(4,1)=2.0_rp*tensor2(4,1)
      tensor2(4,2)=(tensor1(5,1) + tensor1(4,2) + tensor1(3,3) + 2.0_rp*tensor2(4,2))
      tensor2(4,3)=(tensor1(6,1) + tensor1(2,2) + tensor1(4,3) + 2.0_rp*tensor2(4,3))
      
      tensor2(5,1)=(tensor1(5,1) + tensor1(4,2) + tensor1(3,3) + 2.0_rp*tensor2(5,1))
      tensor2(5,2)=2.0_rp*tensor2(5,2)
      tensor2(5,3)=(tensor1(1,1) + tensor1(6,2) + tensor1(5,3) + 2.0_rp*tensor2(5,3))
      
      tensor2(6,1)=(tensor1(6,1) + tensor1(2,2) + tensor1(4,3) + 2.0_rp*tensor2(6,1))
      tensor2(6,2)=(tensor1(1,1) + tensor1(6,2) + tensor1(5,3) + 2.0_rp*tensor2(6,2)) 
      tensor2(6,3)=2.0_rp*tensor2(6,3)
   end subroutine  
   
   
   subroutine GetGrVel2D_ut(vector,tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(2,2), tensor2(3,2),vector(2)
      
      tensor2(1,1)=vector(1)*tensor1(1,1)
      tensor2(1,2)=vector(1)*tensor1(1,2)
      tensor2(2,1)=vector(2)*tensor1(2,1)
      tensor2(2,2)=vector(2)*tensor1(2,2)
      tensor2(3,1)=vector(2)*tensor1(1,1)+vector(1)*tensor1(2,1)
      tensor2(3,2)=vector(2)*tensor1(1,2)+vector(1)*tensor1(2,2)   
   end subroutine  
   
   
   subroutine GetGrVel3D_ut(vector,tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(3,3), tensor2(6,3),vector(3)
      
      tensor2(1,1)=vector(1)*tensor1(1,1)
      tensor2(1,2)=vector(1)*tensor1(1,2)
      tensor2(1,3)=vector(1)*tensor1(1,3)
      tensor2(2,1)=vector(2)*tensor1(2,1)
      tensor2(2,2)=vector(2)*tensor1(2,2)
      tensor2(2,3)=vector(2)*tensor1(2,3)
      tensor2(3,1)=vector(3)*tensor1(3,1)
      tensor2(3,2)=vector(3)*tensor1(3,2)
      tensor2(3,3)=vector(3)*tensor1(3,3)
      
      tensor2(4,1)=vector(3)*tensor1(2,1)+vector(2)*tensor1(3,1)
      tensor2(4,2)=vector(3)*tensor1(2,2)+vector(2)*tensor1(3,2)
      tensor2(4,3)=vector(3)*tensor1(2,3)+vector(2)*tensor1(3,3)
      
          
      tensor2(5,1)=vector(3)*tensor1(1,1)+vector(1)*tensor1(3,1)
      tensor2(5,2)=vector(3)*tensor1(1,2)+vector(1)*tensor1(3,2)
      tensor2(5,3)=vector(3)*tensor1(1,3)+vector(1)*tensor1(3,3)
      
      tensor2(6,1)=vector(2)*tensor1(1,1)+vector(1)*tensor1(2,1)
      tensor2(6,2)=vector(2)*tensor1(1,2)+vector(1)*tensor1(2,2)
      tensor2(6,3)=vector(2)*tensor1(1,3)+vector(1)*tensor1(2,3)
     
   end subroutine  
   
   
   subroutine GetStress_ut(ndime,tensor1,tensor2)
      integer(ip) :: ndime
      real(rp)    :: tensor1(ndime,ndime), tensor2(:,:)
      
      if (ndime == 2) call GetStress2D_ut(tensor1,tensor2)
      if (ndime == 3) call GetStress3D_ut(tensor1,tensor2)
   end subroutine
   
   subroutine AddSigmaGradVel_ut(ndime,tensor1,tensor2)
      integer(ip) :: ndime
      real(rp)    :: tensor1(ndime,ndime), tensor2(:,:)
      
      if (ndime == 2) call AddSigmaGradVel2D_ut(tensor1,tensor2)
      if (ndime == 3) call AddSigmaGradVel3D_ut(tensor1,tensor2)
   end subroutine
   
      subroutine GetGrVel_ut(ndime,vector,tensor1,tensor2)
      integer(ip) :: ndime
      real(rp)    :: vector(ndime,ndime), tensor1(ndime,ndime), tensor2(:,:)
      
      if (ndime == 2) call GetGrVel2D_ut(vector,tensor1,tensor2)
      if (ndime == 3) call GetGrVel3D_ut(vector,tensor1,tensor2)
   end subroutine
   
   
!------------------------------ 
! two vectors uv
!------------------------------   

   subroutine GetVel2D_uv(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(2,2), tensor2(2,2)
       
      tensor2(1,1) = tensor1(2,2)
      tensor2(1,2) = tensor1(2,1)
      tensor2(2,1) = tensor1(1,2)
      tensor2(2,2) = tensor1(1,1)
   end subroutine
   
     subroutine AddSphericalComponent2D_uv(value,tensor)
      implicit none
      real(rp) :: value(2)
      real(rp) :: tensor(2,2)

      tensor(1,1) = tensor(1,1) + value(1)
      tensor(2,2) = tensor(2,2) + value(2)
      
   end subroutine
   
    subroutine GetVel3D_uv(tensor1,tensor2)
      implicit none
      real(rp) :: tensor1(3,3), tensor2(3,3)
       
      tensor2(1,1) = tensor1(2,2) + tensor1(3,3)
      tensor2(1,2) = tensor1(2,1)
      tensor2(1,3) = tensor1(3,1)
      tensor2(2,1) = tensor1(1,2)
      tensor2(2,2) = tensor1(1,1) + tensor1(3,3)
      tensor2(2,3) = tensor1(3,2)
      tensor2(3,1) = tensor1(1,3)
      tensor2(3,2) = tensor1(2,3)
      tensor2(3,3) = tensor1(1,1) + tensor1(2,2)
   end subroutine
   
   subroutine AddSphericalComponent3D_uv(value,tensor)
      implicit none
      real(rp) :: value(3)
      real(rp) :: tensor(3,3)

      tensor(1,1) = tensor(1,1) + value(1)
      tensor(2,2) = tensor(2,2) + value(2)
      tensor(3,3) = tensor(3,3) + value(3)
      
   end subroutine
   
   
    subroutine GetVel_uv(ndime,tensor1,tensor2)
      integer(ip) :: ndime
      real(rp)    :: tensor1(ndime,ndime), tensor2(ndime,ndime)
      
      if (ndime == 2) call GetVel2D_uv(tensor1,tensor2)
      if (ndime == 3) call GetVel3D_uv(tensor1,tensor2)
   end subroutine
   
    subroutine AddSphericalComponent_uv(ndime,value,tensor1)
      integer(ip) :: ndime
      real(rp)    :: tensor1(ndime,ndime), value(ndime)
      
      if (ndime == 2) call AddSphericalComponent2D_uv(value,tensor1)
      if (ndime == 3) call AddSphericalComponent3D_uv(value,tensor1)
   end subroutine
   
   
!----------------------------
! rhu
! ----------------------------

   subroutine GetVel2D_hu(tensor,vector)
      implicit none
      real(rp) :: tensor(2,3), vector(2)
       
      vector(1) = tensor(1,1)+tensor(2,3)
      vector(2) = tensor(2,2)+tensor(1,3)
   end subroutine
   
   subroutine GetVel3D_hu(tensor,vector)
      implicit none
      real(rp) :: tensor(3,6), vector(3)
       
      vector(1) = tensor(1,1)+tensor(2,6)+tensor(3,5)
      vector(2) = tensor(1,6)+tensor(2,2)+tensor(3,4)
      vector(3) = tensor(1,5)+tensor(2,4)+tensor(3,3)
   end subroutine
   
    subroutine GetVel_hu(ndime,tensor,vector)
      integer(ip) :: ndime
      real(rp)    :: tensor(:,:), vector(ndime)
      
      if (ndime == 2) call GetVel2D_hu(tensor,vector)
      if (ndime == 3) call GetVel3D_hu(tensor,vector)
   end subroutine
   
!----------------------------
! rhc
! ----------------------------

   subroutine GetGrSig2D_hc(tensor,vector)
      implicit none
      real(rp) :: tensor(2,3), vector(3)
       
      vector(1) = tensor(1,1)+tensor(2,1)
      vector(2) = tensor(1,2)+tensor(2,2)
      vector(3) = 2.0_rp*(tensor(1,3)+tensor(2,3))
   end subroutine
   
   subroutine GetGrVel2D_hc(tensor,vector1,vector2)
      implicit none
      real(rp) :: tensor(2,2), vector1(3), vector2(3)
       
      vector2(1) = tensor(1,1)*vector1(1)+tensor(1,2)*vector1(3)
      vector2(2) = tensor(2,1)*vector1(3)+tensor(2,2)*vector1(2)
      vector2(3) = tensor(2,1)*vector1(1)+tensor(1,2)*vector1(2)
      !vector2(3) = tensor(2,1)*vector1(1)+tensor(1,2)*vector1(2)+(tensor(1,1)+tensor(2,2))*vector1(3)
   end subroutine
   
   
    subroutine GetStress2D_hc(auxPTT,tensor,vector)
      implicit none
      real(rp) :: tensor(3,3), vector(3)
      integer :: auxPTT
      vector(1) = auxPTT*(tensor(1,1)+tensor(3,3))+(1 - auxPTT)*(tensor(1,1)+tensor(1,2))
      vector(2) = auxPTT*(tensor(3,3)+tensor(2,2))+(1 - auxPTT)*(tensor(1,2)+tensor(2,2))
      vector(3) = 2.0_rp*(tensor(1,3)+tensor(2,3))
   end subroutine
   
    subroutine AddValues2D_hc(value,vector)
      implicit none
      real(rp) :: value(3), vector(3)
      vector(1) = vector(1) + value(1)
      vector(2) = vector(2) + value(2)
      vector(3) = vector(3) + 2.0_rp*value(3)
   end subroutine
   
   
   subroutine GetGrSig3D_hc(tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(3,6), vector(6)
       
      vector(1) = tensor(1,1)+tensor(2,1)+tensor(3,1)
      vector(2) = tensor(1,2)+tensor(2,2)+tensor(3,2)
      vector(3) = tensor(1,3)+tensor(2,3)+tensor(3,3)
      vector(4) = 2.0_rp*(tensor(1,4)+tensor(2,4)+tensor(3,4))
      vector(5) = 2.0_rp*(tensor(1,5)+tensor(2,5)+tensor(3,5))
      vector(6) = 2.0_rp*(tensor(1,6)+tensor(2,6)+tensor(3,6))
   end subroutine
   
   subroutine GetGrVel3D_hc(tensor,vector1,vector2)
      use typre
      implicit none
      real(rp) :: tensor(3,3), vector1(6), vector2(6)
       
      vector2(1) = tensor(1,1)*vector1(1)+tensor(1,2)*vector1(6)+tensor(1,3)*vector1(5)
      vector2(2) = tensor(2,1)*vector1(6)+tensor(2,2)*vector1(2)+tensor(2,3)*vector1(4)
      vector2(3) = tensor(3,1)*vector1(5)+tensor(3,2)*vector1(4)+tensor(3,3)*vector1(3)
      vector2(4) = tensor(3,1)*vector1(6)+tensor(3,2)*vector1(2)+(tensor(3,3)+tensor(2,2))*vector1(4)+vector1(5)*tensor(2,1)+vector1(3)*tensor(2,3)
      vector2(5) = tensor(3,1)*vector1(1)+tensor(3,2)*vector1(6)+(tensor(3,3)+tensor(1,1))*vector1(5)+vector1(4)*tensor(1,2)+vector1(3)*tensor(1,3)
      vector2(6) = tensor(2,1)*vector1(1)+tensor(1,2)*vector1(2)+(tensor(2,2)+tensor(1,1))*vector1(6)+vector1(5)*tensor(2,3)+vector1(4)*tensor(1,3)
      
   end subroutine
   
    subroutine GetStress3D_hc(auxPTT,tensor,vector)
      use typre
      implicit none
      real(rp) :: tensor(6,6), vector(6)
      integer :: auxPTT
      vector(1) = tensor(1,1)+auxPTT*(tensor(6,6)+tensor(5,5))
      vector(2) = tensor(2,2)+auxPTT*(tensor(6,6)+tensor(4,4))
      vector(3) = tensor(3,3)+auxPTT*(tensor(5,5)+tensor(4,4))
      vector(4) = (1 + auxPTT)*(tensor(2,4)+tensor(3,4))+2.0_rp*auxPTT*tensor(6,5)
      vector(5) = (1 + auxPTT)*(tensor(1,5)+tensor(3,5))+2.0_rp*auxPTT*tensor(4,6)
      vector(6) = (1 + auxPTT)*(tensor(1,6)+tensor(2,6))+2.0_rp*auxPTT*tensor(4,5)
   end subroutine
   
   subroutine AddValues3D_hc(value,vector)
      use typre
      implicit none
      real(rp) :: value(6), vector(6)
      vector(1:3) = vector(1:3) + value(1:3)
      vector(4:6) = vector(4:6) + 2.0_rp*value(4:6)
   end subroutine
   
   subroutine GetGrSig_hc(ndime,tensor,vector)
      integer(ip) :: ndime
      real(rp)    :: tensor(:,:), vector(:)
      
      if (ndime == 2) call GetGrSig2D_hc(tensor,vector)
      if (ndime == 3) call GetGrSig3D_hc(tensor,vector)
   end subroutine

    subroutine GetGrVel_hc(ndime,tensor,vector1,vector2)
      use typre
      implicit none
      integer(ip) :: ndime
      real(rp) :: tensor(:,:)
      real(rp) :: vector1(:), vector2(:)
      if (ndime == 2) call GetGrVel2D_hc(tensor,vector1,vector2)
      if (ndime == 3) call GetGrVel3D_hc(tensor,vector1, vector2)
   end subroutine
   
    subroutine GetStress_hc(ndime,auxPTT,tensor,vector)
      use typre
      implicit none
      integer(ip) ::   ndime
      real(rp) :: tensor(:,:), vector(:)
      integer :: auxPTT
      if (ndime == 2) call GetStress2D_hc(auxPTT,tensor,vector)
      if (ndime == 3) call GetStress3D_hc(auxPTT,tensor,vector)
   end subroutine
   
    subroutine AddValues_hc(ndime,value,vector)
      use typre
      implicit none
      integer(ip) ::  auxPTT, ndime
      real(rp) :: value(:), vector(:)
      if (ndime == 2) call AddValues2D_hc(value,vector)
      if (ndime == 3) call AddValues3D_hc(value,vector)
   end subroutine
   
   subroutine GetGrVel_hu(ndime,tensor,vector)
      use typre
      integer(ip) ::  ndime, itens
      real(rp)    ::  vector(ndime), tensor(ndime,ndime)
       vector=0.0_rp
       do itens=1,ndime
            vector(:)=tensor(:,itens)+vector(:)
        end do
   end subroutine   
   
    subroutine GetGrVel(ndime,vector1,vector2)
      use typre
      integer(ip) ::  ndime, itens
      real(rp)    ::  vector1(ndime), vector2(ndime)
       !vector2=0.0_rp
       do itens=1,ndime
            vector2(itens)=vector1(itens)+vector2(itens)
        end do
   end subroutine  
   
   
   subroutine GetPlus(ndime,vector,value)
      use typre
      integer(ip) ::  ndime, itens
      real(rp)    ::  vector(ndime), value
       value=0
       value=sum(vector(1:ndime))
   end subroutine 
   
   
   
!***************************************************************************
! LCR RELATIONS
!***************************************************************************
   
   subroutine PassTensor4ToTensor2(e,tens,Tensor4,Tensor2)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: h,s,t, tens
      real(rp) :: Tensor4(e%ndime,e%ndime,e%ndime,e%ndime)
      real(rp) :: Tensor2(tens,tens)
      real(rp) :: vector1(tens), vector2(tens)
      real(rp) :: tens_stress(e%ndime,e%ndime)
      real(rp) :: tensor_aux2(e%ndime,e%ndime,(e%ndime-1)*(e%ndime-1)+2)
      
      
      do s=1,e%ndime
         do t=1,e%ndime
            tens_stress(:,:)=Tensor4(s,t,:,:)
            call GetStress(e%ndime,tens_stress,vector1)
            tensor_aux2(s,t,:)=vector1(:) 
         end do
      end do 

       do h=1,tens
         tens_stress(:,:)=tensor_aux2(:,:,h)
         call GetStress(e%ndime,tens_stress,vector2)
         Tensor2(:,h)=vector2(:)
       end do
       
   end subroutine
   
   subroutine PassTensor3ToTensor2Left(e,tens,Tensor3,Tensor2)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: h,s,t, tens
      real(rp) :: Tensor3(e%ndime,e%ndime,e%ndime)
      real(rp) :: vector(tens)
      real(rp) :: Tensor2(tens,e%ndime)
      real(rp) :: tens_stress(e%ndime,e%ndime)
      

      do h=1,e%ndime
         tens_stress(:,:)=Tensor3(:,:,h)
         call GetStress(e%ndime,tens_stress,vector)
         Tensor2(:,h)=vector(:)
      end do
         
   end subroutine
   
   
   
   subroutine PassTensor3ToTensor2Left_sym(e,tens,Tensor3,Tensor2)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: h,s,t, tens
      real(rp) :: Tensor3(e%ndime,e%ndime,e%ndime)
      real(rp) :: vector(tens)
      real(rp) :: Tensor2(tens,e%ndime)
      real(rp) :: tens_stress(e%ndime,e%ndime)
      

      do h=1,e%ndime
         tens_stress(:,:)=Tensor3(:,:,h)
!          call GetStress_sym(e%ndime,tens_stress,vector)
         call PassSymMatrixToVector(e%ndime,tens,tens_stress,vector)
         Tensor2(:,h)=vector(:)
      end do
         
   end subroutine
   
    subroutine PassTensor3ToTensor2Right(e,tens,Tensor3,Tensor2)
      implicit none
      class(FiniteElement) :: e
      integer(ip) :: h, tens
      real(rp) :: Tensor3(e%ndime,e%ndime,e%ndime)
      real(rp) :: Tensor2(e%ndime,tens)
      real(rp) :: vector(tens)
      real(rp) :: tens_stress(e%ndime,e%ndime)
      
  
      do h=1,e%ndime
         tens_stress(:,:)=Tensor3(h,:,:)
         call GetStress(e%ndime,tens_stress,vector)
         Tensor2(h,:)=vector(:)
      end do

   end subroutine
   
   
   subroutine PassGradStressToTensor3(e,tens,Grsig,Tensor3)
      implicit none
      class(FiniteElement) :: e
      integer(ip), intent(in) ::tens
      real(rp), intent(in) :: Grsig(tens,e%ndime)
      real(rp) , intent(out) :: Tensor3(e%ndime,e%ndime,e%ndime)
      real(rp) :: tens_stress(e%ndime,e%ndime),vector(tens)
      integer(ip) :: h,i,j
      
      
      do h=1,e%ndime
         vector(:)=Grsig(:,h)
         call sup_SigmaMatrix(e%ndime,tens,vector,tens_stress)
         do concurrent (i=1:e%ndime)
            Tensor3(:,i,h)=tens_stress(:,i)
         end do   
      end do

   end subroutine
   
   
   subroutine PassSymMatrixToVector(ndime,tens,Matrix,Vector)
      implicit none
      integer(ip), intent(in) ::tens,ndime
      real(rp), intent(in) :: Matrix(ndime,ndime)
      real(rp) , intent(out) :: Vector(tens)
      integer(ip) :: h,i,j
      
      
      do i=1,ndime
         Vector(i)=Matrix(i,i)
      end do 
      
      if (ndime==2) then
         Vector(3)=Matrix(1,2)
      else if (ndime==3) then 
         Vector(4)=Matrix(2,3)
         Vector(5)=Matrix(1,3)
         Vector(6)=Matrix(1,2)
      end if
 

   end subroutine
   
   subroutine PassGradStressToTensor3B(ndime,tens,Grsig,Tensor3)
      implicit none
      integer(ip), intent(in) ::tens,ndime
      real(rp), intent(in) :: Grsig(tens,ndime)
      real(rp) , intent(out) :: Tensor3(ndime,ndime,ndime)
      real(rp) :: tens_stress(ndime,ndime),vector(tens)
      integer(ip) :: h,i,j
      
      
      do h=1,ndime
         vector(:)=Grsig(:,h)
         call sup_SigmaMatrix(ndime,tens,vector,tens_stress)
         do concurrent (i=1:ndime)
               Tensor3(:,i,h)=tens_stress(:,i) 
         end do   
      end do

   end subroutine

end module
