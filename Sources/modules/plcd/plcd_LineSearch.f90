subroutine plcd_LineSearch(a,OptimalStepSize)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem), target :: a
   real(rp):: OptimalStepSize
   
   real(rp) :: Residual(2)
   integer(ip) :: LineSearchiterations = 1
   real(rp) :: Chi(2)
   real(rp) :: Alpha_parameter = 0.0_rp
   interface
      subroutine plcd_ComputeResidual(a,Chi,Residual)
         import
         implicit none
         class(PLCDProblem), target :: a
         real(rp) :: Chi,Residual
      end subroutine
   end interface
   
   a%kfl_InternalLineSearchLoop = 1
   a%Displacement(:,:,2) = a%Displacement(:,:,1) 
   
   Chi(1) = 0.0_rp
   Chi(2) = 100.0_rp
   
   call plcd_ComputeResidual(a,Chi(1),Residual(1))
   
   if (Residual(1) .GT. 0.0_rp) then
      Residual(1) = - Residual(1)
      a%unkno = - a%unkno
   endif
   
   call plcd_ComputeResidual(a,Chi(2),Residual(2))
   
   LineSearchiterations = 1
   
   do while (abs(Residual(2)) > a%LineSearchRadiusOfConvergence*abs(Residual(1)) .AND. LineSearchiterations .LE. a%LineSearchMaxIterations)
  
      Alpha_parameter = Residual(1) / Residual(2)
      
      if (Alpha_parameter .LT. 0.0_rp) then
         Chi(2) = Alpha_parameter/2.0_rp + sqrt((Alpha_parameter/2.0_rp)*(Alpha_parameter/2.0_rp)-Alpha_parameter)
      else
         Chi(2) = Alpha_parameter/2.0_rp
      endif
      
      call plcd_ComputeResidual(a,Chi(2),Residual(2))
      
      LineSearchiterations = LineSearchiterations + 1

   enddo
   
   if (LineSearchiterations .GE. a%LineSearchMaxIterations) then
      OptimalStepSize = 1.0_rp
   else   
      OptimalStepSize = Chi(2)
   endif
   a%kfl_InternalLineSearchLoop = 0
   a%Displacement(:,:,1) = a%Displacement(:,:,2)
   a%Displacement(:,:,2) = 0.0_rp
    
 end subroutine

subroutine plcd_ComputeResidual(a,Chi,Residual)
   use typre
   use Mod_PLCD
   implicit none
   class(PLCDProblem), target :: a
   real(rp) :: Chi, Residual
   integer(ip) :: ndime, idime, npoin
   
   interface
      subroutine plcd_ComputeForcesVectors(a)
         use typre
         use Mod_PLCD
         implicit none
         class(PLCDProblem) :: a
      end subroutine
   end interface
   
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   
   a%Displacement(:,:,1) = a%Displacement(:,:,1) + Chi*a%unkno(1:ndime,:)
   call plcd_ComputeForcesVectors(a)
   Residual = 0.0_rp
   do idime = 1,ndime
      Residual = Residual + dot_product(a%unkno(idime,1:npoin),-a%ResidualForcesVector(idime,1:npoin))
   enddo
   a%Displacement(:,:,1) = a%Displacement(:,:,2)

end subroutine
   