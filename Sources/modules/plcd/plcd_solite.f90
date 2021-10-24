subroutine plcd_solite(a,itask)
   use typre
   use Mod_PLCD
   use Mod_plcd_Stages
   implicit none
   class(PLCDProblem), target :: a
   integer(ip) :: itask
   integer(ip) :: ndime
   real(rp) :: OptimalStepSize = 1.0_rp
   
   interface
      subroutine plcd_LineSearch(a,OptimalStepSize)
         import
         implicit none
         class(PLCDProblem), target :: a
         real(rp) :: OptimalStepSize
      end subroutine
   end interface
   
   if (itask == 1) then 
     
   
   elseif (itask == 2) then
      
      if(a%kfl_LineSearch == 1) call plcd_LineSearch(a,OptimalStepSize)
      
      call a%Mesh%GetNdime(ndime)
      a%unkno(1:ndime,:) = OptimalStepSize*a%unkno(1:ndime,:) + a%Displacement(:,:,1)
      
      !UPFormulation
      if (a%UseUPFormulation) a%unkno(ndime+1,:) = OptimalStepSize*a%unkno(ndime+1,:) + a%Pressure(:,1)
      
   endif
  
   
end subroutine plcd_solite
