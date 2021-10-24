subroutine plcd_PoinOpe(b)
   use Mod_plcd_BaseElmope
   use Mod_PLCD
   use Mod_plcd_BMatrixFactory
   use Mod_plcd_Stages
   use Mod_plcd_elmdir
   use Mod_php_AssemblyVectorToSystem
   implicit none
   class(PLCDProblem), target :: b
   
   integer(ip) :: idofn, inode,jnode
   
   a => b
   
   !First of all we assembly to the RHS the vector of the residual forces
   call php_AssemblyVectorToSystem(a,a%ResidualForcesVector)
   
   
   

end subroutine plcd_PoinOpe
