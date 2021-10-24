subroutine plcd_updbcs(a)
   ! DESCRIPTION
   !    This routine prepares for a new time step of the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_PLCD
   use Mod_Mesh   
   use Mod_TimeIntegrator   
   use Mod_plcd_Stages
   use Mod_plcdExacso
   use Mod_plcd_TD_StochasticTopologyOptimization
   implicit none
   class(PLCDProblem), target :: a
   
   integer(ip) :: ndime,ipoin,npoin
   real(rp) :: exvel(3), exveg(3,3)
   type(plcdExacso) :: Exacso
   real(rp), pointer :: coord(:)
   
   
   !If analytical solution, overwrite whatever we have in bvess
   if (a%kfl_exacs > 0) then
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)
      do ipoin = 1,npoin
         call a%Mesh%GetPointCoord(ipoin,coord)
         
         call Exacso%ComputeSolution(ndime,coord,a)
         call Exacso%GetDisplacement(ndime,exvel,exveg)
         
         a%cs%bvess(1:ndime,ipoin) = exvel(1:ndime)
      enddo
   endif
   
  
   
   
end subroutine plcd_updbcs

