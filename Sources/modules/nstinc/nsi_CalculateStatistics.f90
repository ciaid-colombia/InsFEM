module Mod_nsi_Statistics
contains
subroutine nsi_CalculateStatistics(a)

    !-----------------------------------------------------------------------
    !
    ! This routine computes the Statistics: Time Average Velocity, Time average Pressure,  
    ! RMS Velocity, RMS Pressure and Turbulence Intesity
    !
    !-----------------------------------------------------------------------

   use typre
   use Mod_NavierStokes
   use Mod_Mesh

   implicit none
   class(NavierStokesProblem) :: a
   
   integer(ip)                :: idime,ipoin,ndime,npoin

   !Only if they need to be computed      
   if (a%npp_stepi(24) == 0) return   
   if (a%ctime < a%StatisticsStartTime) return
           
        

   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNdime(ndime)
   a%StatisticsCumulativeTime = a%StatisticsCumulativeTime + a%dtime
   
   a%timav(:,:,2) = a%timav(:,:,2) + a%veloc(:,:,1)*a%dtime
   a%timav(:,:,1) = a%timav(:,:,2)/a%StatisticsCumulativeTime
   a%timav(:,:,4) = a%timav(:,:,4) + abs(a%veloc(:,:,1))*a%dtime
   a%timav(:,:,5) = a%timav(:,:,4)/a%StatisticsCumulativeTime
   a%timav(:,:,3) = a%timav(:,:,3) + a%veloc(:,:,1)*a%veloc(:,:,1)*a%dtime
   
   a%timap(:,2) = a%timap(:,2) + a%press(:,1)*a%dtime
   a%timap(:,1) = a%timap(:,2)/a%StatisticsCumulativeTime
   a%timap(:,3) = a%timap(:,3) + a%press(:,1)*a%press(:,1)*a%dtime

   a%TimaTractions(:,:,2) = a%TimaTractions(:,:,2) + a%btraction(:,:)*a%dtime
   a%TimaTractions(:,:,1) = a%TimaTractions(:,:,2)/a%StatisticsCumulativeTime
   
   
   
end subroutine nsi_CalculateStatistics   

subroutine nsi_FinalStatistics(a)
   !Calculate Root Mean Square of Velocity (rmsv) and Pressure (rmsp)
   use typre
   use Mod_NavierStokes

   implicit none
   class(NavierStokesProblem) :: a
   
   integer(ip)             :: idime,ipoin,ndime,npoin
   real(rp)                :: rmsnorm,velnorm
   
  
!    !Only if they need to be computed 

   if (a%npp_stepi(24) == 0) return   
           
        

   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNdime(ndime)
   
   a%rmsp(:) = 0.0_rp
   a%rmsv(:,:) = 0.0_rp
   
   
 
    do ipoin = 1,npoin
       

       a%rmsp(ipoin) = a%rmsp(ipoin)+ a%timap(ipoin,3)+ (a%timap(ipoin,1)*a%timap(ipoin,1)*a%StatisticsCumulativeTime) - (2*a%timap(ipoin,1)*a%timap(ipoin,2))
       a%rmsp(ipoin) = a%rmsp(ipoin)/a%StatisticsCumulativeTime
       if (a%rmsp(ipoin) > 0.0_rp) then
         a%rmsp(ipoin) = sqrt(a%rmsp(ipoin))
       else
         a%rmsp(ipoin) = 0.0_rp
       endif
      
       rmsnorm = 0.0_rp
       velnorm = 0.0_rp
       
       do idime=1,ndime
          

          a%rmsv(idime,ipoin) = a%rmsv(idime,ipoin)+a%timav(idime,ipoin,3)+a%timav(idime,ipoin,1)*a%timav(idime,ipoin,1)*a%StatisticsCumulativeTime-2*a%timav(idime,ipoin,1)*a%timav(idime,ipoin,2)
          a%rmsv(idime,ipoin) = a%rmsv(idime,ipoin)/a%StatisticsCumulativeTime

          if (a%rmsv(idime,ipoin) > 0.0_rp) then
            a%rmsv(idime,ipoin) = sqrt(a%rmsv(idime,ipoin))
          else
            a%rmsv(idime,ipoin) = 0.0_rp
          endif
          velnorm = velnorm + a%timav(idime,ipoin,5)*a%timav(idime,ipoin,5)
          rmsnorm = rmsnorm + a%rmsv(idime,ipoin)*a%rmsv(idime,ipoin)


       end do
 
       !Calculate the Turbulence Intensity (turbi = Urms/Umean)
       rmsnorm = rmsnorm/ndime
       rmsnorm = sqrt(rmsnorm)
       velnorm = sqrt(velnorm)
       
       if(velnorm>epsilon(0.0_rp)) then
          a%turbi(ipoin)=rmsnorm/velnorm
       else
          a%turbi(ipoin)=0.0_rp
       end if
       
    end do
 
 end subroutine
 
end module
