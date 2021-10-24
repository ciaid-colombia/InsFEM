subroutine nsc_CalculateStatistics(a)
   use typre
   use Mod_NSCompressible

   implicit none
   class(NSCompressibleProblem) :: a
   
   integer(ip)             :: idime,ipoin,ndime,npoin
   integer(ip)             :: calculate_statistics

   calculate_statistics=0

   if(a%npp_stepi(15)>0) then     
          calculate_statistics=1
   endif 

   if (calculate_statistics == 1) then 
       call a%Mesh%GetNpoin(npoin) 
       call a%Mesh%GetNdime(ndime)
       do ipoin = 1,npoin
          
         a%timad(ipoin,2) = a%timad(ipoin,2) + a%densf(ipoin,1)*a%dtime
         a%timad(ipoin,1) = a%timad(ipoin,2)/a%ctime
         a%rmsd(ipoin) = a%rmsd(ipoin) +&
                        (a%densf(ipoin,1)-a%timad(ipoin,1))*(a%densf(ipoin,1)-a%timad(ipoin,1))*a%dtime
         do idime = 1,ndime
            a%timav(idime,ipoin,2) =  a%timav(idime,ipoin,2) + a%veloc(idime,ipoin,1)*a%dtime   
            a%timav(idime,ipoin,1) =  a%timav(idime,ipoin,2)/a%ctime   
            a%rmsv(idime,ipoin) = a%rmsv(idime,ipoin) + &
                             (a%veloc(idime,ipoin,1)-a%timav(idime,ipoin,1))*(a%veloc(idime,ipoin,1)-a%timav(idime,ipoin,1))*a%dtime
         end do
         a%timap(ipoin,2) = a%timap(ipoin,2) + a%press(ipoin,1)*a%dtime
         a%timap(ipoin,1) = a%timap(ipoin,2)/a%ctime
         a%rmsp(ipoin) = a%rmsp(ipoin) +&
                        (a%press(ipoin,1)-a%timap(ipoin,1))*(a%press(ipoin,1)-a%timap(ipoin,1))*a%dtime
         a%timat(ipoin,2) = a%timat(ipoin,2) + a%tempe(ipoin,1)*a%dtime
         a%timat(ipoin,1) = a%timat(ipoin,2)/a%ctime
         a%rmst(ipoin) = a%rmst(ipoin) +&
                         (a%tempe(ipoin,1)-a%timat(ipoin,1))*(a%tempe(ipoin,1)-a%timat(ipoin,1))*a%dtime
    
       end do
   end if

end subroutine nsc_CalculateStatistics

subroutine nsc_FinalizeStats(a)
   use typre
   use Mod_NSCompressible

   implicit none
   class(NSCompressibleProblem) :: a
   
   integer(ip)             :: idime,ipoin,ndime,npoin
   real(rp)                :: rms,velnorm

   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNdime(ndime)
   

   do ipoin = 1,npoin
      a%rmsd(ipoin) = a%rmsd(ipoin)/a%ctime
      a%rmsd(ipoin) = sqrt(a%rmsd(ipoin))
      a%rmsp(ipoin) = a%rmsp(ipoin)/a%ctime
      a%rmsp(ipoin) = sqrt(a%rmst(ipoin))
      a%rmst(ipoin) = a%rmst(ipoin)/a%ctime
      a%rmst(ipoin) = sqrt(a%rmst(ipoin))
      rms = 0.0_rp
      velnorm = 0.0_rp
      do idime=1,ndime
         a%rmsv(idime,ipoin) = a%rmsv(idime,ipoin)/a%ctime
         rms = rms + a%rmsv(idime,ipoin)*a%rmsv(idime,ipoin)
         a%rmsv(idime,ipoin) = sqrt(a%rmsv(idime,ipoin))
         velnorm = velnorm + a%timav(idime,ipoin,1)*a%timav(idime,ipoin,1)
      end do
      rms = rms/ndime
      rms = sqrt(rms)
      velnorm = sqrt(velnorm) 
      if(velnorm>epsilon(0.0_rp)) then
         a%turbi(ipoin)=rms/velnorm
      else
         a%turbi(ipoin)=0.0_rp
      end if
   end do
end subroutine nsc_FinalizeStats
