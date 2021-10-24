subroutine nsc_CalculateConservatives(a,ncomp)
   use typre
   use Mod_NSCompressible

   implicit none

   class(NSCompressibleProblem) :: a
   integer(ip),intent(in)   :: ncomp

   integer(ip)             :: iffix_vel,iffix_tem
   integer(ip)             :: idime,ipoin,ndime,npoin,icomp
   real(rp)                :: acvis,actco,accph,accvh


   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)
   
   do icomp = 1,ncomp

      do ipoin = 1,npoin
   
         iffix_tem = a%kfl_fixno(a%ndofbc,ipoin)
         if(iffix_tem==1) then
             a%tempe(ipoin,icomp)= a%bvess(a%ndofbc,ipoin,1)
         end if
      
         ! Ideal Gas State Law (could expand into a module)
         a%densf(ipoin,icomp) = (a%press(ipoin,icomp)+a%relpre)/((accph - accvh) * (a%tempe(ipoin,icomp)+a%reltem))
   
         do idime = 1,ndime
            iffix_vel = a%kfl_fixno(idime+1,ipoin)
            if(iffix_vel==1) then
               a%veloc(idime,ipoin,icomp) = a%bvess(idime+1,ipoin,1)
            end if
            a%momen(idime,ipoin,icomp) = a%veloc(idime,ipoin,icomp) * a%densf(ipoin,icomp)
         end do
   
   
         a%energ(ipoin,icomp) = a%densf(ipoin,icomp)*(accvh*(a%tempe(ipoin,icomp)+a%reltem)+dot_product(a%veloc(:,ipoin,icomp),a%veloc(:,ipoin,icomp))/ 2.0_rp)
      
      enddo
   enddo

end subroutine nsc_CalculateConservatives
