subroutine nsc_CalculatePrimitivesLogic(a)
   use typre
   use Mod_NSCompressible

   implicit none

   class(NSCompressibleProblem) :: a
   
   logical :: isALE
   integer(ip)             :: calculate_primitives

   calculate_primitives=0

   call a%Mesh%GetALE(isALE)
   if (isALE) calculate_primitives=1

   if(a%nptra > 0) calculate_primitives=1

   if (a%npp_stepi(1)>0) then
      if (mod(a%istep,a%npp_stepi(1))==0) calculate_primitives=1
   endif 
   if (a%npp_stepi(2)>0) then
      if (mod(a%istep,a%npp_stepi(2))==0) calculate_primitives=1
   endif 
   if (a%npp_stepi(3)>0) then
      if (mod(a%istep,a%npp_stepi(3))==0) calculate_primitives=1
   endif 
   if (a%npp_stepi(7)>0) then
      if (mod(a%istep,a%npp_stepi(7))==0) calculate_primitives=1
   endif 
   if (a%npp_stepi(15)>0) then
          calculate_primitives=1
   endif 

   if (calculate_primitives == 1) call a%CalculatePrimitives(1)

end subroutine nsc_CalculatePrimitivesLogic

subroutine nsc_CalculatePrimitives(a,ncomp)
   use typre
   use Mod_NSCompressible

   implicit none

   class(NSCompressibleProblem) :: a
   integer(ip),intent(in)   :: ncomp
   
   integer(ip)             :: iffix_vel,iffix_tem
   integer(ip)             :: idime,ipoin,ndime,npoin,icomp
   real(rp)                :: acvis,actco,accph,accvh
   real(rp)                :: invden, auxiliar


   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)

   do icomp = 1,ncomp

      do ipoin = 1,npoin
      
         do idime = 1,ndime
            a%veloc(idime,ipoin,icomp) = a%momen(idime,ipoin,icomp) / a%densf(ipoin,icomp)
            iffix_vel = a%kfl_fixno(idime+1,ipoin)
            if(iffix_vel==1) then
               a%veloc(idime,ipoin,icomp) = a%bvess(idime+1,ipoin,1)
            end if
         end do
   
         ! Ideal Gas State Law (could expand into a module)
         invden = 1.0_rp / a%densf(ipoin,icomp)  
   
         auxiliar = (a%energ(ipoin,icomp) - dot_product(a%momen(:,ipoin,icomp),a%momen(:,ipoin,icomp)) * invden / 2.0_rp )/ accvh
      
         a%tempe(ipoin,icomp) = (invden * auxiliar) - a%reltem
         a%press(ipoin,icomp) = ((accph - accvh) * auxiliar) - a%relpre
      
         iffix_tem = a%kfl_fixno(a%ndofbc,ipoin)
         if(iffix_tem==1) then
             a%tempe(ipoin,icomp)= a%bvess(a%ndofbc,ipoin,1)
         end if
      
      enddo
   enddo

end subroutine nsc_CalculatePrimitives
