subroutine nsc_elmdir_ex(a)
   !------------------------------------------------------------------------
   ! This routine modifies the element stiffness to impose the correct 
   ! boundary conditions
   !------------------------------------------------------------------------
   use typre
   use Mod_NSCompressibleExplicit
   use Mod_NSCompressibleSubroutines
   implicit none
   class(NSCompressibleExplicitProblem)           :: a
   
   integer(ip)                :: npoin,ipoin,idofn,ndime
   integer(ip)                :: iffix_den,iffix_tem,iffix_vel
   real(rp)                   :: energ,acvis,actco,accph,accvh
   
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNdime(ndime)
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)

   !Dirichlet boundary conditions
   do ipoin=1,npoin
      iffix_den = a%kfl_fixno(1,ipoin)
      if(iffix_den==1) then
         a%densf(ipoin,1)= a%bvess(1,ipoin,1)
      end if
      do idofn=1,a%ndofbc-2
         iffix_vel = a%kfl_fixno(idofn+1,ipoin)
         if(iffix_vel==1) then
            a%veloc(idofn,ipoin,1) = a%bvess(idofn+1,ipoin,1)
            a%momen(idofn,ipoin,1) = a%densf(ipoin,1) * a%bvess(idofn+1,ipoin,1)
         end if
      end do
      iffix_tem = a%kfl_fixno(a%ndofbc,ipoin)
      if(iffix_tem==1) then
         a%tempe(ipoin,1)= a%bvess(a%ndofbc,ipoin,1)
         call nsc_ComputeEnergy(ndime,accvh,a%tempe(ipoin,1),a%veloc(:,ipoin,1),energ)
         a%energ(ipoin,1) = a%densf(ipoin,1)*energ
      end if
   end do

end subroutine nsc_elmdir_ex
