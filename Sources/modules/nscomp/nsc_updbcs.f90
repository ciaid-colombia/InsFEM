subroutine nsc_updbcs(a)
!-----------------------------------------------------------------------
!****f* Nscomp/nsc_updbcs
! NAME 
!    nsc_updbcs
! DESCRIPTION
!    This routine updates the velocity boundary conditions:
!    1. Before a time step begins
!    2. Before a global iteration begins
!    3. Before an inner iteration begins
! USED BY
!    nsc_begste
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_Element
   use Mod_memor   
   implicit none
   class(NSCompressibleProblem) :: a
   
   integer(ip)        :: idime,ipoin,ndime,npoin
   logical            :: isALE  = .FALSE.
   logical            :: funbcs = .FALSE.
   real(rp), pointer  :: meshve(:,:,:) => NULL()


   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNdime(ndime)
   
   call a%Mesh%GetALE(isALE)
 
   if (isALE) then
      call a%Mesh%GetMeshVeloc(meshve)
      do ipoin=1,npoin
         do idime=1,ndime
            if(a%kfl_fixno(idime+1,ipoin)==1) then
               a%bvess(idime+1,ipoin,1)=a%bvess(idime+1,ipoin,2) + meshve(idime,ipoin,1)
            end if
         end do
      end do 
   end if
 
   if(a%kfl_conbc /=1) then
        if(maxval(a%kfl_funty(:,1))==9) funbcs = .TRUE.
   end if
 
   if (funbcs) then
      do ipoin = 1,npoin
         if(a%kfl_fixno(1,ipoin)==1) then
            if(a%kfl_funno(ipoin)>0) then
               if (a%kfl_funty(1,a%kfl_funno(ipoin)) == 9) then
                  a%bvess(1,ipoin,1)=a%bvess(1,ipoin,2)
               endif
            end if
         end if
         if(a%kfl_fixno(a%ndofbc,ipoin)==1) then
            if(a%kfl_funno(ipoin)>0) then
               if (a%kfl_funty(1,a%kfl_funno(ipoin)) == 9) then
                  a%bvess(a%ndofbc,ipoin,1)=a%bvess(a%ndofbc,ipoin,2)
               endif
            end if
         end if
      end do 
   end if
       
end subroutine nsc_updbcs
