subroutine sldale_begite(a)
   use typre
   use Mod_sldAlemov
   implicit none
   class(sldAlemovProblem), target :: a
   real(rp), pointer     :: r_i1(:,:) => NULL(),r_i21(:,:) => NULL()
   integer(ip) :: ipoin,npoin,idime,ndime
   real(rp) :: w_r

   !For solid interaction
   if (associated(a%sldDisp)) then
      call a%Mesh%GetNpoin(npoin)
      do ipoin=1,npoin
          do idime=1,a%ndofbc
              if(a%kfl_fixno(idime,ipoin)==4) then
                  !d_i1(:,:)  = d_i1(:,:)
                  a%bdisp(idime,ipoin,2) = a%bvess(idime,ipoin,1) 
                  !d_i2_guess(:,:)  = d_i2_guess(:,:)
                  a%bdisp(idime,ipoin,1) = a%sldDisp(idime,ipoin)
              end if
          end do
      end do  

      !r_i1(:,:)    = r_i2_past(:,:)
      a%bres(:,:,2) = a%bres(:,:,1)

      !r_i2(:,:)    = d_i2_guess(:,:)- d_i1(:,:)
      a%bres(:,:,1) = a%bdisp(:,:,1) - a%bdisp(:,:,2) 

      !r_i21(:,:)    = r_i2(:,:)    - r_i1(:,:)
      a%bres(:,:,3) = a%bres(:,:,1) - a%bres(:,:,2)

      if(a%kfl_doAitken) then
          r_i1   => a%bres(:,:,2)
          r_i21  => a%bres(:,:,3)
          call a%relaxVector('AITKE',a%sldale_dispRelax,a%sldale_dispRelax_max,r_i1,r_i21)
      endif

      w_r = a%sldale_dispRelax
      do ipoin=1,npoin
          do idime=1,a%ndofbc
              if(a%kfl_fixno(idime,ipoin)==4) then
                  !d_i2                 = w_r*d_i2_guess             + (1.0_rp-w_r)*d_i1
                  a%bvess(idime,ipoin,1)= w_r*a%bdisp(idime,ipoin,1) + (1.0_rp-w_r)*a%bdisp(idime,ipoin,2)
              end if
          end do
      end do  
   endif

end subroutine
