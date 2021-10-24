subroutine sld_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* lmachn/sld_endite
   ! NAME 
   !    sld_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem),target :: a
   integer(ip) :: itask
   real(rp), pointer     :: r_i1(:,:) => NULL(),r_i21(:,:) => NULL()
   real(rp) :: w
   
   select case(itask)

   case(0)

      !w = a%sld_dispRelax
      !a%disp(:,:,1) = w*a%unkno(:,:) + (1.0_rp - w)*a%disp(:,:,1)
      a%disp(:,:,1) = a%unkno(:,:)

   case(1)

      if (a%kfl_timei == 1)     call a%sld_calcVelandAccel(a%dtinv)

      if(a%sld_type== 'NONLI' ) call a%EndElmope('Endste')

  case(2)

      !Update u(n,i-1,*) <-- u(n,i,*)
      a%disp(:,:,2) = a%disp(:,:,1)

  case(4)

      !Update coupling values
      w = a%sld_dispRelax
      a%disp_cp(:,:) = w*a%disp(:,:,1) + (1.0_rp - w)*a%disp_cp(:,:)

  end select  

end subroutine sld_endite
