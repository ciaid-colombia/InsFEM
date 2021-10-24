subroutine sldup_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* solids/sldsup_endite
   ! NAME 
   !    sldsup_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem),target :: a
   integer(ip) :: itask,ndime,np,idime,ipoin,icomp
   real(rp) :: w,w_u,w_p
   integer(ip) :: u1,uf,s1,sf,p1,bc
   
   select case(itask)

   case(0)

       !Update unknowns
       call a%Mesh%Getnpoin(np)
       call a%Mesh%GetNdime(ndime)

       call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

       w = a%sld_dispRelax

       w_u = w
       w_p = w
       a%disp (:,:,1) = a%unkno(u1:uf,:)*w_u + a%disp(:,:,1) *(1.0_rp-w_u) 
       a%press(:,1)   = a%unkno(p1,:)*w_p    + a%press(:,1)  *(1.0_rp-w_p) 

       !a%disp (:,:,1) = a%unkno(u1:uf,:) + a%disp(:,:,1)
       !a%press(:,1)   = a%unkno(p1,:)    + a%press(:,1) 

   case(1)

      if (a%kfl_timei == 1)     call a%sld_calcVelandAccel(a%dtinv)

      call a%EndElmope('Endite')

   case(2)

       !Update u(n,i-1,*) <-- u(n,i,*)
       a%disp(:,:,2)  = a%disp(:,:,1)
       !Update p(n,i-1,*) <-- p(n,i,*)
       a%press(:,2)   = a%press(:,1)

   case(4)

       !Update coupling values
       w = a%sld_dispRelax
       a%disp_cp(:,:)  = w*a%disp(:,:,1)  + (1.0_rp - w)*a%disp_cp(:,:)
       a%press_cp(:)   = w*a%press(:,1)   + (1.0_rp - w)*a%press_cp(:)

  end select  

end subroutine sldup_endite
