subroutine nsi_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* Nstinc/nsi_endite
   ! NAME 
   !    nsi_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip)   :: itask,ndime
   
   select case(itask)
   case(0)

      !Update unknowns
      call a%Mesh%GetNdime(ndime)

      !Newton Raphson, update the velocity
      if (a%kfl_linea==2 .and. a%istep>1) then
         a%unkno(1:ndime,:) = a%unkno(1:ndime,:)*0.5 + a%veloc(:,:,1)*(1.0_rp-0.5)
      endif
      
      !Subrelaxation only for the velocity
      a%veloc(:,:,1) = a%unkno(1:ndime,:)*a%subrelax + a%veloc(:,:,1)*(1.0_rp-a%subrelax)
      a%press(:,1)   = a%unkno(ndime+1,:)*a%subrelax + a%press(:,1)*(1.0_rp-a%subrelax)

   case(1)
      !Update the subscale and/or compute residual projection
      call a%EndElmope('Endite')
      call a%EndBouope('Endite')
      
   case(2)
      !Update u(n,i-1,*) <-- u(n,i,*)
      a%veloc(:,:,2) = a%veloc(:,:,1)
      ! Update p(n,i-1,*) <-- p(n,i,*)
      a%press(:,2) = a%press(:,1)

   case(4)
       !Update coupling values
       a%veloc_cp(:,:) = a%veloc(:,:,1)
       a%press_cp(:) = a%press(:,1)
   end select  

end subroutine nsi_endite
