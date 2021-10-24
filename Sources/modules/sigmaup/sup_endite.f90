subroutine sup_endite(a,itask)
   !-----------------------------------------------------------------------
   !****f* ThreeField/sup_endite
   ! NAME 
   !    sup_endite
   ! DESCRIPTION 
   !    This routine checks convergence and performs updates at:
   !    - itask=1 The end of an internal iteration
   !    - itask=2 The end of the external loop iteration
   !-----------------------------------------------------------------------
   use typre
   use def_parame
   use Mod_ThreeField  
   implicit none
   class(ThreeFieldNSProblem) :: a
   integer(ip) :: itask
   integer(ip) :: ndime,ntens
   real(rp)    :: subres,subreu,subrep,subre

   select case(itask)

   case(0)
      !Update unknowns
      call a%Mesh%GetNdime(ndime)
      ntens=(ndime-1)*(ndime-1)+2
      subre= a%subrelax
      subres = subre
      subreu = subre
      subrep = subre
      a%sigma(:,:,1) = a%unkno(1:ntens,:)*subres + a%sigma(:,:,1)*(1.0_rp-subres)
      a%veloc(:,:,1) = a%unkno(ntens+1:ntens+ndime,:)*subreu + a%veloc(:,:,1)*(1.0_rp-subreu) 
      a%press(:,1)   = a%unkno(ntens+ndime+1,:)*subrep + a%press(:,1)*(1.0_rp-subrep) 
   case(1)
      !Update the subscale and/or compute residual projection
      call a%EnditeElmope('Endite')
      call a%EndBouope('Endite')

   case(2)
      !Update s(n,i-1,*) <-- s(n,i,*)
      a%sigma(:,:,2) = a%sigma(:,:,1)
      !Update u(n,i-1,*) <-- u(n,i,*)
      a%veloc(:,:,2) = a%veloc(:,:,1)
      ! Update p(n,i-1,*) <-- p(n,i,*)
      a%press(:,2) = a%press(:,1)

   case(4)
      !Update coupling values
      a%sigma_cp(:,:) = a%sigma(:,:,1)    
      a%veloc_cp(:,:) = a%veloc(:,:,1)
      a%press_cp(:)   = a%press(:,1)
   end select  

end subroutine sup_endite
