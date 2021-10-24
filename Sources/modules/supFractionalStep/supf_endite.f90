subroutine supf_endite(a,itask)
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
   use Mod_SUPFractionalStep  
   implicit none
   class(SUPFractionalStepProblem) :: a
   integer(ip) :: itask   
   integer(ip) :: ndime,ntens,aux,aux1
   real(rp) :: subres,subreu,subrep,subre
   
   
   call a%Timer%Endite%Tic   

   select case(itask)
  
   case(2)
   
      !Compute convergence residual of the external iteration
      call a%Cvgunk(two)
      
      
      !En esta parte se puede usar en el extrapolador de alto orden
      
      !Update s(n,i-1,*) <-- s(n,i,*)
      a%sigma(:,:,2) = a%sigma(:,:,1)
      !Update u(n,i-1,*) <-- u(n,i,*)
      a%veloc(:,:,2) = a%veloc(:,:,1)
      ! Update p(n,i-1,*) <-- p(n,i,*)
      a%press(:,2) = a%press(:,1)

   end select  
   
   call a%Timer%Endite%Toc

end subroutine supf_endite
