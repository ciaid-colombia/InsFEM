subroutine nsc_pr_endite(a,itask)
   use typre
   use Mod_NSCompressiblePrimitive
   
   implicit none
   class(NSCompressiblePrimitiveProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: ndime

   select case(itask)

   case(0)

      !Update unknowns
      call a%Mesh%GetNdime(ndime)

      !Subrelaxation
      a%press(:,1)   = a%unkno(1,:)        *a%subrelax + a%press(:,1)  *(1.0_rp-a%subrelax)    
      a%veloc(:,:,1) = a%unkno(2:ndime+1,:)*a%subrelax + a%veloc(:,:,1)*(1.0_rp-a%subrelax)  
      a%tempe(:,1)   = a%unkno(a%ndofn,:)  *a%subrelax + a%tempe(:,1)  *(1.0_rp-a%subrelax)    

   case(1)

      call a%EndElmope('Endite')
      
   case(2)

      !Update var(n,i-1,*) <-- var(n,i,*)
      a%press(:,2) = a%press(:,1)
      a%veloc(:,:,2) = a%veloc(:,:,1)
      a%tempe(:,2) = a%tempe(:,1)

   end select  

end subroutine nsc_pr_endite
