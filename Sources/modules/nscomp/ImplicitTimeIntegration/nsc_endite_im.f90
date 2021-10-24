subroutine nsc_endite_im(a,itask)
   use typre
   use def_parame
   use Mod_NSCompressibleImplicit
   
   implicit none
   class(NSCompressibleImplicitProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: ndime

   real(rp), parameter   :: &
     zensi = epsilon(0.0_rp)        ! zero

   select case(itask)

   case(0)

      !Update unknowns
      call a%Mesh%GetNdime(ndime)

      !Subrelaxation
      a%densf(:,1)   = a%unkno(1,:)        *a%subrelax + a%densf(:,1)  *(1.0_rp-a%subrelax)    
      a%momen(:,:,1) = a%unkno(2:ndime+1,:)*a%subrelax + a%momen(:,:,1)*(1.0_rp-a%subrelax)  
      a%energ(:,1)   = a%unkno(a%ndofn,:)  *a%subrelax + a%energ(:,1)  *(1.0_rp-a%subrelax)    

   case(1)

      call a%EndElmope('Endite')
      
      if(minval(a%densf(:,1))<(-zensi)) then
         call runend('Nsc_EndIte: Density is zero')
      end if

   case(2)

      !Update var(n,i-1,*) <-- var(n,i,*)
      a%densf(:,2) = a%densf(:,1)
      a%momen(:,:,2) = a%momen(:,:,1)
      a%energ(:,2) = a%energ(:,1)

   end select  

end subroutine nsc_endite_im
