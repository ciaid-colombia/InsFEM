subroutine supf_pbcpre(a)
   use typre
   use Mod_SUPFractionalStep   
   use Mod_Element
   implicit none
   class(SUPFractionalStepProblem) :: a
   
   integer(ip) :: iffix_aux,iffix_aux2,pelty,iboun,idime,inodb,ipoin,nboun,ndime,ipoin2
   class(FiniteElement), pointer :: e => NULL()
   
   integer(ip) :: kfl_perio
   real(rp), pointer :: exnor(:,:) => NULL()
   integer(ip) :: ibopo,npoin,auxtens

   real(rp)    :: raux
   
   a%kfl_fixpr = 0
   
   !Only for open flow
   if (a%kfl_confi <= 0) then

      call a%Mesh%GetNboun(nboun)
      call a%Mesh%GetNdime(ndime)
      
      !Memory allocation
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsf_pbcpre')
      !Loop over boundaries
      boundaries: do iboun=1,nboun
         !Load Element
         call a%Mesh%BoundaryLoad(iboun,e)
         
         if (a%kfl_fixbo(iboun) == 2) then
            a%kfl_fixpr(e%lnodb(1:e%pnodb)) = 2
            a%bvpress(e%lnodb(1:e%pnodb)) = a%bvnat(iboun)%a(1)
            
         elseif (a%kfl_fixbo(iboun) == 0) then
         
            iffix_aux2 = 0
            do inodb = 1,e%pnodb
               ipoin = e%lnodb(inodb)
               
               iffix_aux = 1
               do idime = 1,ndime
               auxtens=(e%ndime-1)*(e%ndime-1)+2
                  if (a%kfl_fixno(auxtens+idime,ipoin) == 1) iffix_aux = 0
               enddo

               if (iffix_aux == 1) then
                  iffix_aux2 = 1
               endif
            enddo
            
            if (iffix_aux2 == 1) then
               !Compute the boundary normal
               e%igaub = 1
               call e%bounor
            
               do inodb = 1,e%pnodb
                  ipoin = e%lnodb(inodb)
                  
                  if (a%kfl_fixpr(ipoin) == 0) then
                     !Make sure that the point normal is not orthogonal to the 
                     !face normal (this would cause false positives in periodic BC, for instance)
                     call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)        
                     if (ibopo /= 0) then
                        raux = dot_product(e%baloc(1:e%ndime,e%ndime),exnor(1:e%ndime,1))
                        if (abs(raux) > 1e-6) then
                           a%kfl_fixpr(ipoin) = 1
                        endif
                     endif
                  endif   
               enddo
            endif
         elseif (a%kfl_fixbo(iboun) == 1) then
            !Do nothing, dirichlet
         elseif (a%kfl_fixbo(iboun) == 3) then
            !Do nothing, wall law
         elseif (a%kfl_fixbo(iboun) == 4) then
            !Do nothing, slip   
         else
            !Do nothing, free traction
            !write(*,*) 'kfl_fixbo(iboun) = ',a%kfl_fixbo(iboun)
            !call runend('nsf_pbcpre: kfl_fixbo(iboun) =  this boundary condition on boundaries is not implemented for fractional step')
         endif   
      enddo boundaries
      call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsf_pbcpre')       
   endif   

  
   
end subroutine