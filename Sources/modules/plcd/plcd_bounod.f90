subroutine plcd_bounod(a,itask)
   !This subroutine passes boundary conditions from boundary to nodes
   use typre
   use Mod_PLCD
   implicit none

   class(PLCDProblem) :: a
   integer(ip) :: itask, nboun, iboun, istage, pnodb, mnodb, inodb, ipoin
   integer(ip), allocatable :: lnodb(:)
   
   if (itask == 3) then
   
      call a%Mesh%GetNboun(nboun)
      call a%Mesh%GetMnodb(mnodb)
   
      call a%Memor%alloc(mnodb,lnodb,'lnodb','plcd_bounod')
   
      ! Dirichlet
      do istage = 1,a%NumberOfStages
         do iboun=1,nboun
            if(a%Stages(istage)%kfl_fixbo(iboun)==1) then
               a%Stages(istage)%kfl_fixbo(iboun) = 0
               call a%Mesh%GetLnodb(iboun,pnodb,lnodb)
               do inodb=1,pnodb
                  ipoin = lnodb(inodb)
                  if(maxval(a%Stages(istage)%kfl_fixno(:,ipoin)) == -1) then
                     a%Stages(istage)%kfl_fixno(:,ipoin)=1
                     a%Stages(istage)%bvess(:,ipoin)=a%Stages(istage)%bvnat(iboun)%a(:)
                  endif
               enddo
            endif
         enddo
      enddo
      
      call a%Memor%dealloc(mnodb,lnodb,'lnodb','plcd_bounod')

   end if

end subroutine
