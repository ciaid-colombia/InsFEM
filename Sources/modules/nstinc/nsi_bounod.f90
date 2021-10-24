subroutine nsi_bounod(a,itask)
   !This subroutine passes boundary conditions from boundary to nodes
   use typre
   use Mod_NavierStokes
   implicit none

   class(NavierStokesProblem) :: a
   integer(ip) :: itask
   
   integer(ip) :: lnodb(a%gmnodb),pnodb,inodb,ipoin,ndime,idime
   
   interface 
      subroutine nsi_BCNoSlipToNodes(a,iboun)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: iboun
      end subroutine
   end interface
      
   
   if (itask == 1) then
      a%kfl_fixrs(a%gipoin) = a%kfl_bours(a%giboun)
      
   elseif(itask == 2) then

     call a%Mesh%GetLnodb(a%giboun,pnodb,lnodb)

     ! No slip wall
     if(a%kfl_fixbo(a%giboun)==7) then
        !We pass noslip boundary conditions from boundaries to nodes
        call nsi_BCNoSlipToNodes(a,a%giboun)
        
     !Wall law / Slip-wall / Symmetry   
     else if((a%kfl_fixbo(a%giboun)==3).OR.(a%kfl_fixbo(a%giboun)==4)) then
        do inodb=1,pnodb
           ipoin = lnodb(inodb)
           if(maxval(a%kfl_fixno(:,ipoin)) <= 0) then
              a%kfl_fixno(1,ipoin)=1
              a%kfl_fixno(2:a%ndofbc,ipoin) = 0
              a%bvess(1:a%ndofbc,ipoin,1)=0.0_rp
              a%kfl_fixrs(ipoin)=-1
              if(a%kfl_conbc==0) a%kfl_funno(ipoin)=0
           end if
        end do        
     end if
   endif
end subroutine

subroutine nsi_BCNoSlipToNodes(a,iboun)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   integer(ip) :: iboun
   
   integer(ip) :: lnodb(a%gmnodb),pnodb,inodb,ipoin
   
   call a%Mesh%GetLnodb(iboun,pnodb,lnodb)
   do inodb=1,pnodb
      ipoin = lnodb(inodb)
      !if(maxval(a%kfl_fixno(:,ipoin)) == -1) then
         a%kfl_fixno(1:a%ndofbc,ipoin)=1
         a%bvess(1:a%ndofbc,ipoin,1)=0.0_rp
         a%kfl_fixrs(ipoin)=0
         if(a%kfl_conbc==0) a%kfl_funno(ipoin)=0
      !end if
   end do  
end subroutine
