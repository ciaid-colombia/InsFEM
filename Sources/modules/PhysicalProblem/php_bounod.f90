subroutine php_bounod(a)
!-----------------------------------------------------------------------
!    Passes conditions on boundaries to nodal conditions.
!-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip) :: nboun,pnodb,ipoin,idime,iboun,mnodb,inodb
   integer(ip), allocatable :: lnodb(:)
   character(150) :: outstr
   
   call a%Mesh%GetMnodb(mnodb)
   call a%Mesh%GetNboun(nboun)
   a%gmnodb = mnodb
   
   !Output
   outstr = adjustl(trim(a%exmod))//'_bounod'
   
   call a%Memor%alloc(mnodb,lnodb,'lnodb',outstr)
   
   ! Loop over boundaries
   do iboun=1,nboun
      !For Specific call
      a%giboun = iboun
      
      ! Dirichlet
      if(a%kfl_fixbo(iboun)==1) then
         a%kfl_fixbo(iboun) = 0
         call a%Mesh%GetLnodb(iboun,pnodb,lnodb)
         do inodb=1,pnodb
            ipoin = lnodb(inodb)
            if(maxval(a%kfl_fixno(:,ipoin)) == -1) then
               a%kfl_fixno(:,ipoin)=1
               do idime=1,a%ndofbc
                  a%bvess(idime,ipoin,1)=a%bvnat(iboun)%a((inodb-1)*a%ndofbc+idime)
               end do
               if(a%kfl_conbc==0) then
                  a%kfl_funno(ipoin)=a%kfl_funbo(iboun)
                  do idime=1,a%ndofbc
                     a%bvess(idime,ipoin,2)=a%bvess(idime,ipoin,1)
                  end do
               end if
               !Specific call
               a%gipoin = ipoin
               call a%SpecificBounod(1)
            end if
         end do
      endif
      
      call a%SpecificBounod(2)
      
   end do
   
   call a%SpecificBounod(3)
   
   call a%Memor%dealloc(mnodb,lnodb,'lnodb',outstr)
   
end subroutine php_bounod
