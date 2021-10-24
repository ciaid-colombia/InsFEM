subroutine SetupDopost(itask,istep,npost,npp_inits,npp_stepi,pos_alrea,dopost)
   use typre
   implicit none
   integer(ip) :: itask,istep,npost,npp_inits,npp_stepi(npost),dopost(npost),pos_alrea(npost)
   
   integer(ip) :: ipost
   
   dopost = 0
   
   if (itask == 0) then
      pos_alrea = 0
      do ipost = 1,npost
         if(istep>=npp_inits.and.npp_stepi(ipost)>0) then     
            if(mod(istep,npp_stepi(ipost))==0) then
               pos_alrea(ipost)=1
               dopost(ipost) = 1
            endif
         endif
      enddo
      
   elseif (itask == 1) then
      do ipost = 1,npost
         if(npp_stepi(ipost)/=0) then
            if (pos_alrea(ipost) == 0) dopost(ipost) = 1
         end if
      enddo
      
   endif

end subroutine