subroutine ale_reaous(a,itask)
   use typre
   use def_parame
   use Mod_iofile
   use Mod_Alemov
   use Mod_sldAlemov
   implicit none

   class(AlemovProblem)  :: a
   integer(ip)           :: itask
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
 
   !Initializations
   if (itask == 0) then

      
       !Output section   
   elseif(itask == 1) then

       !Step or time for post-process
       if(words(1)=='POSTP') then
           if(words(2)=='ALEDI') then  
               if(words(3)=='STEPS') then
                   a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for displacement')
                   if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
               end if
           end if

           if(words(2)=='PSTRA') then  
               if(words(3)=='STEPS') then
                   a%npp_stepi(2) = a%Listener%getint('STEPS',1,'#Postprocess step interval for ppal strains - sldale')
                   if(a%Listener%exists('ATTIM')) a%pos_times(1:10,2)=param(4:13)
                   select type (a)
                   type is (sldAlemovProblem)
                       a%kfl_printPrincipalStresses=.true.
                   end select
               end if
           end if
       end if
  !Finalizations   
  elseif (itask == 100) then 

      if (a%nptra/=0) then
          if (a%MPIrank == a%MPIroot) then
              a%fil_trap = adjustl(trim(a%namda))//'.ale.tp'
              call iofile(zero,a%lun_trap,a%fil_trap,'ALE TRACKING OF POINTS')
          endif
      endif 


  end if

end subroutine
