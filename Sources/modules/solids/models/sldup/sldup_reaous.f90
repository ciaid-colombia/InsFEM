subroutine sldup_reaous(a,itask)
   use typre
   use Mod_UPSolids
   implicit none

   class(UPSolidsProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   integer(ip) :: ndime,ipara
   
   
   !Initializations
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then  

   elseif (itask == 1) then  

      ! Step or time for post-process
      if(words(1)=='POSTP') then

         if(words(2)=='J2STR') then
             a%kfl_printJ2Stresses = .true.
             if(words(3)=='STEPS') then
                a%npp_stepi(15) = a%Listener%getint('STEPS',1,'#Postprocess step interval for J2 stresses')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,15)=param(4:13)
             end if
         elseif(words(2)=='RESID') then
             a%kfl_printResiduals = .true.
             if(words(3)=='STEPS') then
                 a%npp_stepi(16) = a%Listener%getint('STEPS',1,'#Postprocess step interval for Residual')
                 if(a%Listener%exists('ATTIM')) a%pos_times(1:10,16)=param(4:13)
             endif  
         elseif(words(2)=='VESGS') then
             if(words(3)=='STEPS') then
                a%npp_stepi(17) = a%Listener%getint('STEPS',1,'#Postprocess step interval for vesgs ')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,17)=param(4:13)
             end if
         elseif(words(2)=='PRSGS') then
             if(words(3)=='STEPS') then
                a%npp_stepi(19) = a%Listener%getint('STEPS',1,'#Postprocess step interval for prsgs ')
                if(a%Listener%exists('ATTIM')) a%pos_times(1:10,19)=param(4:13)
             end if
          end if

      end if


   elseif (itask == 100) then   

      if (a%kfl_repro == 0 .and. a%npp_stepi(16)>0) then
         a%npp_stepi(16) = 0
         write(*,*) 'WARNING: Pick OSS residual to print residual projection'
      end if         
      

   endif

   call a%ReaousExtend(itask)

   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine sldup_reaous

