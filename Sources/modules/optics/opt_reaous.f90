subroutine opt_reaous(a,itask)
   use typre
   use Mod_Listen
   use Mod_Optics
   implicit none

   class(OpticsProblem) :: a
   integer(ip) :: itask
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   integer(ip) :: istat,myIOSTat1,ndime,ibeam
   real(rp) :: r,vnor
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   if (itask == 0) then         !Initializations.
      !Default wavelength is 500 nanometers
      a%lambda = 500e-9
      
      !Default is do not do the average in 1 dimension
      a%kfl_Avg1DCn2 = 0
     
   elseif (itask == 1) then  

      ! Step or time for post-process
      if(words(1)=='POSTP') then
         if(words(2)=='CN2  ') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(1) = a%Listener%getint('STEPS',1,'#Postprocess step interval for CN2  ')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,1)=param(4:13)
            end if
         else if(words(2)=='CT2  ') then  
            if(words(3)=='STEPS') then
               a%npp_stepi(2) = a%Listener%getint('STEPS',1,'#Postprocess step interval for CT2  ')
               if(a%Listener%exists('ATTIM')) a%pos_times(1:10,2)=param(4:13)
            end if
         end if
      else if(words(1)=='NUMBE') then
         call a%Mesh%GetNdime(ndime)
         a%nbeams = int(a%Listener%param(1))
         allocate(a%beams(a%nbeams),STAT=istat)
         call a%Memor%allocObj(istat,'beams','opt_reaous',6*rp*a%nbeams)
         do ibeam = 1,a%nbeams
            call a%Listener%listen('opt_reaous')
            call RANDOM_NUMBER(r);
            r = (r -0.5)*1e-6;
            a%beams(ibeam)%origin(1:ndime) = a%Listener%param(1:ndime)+r
            a%beams(ibeam)%direct(1:ndime) = a%Listener%param(ndime+1:2*ndime)
            call vecnor(a%beams(ibeam)%direct,ndime,vnor,2_ip)
            a%beams(ibeam)%direct = a%beams(ibeam)%direct/vnor
            a%beams(ibeam)%length = a%Listener%param(2*ndime+1)
            !read(a%Listener%nunit,*,IOStat=myIOStat1)  a%beams(ibeam)%origin(1:ndime), a%beams(ibeam)%direct(1:ndime)      
         enddo 
      elseif(words(1)=='WAVEL') then
         a%lambda = a%Listener%param(1)
      elseif(words(1)=='AVG1D') then
         if(a%Listener%exists('ON   ')) then
            a%kfl_Avg1DCn2 = 1
            if (a%Listener%exists('X    ')) then
               a%Avg1DIdime = 1
            elseif (a%Listener%exists('Y    ')) then
               a%Avg1DIdime = 2
            elseif (a%Listener%exists('Z    ')) then
               a%Avg1DIdime = 3
            endif 
         endif
            
      end if
      
      
   elseif (itask == 100) then   

   endif


   1 format(a)
   5 format(a,1x,i2)

   101 format(5x,'WARNING: ',a)

end subroutine opt_reaous
