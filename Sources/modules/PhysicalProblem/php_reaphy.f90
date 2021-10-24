subroutine php_reaphy(a)
   use typre
   use Mod_PhysicalProblem
   implicit none
   
   class(PhysicalProblem) :: a
   
   integer(ip) :: istat,npara,ipara
   real(rp)    :: dummr,gnorm
   
   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   character(150) :: outstr
   
   call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   
   outstr = adjustl(trim(a%exmod))//'_REAPHY'
   
   ! Initializations (defaults)
   a%kfl_stead = 0_ip                         ! Steady-state
   a%kfl_timei = 0_ip                         ! Time integration off
   a%dtinv     = 1.0_rp                       ! 1/dt=1.0
   
   call a%SpecificReaphy(0)
   
   call a%Listener%rewind
   call a%Listener%listen(outstr)
   do while(words(1)/='PHYSI')
      call a%Listener%listen(outstr)
   end do
   
   !Begin to read data
   do while(words(1)/='ENDPH')
      call a%Listener%listen(outstr)
      if(words(1)=='PROBL') then        ! Problem definition data
         call a%Listener%listen(outstr)
         do while(words(1)/='ENDPR')
            if(words(1)=='TEMPO') then                    ! Temporal evolution
               if(a%Listener%exists('ON   ')) a%kfl_timei = 1 
            
            !Analytical Solution
            else if(words(1)=='ANALY') then
               a%kfl_exacs = a%Listener%getint('ANALY',0,'#Exact solution')  
               a%expar = param(2:11)
            endif
            
            !Problem data
            call a%SpecificReaphy(1)
            
             call a%Listener%listen(outstr)
         end do
      else if(words(1)=='PROPE') then   ! Properties
         call a%Listener%listen(outstr)
         do while(words(1)/='ENDPR')
         
            !Properties
            call a%SpecificReaphy(2)
            
            call a%Listener%listen(outstr)
         end do
      end if
   end do
   
   if(a%kfl_timei==0) then
     a%dtinv=1.0_rp
   else
     a%kfl_timei=1 
   end if
   
   call a%SpecificReaphy(100)
   
end subroutine
