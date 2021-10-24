subroutine php_begste(a)
   !-----------------------------------------------------------------------
   !    This routine prepares for a new time step 
   !-----------------------------------------------------------------------
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   integer(ip) :: icomp,ndime,npoin
   character(5) :: kfl_tsche_1st_tmp,kfl_tsche_2nd_tmp
   if(a%kfl_stead==1) return
   call a%Timer%Total%Tic
   call a%Timer%Begste%Tic
   
   if(a%kfl_timei==1) then
      a%istep = a%istep+1
   endif
   
   a%kfl_elmat_current = a%kfl_elmat_datafile
   
   a%kfl_tsche_change = 0
   kfl_tsche_1st_tmp = a%kfl_tsche_1st_current
   kfl_tsche_2nd_tmp = a%kfl_tsche_2nd_current
   
   if (a%istep<=a%neule.and.a%kfl_timei==1.and.a%kfl_rstar==0) then
      a%kfl_tsche_1st_current = 'BDF1 '
   else
      a%kfl_tsche_1st_current = a%kfl_tsche_1st_datafile
   endif
   
   if (a%istep<=a%neule_2nd.and.a%kfl_timei==1.and.a%kfl_rstar==0) then
      a%kfl_tsche_2nd_current = 'BDF1 '
   else
      a%kfl_tsche_2nd_current = a%kfl_tsche_2nd_datafile
   endif
   
   !call runend('php_begste: provisional stop for testing')
   
   !Set Time step
   if(a%kfl_stead==1) then
      a%dtinv  = 0.0_rp
      a%dtinv2 = 0.0_rp
   end if
   
   !Update boundary conditions
   call a%Updbcs
   
   call a%SpecificBegste
   
   if(a%kfl_tsche_1st_current=='CNOBS') then
      a%dtinv = a%dtinv*2.0_rp
      a%bctime = a%ctime - a%dtime*0.5_rp
   elseif (a%kfl_tsche_1st_current=='CN   ') then
      !we do not modify dtinv because we now use coefficients for new crank-nicolson
      a%bctime = a%ctime - a%dtime*0.5_rp
   end if
   
   if(a%kfl_timei==0) then
      a%dtinv  = 0.0_rp
      a%dtinv2 = 0.0_rp
   end if
   if(a%MPIrank == a%MPIroot .and. a%kfl_timei/=0.and.a%kfl_stead/=1) then
       if(a%dtcri < 1e-6) then
           write(a%lun_outpu,110) a%namod, a%dtcri, 0.0
       else
           write(a%lun_outpu,110) a%namod, a%dtcri,a%dtime/a%dtcri
       end if
   end if
   
   if ((kfl_tsche_1st_tmp/=a%kfl_tsche_1st_current).OR.(kfl_tsche_2nd_tmp/=a%kfl_tsche_2nd_current)) then
      a%kfl_tsche_change = 1
   end if
   
   if ((a%kfl_tsche_change==1).OR.(a%istep==1)) then
      a%kfl_elmat_current = 1
   end if

   call a%Timer%Begste%Toc
   call a%Timer%Total%Toc
   
   ! Formats.
   100 format(//,'============================================================', &
         //,5x,'>>>  TIME STEP NUMBER: ',i5)
         
   110 format(&
         5x,5x,a6,':',/,&
         5x,'     - Critical Time step dtc: ',e12.6,/,&
         5x,'     - Ratio dt/dtc:           ',e12.6 )
   120 format(&
         5x,5x,a6,':',/,&
         5x,'     - Critical Time step dtc: ',e12.6,/,&
         5x,'     - Using local time step')   

end subroutine php_begste
