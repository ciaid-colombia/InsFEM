module Mod_caseSetste
   use def_parame
   use MPI
   use typre
   use Mod_GeneralCase
   use Mod_caseVariables
   use Mod_DistributedContainer
   use Mod_DC_Driver
   use Mod_DriverInterface
   implicit none
   
   !Needs to be save, in this way the internal subroutines can also access it!
   type(caseVariables), pointer :: c => NULL()
   type(masterVariables), pointer :: m => NULL()
   
contains
   subroutine LoopGetTimeStep(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%GetTimeStep(c)
   end subroutine
   
   subroutine LoopSetTimeStep(myDC)
      class(DistributedContainer), pointer :: myDC
      class(DriverInterface), pointer :: myDriver => NULL()
      
      call ExtractDriver(myDC,myDriver)
      call myDriver%SetTimeStep(c)
   end subroutine
   
   !This subroutine sets a variable time step in time   
   subroutine case_Setgts(a)
      class(GeneralCase), target :: a

      real(rp)    :: timea

      if(m%kfl_timco==1) then
         !Use minimum of critical time steps
         if(m%dtinv>1.0d-10) m%dtime = 1.0_rp/m%dtinv
      else if(m%kfl_timco==3) then
         !Time step as a function of the step
         if(nint(m%timfu(1))<=m%istep.and.m%istep<=nint(m%timfu(2))) then 
           m%dtime=m%timfu(3)*real(m%istep*m%istep)+m%timfu(4)*real(m%istep)+m%timfu(5)
         else if (m%istep>nint(m%timfu(2))) then
           timea=m%timfu(2)
           m%dtime=m%timfu(3)*timea*timea+m%timfu(4)*timea+m%timfu(5)
         else if (m%istep<nint(m%timfu(1))) then
           timea=m%timfu(1)
           m%dtime=m%timfu(3)*timea*timea+m%timfu(4)*timea+m%timfu(5)
         end if
         if(m%dtime>0.0) then
           m%dtinv = 1.0_rp/m%dtime
         else
            call runend('SETGTS: TIME STEP IS <= 0 WRONG FUNCTION?')
         end if
      else if(m%kfl_timco==4) then
         !Time step for predictor corrector
         call runend('Setgts: kfl_timco=4')
         !dtime = muldt/dtin0
         !dtinv = 1.0_rp/dtime
      end if
      if (m%MPIrank == m%MPIroot) write(m%lun_outpu,100) m%dtime
     
      100 format(5x,'     CURRENT TIME STEP dt = ',e12.6)
     
   end subroutine 

end module

subroutine case_Setste(a)
   !-----------------------------------------------------------------------
   !****f* master/Setste
   ! NAME
   !    Setste
   ! DESCRIPTION
   !    a routine computes the global time step when two or more problems
   !    are running simultaneously. Each module is called to determine its
   !    own critical time step and the minimum over the modules is selected.
   !-----------------------------------------------------------------------
   use Mod_caseSetste
   implicit none
   class(GeneralCase), target :: a

   c => a%caseVars
   m => a%caseVars%masterVars
   
   m%istep = m%istep + 1
   if(m%kfl_timco==1) m%dtinv = 0.0_rp

   m%dtinv = 1.0_rp/m%dtime
   
   !Loop through drivers and get their critical time step
   call a%DriverList%LoopList(LoopGetTimeStep)
      
   !Choose the next time step and update current time
   
   call case_Setgts(a)
   m%ctime = m%ctime + 1.0_rp/m%dtinv
   
   if (m%MPIrank == m%MPIroot) write(m%lun_outpu,100) m%istep,m%ctime

   !Loop through drivers and tell them the new time step
   call a%DriverList%LoopList(LoopSetTimeStep)

   !Formats.
   100 format(/,5x,&
         '>>>  TIME STEP NUMBER: ',i5,&
         '     STARTED AT t = ',e12.6)

end subroutine 
