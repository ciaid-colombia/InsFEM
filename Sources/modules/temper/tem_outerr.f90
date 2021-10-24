subroutine tem_outerr(a)
   use typre
   use Mod_Temperature
   use def_parame
   use Mod_int2str
   implicit none
   class(TemperatureProblem) :: a
   
   integer(ip) :: ierro = 0

   if (a%kfl_trasg /= 0 .and. 'ForceClosedRule' .eq. a%EndLoopQuadrature) then
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_outpu,*)
         write(a%lun_outpu,*) 'TemperatureProblem: If tracking of subgrid scales is used, quadratures at the end of step and end of iteration must be default'
         write(a%lun_outpu,*) 'Setting EndLoopQuadrature to Default'
         write(a%lun_outpu,*)
      endif
      a%EndLoopQuadrature = 'DefaultRule'
      !a%kfl_ProjectionType = 1
   endif
   
   


   !Write Errors and Warnings
   if(ierro==1) then
      call runend(adjustl(trim(int2str(ierro)))//' ERROR HAS BEEN FOUND')
   else if(ierro>=2) then
      call runend(adjustl(trim(int2str(ierro)))//' ERRORS HAVE BEEN FOUND')
   end if

! Formats
100 format(5x,'ERROR:   ',a)
101 format(5x,'WARNING: ',a)



end subroutine
