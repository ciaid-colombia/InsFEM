subroutine php_endste(a,kfl_gotim)
   !This routine ends a time step.
   use typre
   use def_parame
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   integer(ip), intent(out) :: kfl_gotim
   
   call a%Timer%Total%Tic
   call a%Timer%Endste%Tic
   
   if(a%kfl_stead==0) then
      !Compute convergence residual of the time evolution 
      call a%Cvgunk(three)
   end if
   
   !Specific Pre-CrankNicolson
   call a%SpecificEndste(one)
   
   !CrankNicolson, End of Step Update
   if(a%kfl_stead==0.and.a%kfl_timei==1) then
      if((a%kfl_tsche_1st_current=='CNOBS').OR.(a%kfl_tsche_1st_current=='CN   ').OR.(a%kfl_tsche_1st_current=='TR   ')) then
         call a%SpecificCrankNicolsonEndste
      end if
   end if
   
   !Specific Post-CrankNicolson
   call a%SpecificEndste(two)
   
   !Output results corresponding to the end of a time step
   !call a%Timer%Output%Tic
   !Set already postprocessed to 0 (only at the end of a time step)
   a%pos_alrea = 0
   !call a%Output(zero)
   !call a%Timer%Output%Toc
   
   ! If not steady, go on
   if(a%kfl_stead==0.and.a%kfl_timei==1) kfl_gotim = 1

   call a%Timer%Endste%Toc
   call a%Timer%Total%Toc

      
   call a%OutputTimes


end subroutine php_endste

subroutine php_Output(a)
    !This routine ends a time step.
    use typre
    use def_parame
    use Mod_PhysicalProblem
    implicit none
    class(PhysicalProblem) :: a

   !Output results corresponding to the end of a time step
   call a%Timer%Output%Tic
   call a%Output(zero)
   call a%Timer%Output%Toc

end subroutine

subroutine php_OutputTimes(a)
   use typre
   use Mod_PhysicalProblem
   class(PhysicalProblem) :: a
   
   integer(ip) :: itime
   
   real(rp)    :: cpu_tim(14)
   
   !Brief computing times
   call php_TimesToCpuModulArray(a,cpu_tim)
   if (a%MPIrank == a%MPIroot) write(a%lun_outpu,700) a%namod, (cpu_tim(itime),itime=1,14)

700 format(  10x,a6,':',14(2x,e9.3))
   
end subroutine
