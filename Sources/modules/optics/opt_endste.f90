subroutine opt_endste(a,kfl_gotim)
   !This routine ends a time step.
   use typre
   use def_parame
   use Mod_Timer
   use Mod_Memor
   use Mod_Listen
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_Optics
   implicit none
   class(OpticsProblem) :: a
   integer(ip), intent(out) :: kfl_gotim
   
   real(rp)    :: cputim1, cputim2, cputim3, cputim4

   interface
      subroutine opt_coupling_outerr(a)
         use Mod_Optics
         implicit none
         class(OpticsProblem) :: a
      end subroutine
      
      subroutine opt_elmope(a)
         use Mod_Optics
         implicit none
         class(OpticsProblem) :: a
      end subroutine
      
      subroutine opt_EndsteComputations(a)
         use Mod_Optics
         implicit none
         class(OpticsProblem) :: a
      end subroutine
   
   end interface
   
   call a%Timer%Total%Tic
   call a%Timer%Endste%Tic
   
   !Update the step count (usually done in begste, but here it is NULLSUB
   a%istep = a%istep+1
   
   !!Check if the coupling between modules is done
   call opt_coupling_outerr(a)
   
   !Compute nodal cn2 values
   call opt_elmope(a)
   
   !Output results corresponding to the end of a time step
   !call a%Timer%Output%Tic
   !call a%Output(zero)
   !call a%Timer%Output%Toc
   
   !Compute and output ray quadratures
   call opt_EndsteComputations(a)
   
   !Do not modify kfl_gotim
   kfl_gotim = kfl_gotim

   call a%Timer%Endste%Toc
   call a%Timer%Total%Toc
   
   call a%OutputTimes   
   


end subroutine opt_endste
