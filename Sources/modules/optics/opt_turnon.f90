subroutine opt_Turnon(a)
   use typre
   use Mod_Optics
   implicit none

   class(OpticsProblem) :: a
   
   call a%Timer%Total%Tic
   call a%Timer%Turnon%Tic
   
   call a%SetExmod
   call a%OpenFiles
   if (a%MPIrank == a%MPIroot) then
         
      !Read Physical Problem
      call a%Reaphy
      
      !Read Numerical Treatment
      call a%Reanut
      
      !Read Output Strategy
      call a%Reaous
         
   endif
   call a%ReaMPI
  
   !Warnings and errors
   call a%Outerr(1)

   !Memory allocation
   call a%Memall
   
   !Warnings and errors
   call a%Outerr(2)
  
   call a%CloseDataFile
   
   call a%Timer%Turnon%Toc
   call a%Timer%Total%Toc


end subroutine opt_Turnon
