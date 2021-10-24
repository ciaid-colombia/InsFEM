subroutine php_doiter(a)
   !This routine solves an iteration of the physical problem equations
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   
   call a%Timer%Total%Tic
   call a%Timer%Doiter%Tic
   
   if(a%kfl_stead==0) then
      call a%Begite
      do while(a%kfl_goite==1)
         call a%Solite
         call a%Endite(0)
         call a%Endite(1)
      enddo   

      call a%Endite(2)
      if(a%kfl_docoupconv) call a%Endite(4)

   end if
   
   call a%Timer%Doiter%Toc
   call a%Timer%Total%Toc
end subroutine php_doiter

