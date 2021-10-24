subroutine LocalGraph(a)
   use Mod_Mesh
   use Mod_BuildGraph
   implicit none
   class(FemMesh) :: a
   
   call a%Timer%LocalGraph%Tic
   
   !Usual case
   if (a%kfl_HangingNodes .eqv. .false.) then
      call BuildElpo(a%pelpo,a%lelpo,a%npoin,a%nelem,a%pnods,a%lnods,a%nelty,a%nnode,a%Memor)     !Compute glelpo and gpelpo
      call BuildGraph(a%ia,a%ja,a%pelpo,a%lelpo,a%npoin,a%nelem,a%pnods,a%lnods,a%nelty,a%nnode,a%Memor)
   
   !Hanging nodes
   else
      call a%BuildHangingGraph
   endif
   
   call a%Timer%LocalGraph%Toc
   
end subroutine 
