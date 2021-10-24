subroutine SolveSystem(n,Matrix,RHS,Solution)
   use typre
   implicit none
   integer(ip) :: n
   real(rp) :: Matrix(n,n), RHS(n), Solution(n)

   integer(ip) :: info
   !Solution is auxiliary for dgesv
   !The solution comes out in RHS
   call DGESV( n, 1_ip, Matrix, n, Solution, RHS, n, INFO )
   Solution = RHS;
end subroutine