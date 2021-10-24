subroutine ale_cvgunk(a,itask)
   use MPI
   use typre
   use def_parame
   use Mod_Alemov
   implicit none
   class(AlemovProblem) :: a
   integer(ip) :: itask

   integer(ip), save       :: ipasscp=0
   real(rp)                :: riale
   integer(ip)             :: npoinLocal,ndime,ipoin
   integer                 :: ierr
   
   select case(itask)
   
   !Inner iterations
   case(1)
            
         a%kfl_goite=0
   
   !Outer iterations
   case(2)
      
   
   !End of step convergence
   case(3)

   !Coupling convergence
   case(4)

       a%kfl_coupconv= .true.
       CALL MPI_BCAST(a%kfl_coupconv, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)

  end select

  !Formats
  100 format('$ ','       Time','      Inner',&
          &      '       Current','      Displacement',/,&
          & '$ ','       step','  iteration',&
          &      '          time','      residual')
  101 format(4x,i9,2x,i9,20(2x,e12.6))
  102 format('$ >>>  ALEMOV IS STATIONARY AT TIME STEP ',i5)

end subroutine
