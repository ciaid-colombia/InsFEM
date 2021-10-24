module Mod_sup_Stats
   use typre
   use Mod_ThreeField
   use MPI
   implicit none

contains   


   subroutine sup_Stats(a)
      class(ThreeFieldNSProblem) :: a
      integer(ip) :: itask
     
      integer :: ierr
      real(rp)                :: tamin, tamax, tamea
      real(rp)                :: tasigmin, tasigmax, tasigmea
      real(rp)                :: remin, remax, remea

      
      !Endste   
      call MPI_REDUCE( a%tamea, tamea, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( a%tamin, tamin, 1, MPI_REAL8, MPI_MIN, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( a%tamax, tamax, 1, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
      
      call MPI_REDUCE( a%tasigmea, tasigmea, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( a%tasigmin, tasigmin, 1, MPI_REAL8, MPI_MIN, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( a%tasigmax, tasigmax, 1, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
      
      call MPI_REDUCE( a%remea, remea, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( a%remin, remin, 1, MPI_REAL8, MPI_MIN, a%MPIroot,a%MPIcomm, ierr )
      call MPI_REDUCE( a%remax, remax, 1, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

      if (a%MPIrank == a%MPIroot) then
         a%tamea = tamea/a%MPIsize
         a%tamax = tamax
         a%tamin = tamin
         
         a%tasigmea = tasigmea/a%MPIsize
         a%tasigmax = tasigmax
         a%tasigmin = tasigmin
         
         a%remea = remea/a%MPIsize
         a%remax = remax
         a%remin = remin
      endif
         
   end subroutine   

end module


subroutine supm_InitStats(a)
   use typre
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
   
   !Write info
   a%tamin=1.0e20_rp                  ! Minimum tau
   a%tamax=0.0_rp                     ! Maximum tau
   a%tamea=0.0_rp                     ! Mean tau
   a%remin=1.0e20_rp                  ! Minimum Re
   a%remax=0.0_rp                     ! Maximum Re
   a%remea=0.0_rp                     ! Mean Re
   a%nmean=0
   a%tasigmin=1.0e20_rp                  ! Minimum tau
   a%tasigmax=0.0_rp                     ! Maximum tau
   a%tasigmea=0.0_rp                     ! Mean tau
end subroutine


subroutine supm_InGaussStats_sup(a,acden,acvis,gpvno,chale,timom,tisig)
   use typre
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
   real(rp) :: acden,acvis,gpvno,chale(2),timom,tisig
   real(rp) :: reyno
      
   ! Compute cell reynolds number and store values.
   a%nmean = a%nmean+1
   reyno = acden*gpvno*chale(2)/acvis
   a%remin = min(a%remin,reyno)
   a%remax = max(a%remax,reyno)
   a%remea = a%remea + reyno

   a%tamin = min(a%tamin,timom)
   a%tamax = max(a%tamax,timom)
   a%tamea = a%tamea + timom
   
  
   a%tasigmin = min(a%tasigmin,tisig)
   a%tasigmax = max(a%tasigmax,tisig)
   a%tasigmea = a%tasigmea + tisig
end subroutine


subroutine supm_FinalizeStats(a)
   use typre
   use Mod_ThreeField
   implicit none
   class(ThreeFieldNSProblem) :: a
   
   a%tamea=a%tamea/real(a%nmean,rp)
   a%tasigmea=a%tasigmea/real(a%nmean,rp)
   a%remea=a%remea/real(a%nmean,rp)
   
end subroutine


   

