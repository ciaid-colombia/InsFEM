module Mod_timer
   use typre
   use MPI
   implicit none
   private
   public :: timer,cputim

   type Timer
      private
      real(rp) :: TimeCount = 0.0_rp
      real(rp) :: TimeInitial
      real(rp) :: TimeLast = 0.0_rp
      
contains
     procedure :: ToZero
     procedure :: Tic
     procedure :: Toc
     procedure :: Tic2
     procedure :: Toc2
     procedure :: GetValue
     procedure :: GetLast
   end type
  
contains

   subroutine ToZero(a)
      implicit none
      class(timer) :: a
   
      a%TimeCount = 0.0_rp
   end subroutine
   
   subroutine Tic2(a)
      implicit none
      class(timer) :: a
      
      !Several timers need to be added (openmp)
      call cputim2(a%TimeInitial)
   end subroutine
   
   subroutine Tic(a)
      implicit none
      class(timer) :: a
      
      !Several timers need to be added (openmp)
      call cputim(a%TimeInitial)
   end subroutine
   
   subroutine Toc2(a)
      implicit none
      class(timer) :: a
      real(rp) :: tfin
      
      call cputim2(tfin)
      a%TimeLast  = tfin - a%TimeInitial
      a%TimeCount = a%TimeCount + a%TimeLast
   end subroutine
   
   subroutine Toc(a)
      implicit none
      class(timer) :: a
      real(rp) :: tfin
      
      call cputim(tfin)
      a%TimeLast  = tfin - a%TimeInitial
      a%TimeCount = a%TimeCount + a%TimeLast
   end subroutine
   
   subroutine GetValue(a,TimeCount)
      implicit none
      class(timer) :: a
      real(rp)    :: TimeCount
      
      TimeCount = a%TimeCount
   end subroutine
   
   subroutine GetLast(a,TimeLast)
      implicit none
      class(timer) :: a
      real(rp)    :: TimeLast
      
      TimeLast = a%TimeLast
   end subroutine
   
   subroutine cputim(time)
      implicit none
      real(rp) :: time
      time = MPI_WTIME()
   end subroutine
   
   subroutine cputim2(time)
      implicit none
      real(rp) :: time
      call cpu_time(time)
   end subroutine
   
end module
