module Mod_NaivePartition
   use typre
   implicit none
   private
   public NaivePartition
   
   type :: NaivePartition
      integer(ip) :: npoinLocal
      integer(ip) :: gnpoin
      integer(ip) :: MPIsize
   
contains
      procedure :: Initialize
      procedure :: GetIniEndPoint
      procedure :: GetPointRank
   end type
   
contains   
   subroutine Initialize(a,gnpoin,MPIsize)
      use typre
      implicit none
      class(NaivePartition) :: a
      integer(ip) :: irank,MPIsize,gnpoin
      
      a%MPIsize = MPIsize
      a%gnpoin = gnpoin
      
      a%npoinLocal = ceiling(real(a%gnpoin)/real(a%MPIsize))
   
   end subroutine
   
   !Used by many subroutines
   subroutine getiniendpoint(a,irank,inipoin,endpoin)
      use typre
      implicit none
      class(NaivePartition) :: a
      integer(ip) :: irank,inipoin,endpoin,npoin

      inipoin = a%npoinLocal*(irank)+1
      if (irank < a%MPIsize-1) then
         endpoin = a%npoinLocal*(irank+1)
      else
         endpoin = a%gnpoin
      endif
      
      if (inipoin > a%gnpoin) then
         inipoin = a%gnpoin+1
      endif
      if (endpoin > a%gnpoin) then
         endpoin = a%gnpoin
      endif
      
   end subroutine

   subroutine getPointRank(a,ipoin,irank)
      use typre
      implicit none
      class(NaivePartition) :: a
      integer(ip) :: ipoin,irank

      irank = (ipoin-1)/(a%npoinLocal)
      if (irank == a%MPIsize) irank = a%MPIsize-1
   
   end subroutine

end module
