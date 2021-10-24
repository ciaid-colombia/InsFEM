module Mod_HeightGauge
   use typre
   use Mod_CutMesh
   use MOd_MPIObject
   use MPI
   implicit none
   
   type, extends(MPIObject) :: HeightGauge 
      integer(ip) :: ndime
      real(rp) :: Origin(3), DirectionVector(3)
      real(rp) :: height = -1e12
contains
      procedure :: ResetGauge
      procedure :: SetNdime
      procedure :: CheckHeight
      procedure :: CollectHeight
   end type

contains
   
   
   subroutine ResetGauge(a)
      class(HeightGauge) :: a
      
      a%height = -1e12
   end subroutine
   
   subroutine SetNdime(a,ndime)
      class(HeightGauge) :: a
      integer(ip) :: ndime
      
      a%ndime = ndime
   end subroutine
   
   subroutine CheckHeight(a,ninters,xglob)
      class(HeightGauge) :: a
      integer(ip) :: ninters
      real(rp) :: xglob(a%ndime,ninters)
      
      real(rp) :: auxxglob(3,3),distance
      logical :: isfound
      
      real(rp) :: lambda2
      
      if (a%ndime == 3) then
         if (ninters == 3) then
         
            call MollerTrumbore(a%DirectionVector,a%Origin,xglob,isfound,distance)
            if (isfound) a%height = distance
         
         elseif (ninters == 4) then
            !We just check the 4 possible triangles
            !Not efficient (twice the real effort) but only in the surface elements
            
            auxxglob(1:3,1) = xglob(1:3,1)
            auxxglob(1:3,2) = xglob(1:3,2)
            auxxglob(1:3,3) = xglob(1:3,3)
            call MollerTrumbore(a%DirectionVector,a%Origin,auxxglob,isfound,distance)
            if (isfound) a%height = distance
            
            auxxglob(1:3,1) = xglob(1:3,1)
            auxxglob(1:3,2) = xglob(1:3,2)
            auxxglob(1:3,3) = xglob(1:3,4)
            call MollerTrumbore(a%DirectionVector,a%Origin,auxxglob,isfound,distance)
            if (isfound) a%height = distance
            
            auxxglob(1:3,1) = xglob(1:3,2)
            auxxglob(1:3,2) = xglob(1:3,3)
            auxxglob(1:3,3) = xglob(1:3,4)
            call MollerTrumbore(a%DirectionVector,a%Origin,auxxglob,isfound,distance)
            if (isfound) a%height = distance
            
            auxxglob(1:3,1) = xglob(1:3,1)
            auxxglob(1:3,2) = xglob(1:3,3)
            auxxglob(1:3,3) = xglob(1:3,4)
            call MollerTrumbore(a%DirectionVector,a%Origin,auxxglob,isfound,distance)
            if (isfound) a%height = distance
         else
            call runend('CheckHeight: Wrong number of intersections')
         endif
      
         
      else
         if (ninters == 2) then
            lambda2 = ( xglob(2,2)*a%DirectionVector(1) - xglob(1,1)*a%DirectionVector(2) + a%Origin(1)*a%DirectionVector(2) - a%Origin(2)*a%DirectionVector(1))/ &
            (xglob(2,1)*a%DirectionVector(1) - xglob(2,2)*a%DirectionVector(1) - xglob(1,1)*a%DirectionVector(2) + xglob(1,2)*a%DirectionVector(2))
            
            if (lambda2 >= -1e-6 .and. lambda2 <= 1.000001) then
               if (a%DirectionVector(1) /= 0.0_rp) then
                  distance = (xglob(1,1)-a%Origin(1) - xglob(1,1)*lambda2 + xglob(1,2)*lambda2) / a%DirectionVector(1)
               else
                  distance = (xglob(2,1)-a%Origin(2) - xglob(2,1)*lambda2 + xglob(2,2)*lambda2) / a%DirectionVector(2)
               endif   
               a%Height = distance
            endif
            
         else
            call runend('CheckHeight: Wrong number of intersections')
         endif
      endif

   
   end subroutine
   
   subroutine MollerTrumbore(DirectionVector,Origin,xglob,isfound,distance)
      implicit none
      real(rp) :: DirectionVector(3), Origin(3), xglob(3,3)
      logical  :: isfound
      real(rp) :: distance
      
      real(rp) :: sysmat(3,3), sysrhs(3), syssol(3)
      
      isfound = .false.
      !Moller Trumbore algorithm
      sysmat(1:3,1) = -DirectionVector(1:3)
      sysmat(1:3,2) = xglob(1:3,2)-xglob(1:3,1)
      sysmat(1:3,3) = xglob(1:3,3)-xglob(1:3,1)
      
      sysrhs(1:3) = Origin(1:3)-xglob(1:3,1)
      call SolveSystem(3,sysmat,sysrhs,syssol)
      
      if ((syssol(2) >= 0) .and.(syssol(3) >= 0) .and. (syssol(2)+syssol(3) <= 1) ) then
         isfound = .true.
         distance = syssol(1)
      endif
   end subroutine
            
   
   subroutine CollectHeight(a)
      class(HeightGauge) :: a
      
      integer(ip) :: ierr
      real(rp) :: height
      
      call MPI_REDUCE(a%height, height, 1, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
      if (a%MPIrank == a%MPIroot) a%height = height
      
      
   end subroutine
   
   
end module