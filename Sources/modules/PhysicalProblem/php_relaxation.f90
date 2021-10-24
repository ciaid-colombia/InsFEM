subroutine php_relaxVector(a,itask,w_r,w_rmax,r_i1,r_i21)
   !This routine solves an iteration of the physical problem equations
   use MPI
   use typre
   use Mod_PhysicalProblem
   implicit none
   class(PhysicalProblem) :: a
   character(5) :: itask
   real(rp) :: vnor,prod,part_prod,w,w_r,w_rmax
   integer  idime,ndime,inod,npoin,npoinLocal,ierr
   real(rp), intent(in) :: r_i1(:,:),r_i21(:,:) 

   w = 0.0_rp

   select case(itask)

   case('AITKE')

      !Calculating Aitken relaxation variables, see kuttler2008 for notations

      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNpoinLocal(npoinLocal)     
      call a%Mesh%GetNdime(ndime)

      !So that Aitken relaxes
      !r_i+2(:,:)  = x_guess(i+2,:)   - x_(i+1,:)

      !r_i+1(:,:)  = x_guess(i+1,:)   - x_(i,:)

      !r_i21(:,:) = r_i+2     - r_i+1

      part_prod = 0.0_rp
      prod = 0.0_rp
      do inod = 1, npoin
          do idime =1,ndime
              part_prod = part_prod + r_i1(idime,inod)*r_i21(idime,inod)
          end do
      end do

      call MPI_REDUCE(part_prod,prod,1,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)

      call vecnorGenericMPI(r_i21,npoinLocal,ndime,vnor,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)

      if (a%MPIrank == a%MPIroot) then

          if(vnor > 0.0_rp) then
              if(a%cpiter .le. a%kfl_initRelax) then !We allow two iters to have residuals
                  w_r = min(w_r,w_rmax)
              else
                  w_r = -w_r*(prod/(vnor**2))
              endif
          endif

      endif

      CALL MPI_BCAST(w_r,1,MPI_REAL8,a%MPIroot,a%MPIcomm,ierr)

   end select

end subroutine php_relaxVector


