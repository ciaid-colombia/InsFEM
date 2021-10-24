module Mod_DoubleIndexMatmul

contains

function DoubleIndexMatmul(elmat,elunk)
   use typre
   implicit none
   real(rp), contiguous, target :: elmat(:,:,:,:), elunk(:,:)
   
   real(rp),allocatable, target :: DoubleIndexMatmul(:,:)
   
   real(rp), pointer :: auxelmat(:,:) => NULL(), auxelunk(:) => NULL(), auxelrhs(:) => NULL()
   integer(ip) :: matsize
   
   allocate(DoubleIndexMatmul(size(elunk,1),size(elunk,2)))
   
   matsize = size(elmat,1)*size(elmat,2)
   
   
   auxelmat(1:matsize,1:matsize) => elmat
   auxelrhs(1:matsize) => DoubleIndexMatmul
   auxelunk(1:matsize) => elunk
   
   auxelrhs = matmul(auxelmat,auxelunk)
   
   
end function

end module