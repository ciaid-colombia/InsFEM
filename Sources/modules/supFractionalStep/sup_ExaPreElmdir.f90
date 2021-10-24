subroutine supf_ExaPreElmdir(a,e,elmat,elrhs)
   use typre
   use Mod_Mesh
   use Mod_SupExacso       
   use Mod_ThreeField
   use Mod_SUPFractionalStep
   use Mod_Element
   implicit none 
   class(SUPFractionalStepProblem) :: a
   class(FiniteElement) :: e
   type(SupExacso)  :: exacso      
   
   real(rp) :: elmat(1,e%mnode,1,e%mnode), elrhs(1,e%mnode)
     
   real(rp) :: adiag
   integer(ip) :: inode,ipoin,jnode
   !Exact Values
   real(rp)                   :: exprg(e%ndime)
   real(rp)                   :: expre
   real(rp), pointer          :: coord(:) => NULL()   
   
   
   real(rp)    :: prepr
     
   do inode = 1,e%pnode
      ipoin = e%lnods(inode)
      
      call a%Mesh%GetPointCoord(ipoin,coord)
      
      !Exact solution
      call exacso%sup_ComputeSolution(e%ndime,coord,a%ctime,a%LogFormulation,a)
      call exacso%sup_GetPressure(e%ndime,expre,exprg)        

      prepr = expre

      adiag=elmat(1,inode,1,inode)
      elmat(1,inode,1,1:e%pnode)=0.0_rp
      do jnode = 1,e%pnode
         elrhs(1,jnode) = elrhs(1,jnode)- elmat(1,jnode,1,inode)*prepr
      enddo
      elmat(1,1:e%pnode,1,inode) = 0.0_rp
      elmat(1,inode,1,inode) = adiag
      elrhs(1,inode) = adiag*prepr
      
   enddo
end subroutine 
