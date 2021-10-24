
subroutine supm_elmdir(a,e,elmat,elrhs)
   use typre
   use Mod_nsm_elmdir
   use Mod_ThreeField
   use Mod_Element
   use Mod_Mesh
   use Mod_sup_elmdir
   use Mod_php_elmdir
   implicit none
   class(ThreeFieldNSProblem)             :: a
   class(FiniteElement), intent(in)        :: e
   real(rp), intent(inout)                :: elmat(a%ndofn,e%mnode,a%ndofn,e%mnode), elrhs(a%ndofn,e%mnode)
   integer(ip)               :: aux,aux1
   integer(ip) :: currentbvess,bcstar,auxtens,ndime 
   !To do multy materials
   integer(ip) :: imat=1
      
         
   call nsm_pressdir(a,e,a%ndofn,elmat,elrhs)
   call nsm_rotdir(a,e,a%ndofn,elmat,elrhs)
   call a%Mesh%Getndime(ndime)
   auxtens=(ndime-1)*(ndime-1)+2  
      
   !------------------------------------------------------------------------
   if(a%kfl_bc_number>0.and.a%kfl_confi==1.and.a%MatProp(imat)%lawvi<0)then
      call sup_Exaelmdir(a,e,a%ndofn,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
   else
      call php_elmdir(a,e,a%ndofn,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
   end if   
       
   call nsm_signdir(a,e,a%ndofn,elmat,elrhs)
end subroutine      
         



