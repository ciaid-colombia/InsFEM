subroutine plcd_PostprocessMaterialData(b,iiter_char)
   use Mod_plcd_BaseElmope
   use Mod_PLCD
   use Mod_plcd_BMatrix
   use Mod_plcd_BMatrixFactory 
   implicit none
   class(PLCDProblem), target :: b
   character(5) :: iiter_char
   
   integer(ip) :: imaterial
   
   a=>b
   
   
   !Setup Postprocess
   do imaterial = 1,a%NumberOfMaterials
      call a%Materials(imaterial)%p%SetupPostprocess(a%Mesh,a%Memor)
   enddo
   
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_EnditeElmope')
   call a%Mesh%GetNelem(nelem)
   elements : do ielem = 1,nelem
      ElementMatData => a%ElementMaterialsData(ielem)%p

      call ElementMatData%ContributeToPostprocess(ielem)
      

   enddo elements

   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Elmope')

   !DoPostprocess
   do imaterial = 1,a%NumberOfMaterials
      call a%Materials(imaterial)%p%DoPostprocess(a%istep,a%ctime,a%FilePostpr,a%Mesh,a%Memor,iiter_char)
   enddo
   
end subroutine
