subroutine opt_EndsteComputations(a)
   use typre
   use Mod_Optics
   use Mod_Element
   use MPI
   use def_parame
   use Mod_int2str
   use Mod_iofile
   implicit none
   class(OpticsProblem) :: a
   
   integer(ip) :: ierr

   integer(ip) :: ibeam,ielem,nelem
   real(rp)    :: intersection(3,2)
   logical     :: kfl_intersect
   
   integer(ip) :: iinter
   
   real(rp), allocatable :: elcn2(:),elcn2u53(:)
   real(rp) :: gpcn2(1),gpcn2u53(1)
   
   integer(ip) :: npoinLocal,inode,istat
   real(rp)    :: weighting
   
   real(rp) :: length
   real(rp) :: length_vec(3)
   real(rp) :: xloc(3)
   
   real(rp) :: cn2_quadrature(a%nbeams), length_quadrature(a%nbeams),cn2u53_quadrature(a%nbeams)
   real(rp) :: r0(a%nbeams), fg(a%nbeams),vwind(a%nbeams)
  
   
   
   interface
      subroutine opt_RayQuadratures(a,cn2,cn2u53,cn2_quadratureb,cn2u53_quadratureb,length_quadratureb,r0,vwind,fg)
         use typre
         use Mod_Optics
         implicit none
         class(OpticsProblem) :: a
         real(rp) :: cn2(*), cn2u53(*)
         real(rp) :: cn2_quadratureb(a%nbeams), length_quadratureb(a%nbeams),cn2u53_quadratureb(a%nbeams)
         real(rp) :: r0(a%nbeams), fg(a%nbeams),vwind(a%nbeams)
      end subroutine   
   end interface
   
   !Compute quadrature values
   call opt_RayQuadratures(a,a%cn2(:),a%cn2u53(:),cn2_quadrature,cn2u53_quadrature,length_quadrature,r0,vwind,fg)
   
   !Postprocess
   if (a%MPIrank == a%MPIroot) then
      if (a%istep == 1) then
         write(a%lun_outres,101)  (ibeam,ibeam = 1,a%nbeams)
         write(a%lun_outres,103,advance='no') 'Step    ',('r0           vwind        fg           cn2_quad     cn2u53_quad  ',ibeam = 1,a%nbeams)
         write(a%lun_outres,104) ' '
      endif
   
      write(a%lun_outres,100) a%istep, (r0(ibeam), vwind(ibeam), fg(ibeam), cn2_quadrature(ibeam),cn2u53_quadrature(ibeam), ibeam=1,a%nbeams)
      if (a%kfl_flush == 1) call flush(a%lun_outres)
   endif
   
   
   
   

   100 format(i6,3x,200(e12.6,x))
   101 format('Beam Nº: ',i1,i65,100(i65))
   103 format(a8,100(a70))
   104 format(a1)

end subroutine