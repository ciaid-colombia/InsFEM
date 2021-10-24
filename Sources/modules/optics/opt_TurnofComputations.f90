subroutine opt_TurnofComputations(a)
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
   
   real(rp) :: ct2
   integer(ip) :: ipoin,npoin
  
   
   
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
   
  
   !---------------------------------------------------------------
   !Compute the cn2 mean value
   a%avg_cn2(:,1) = a%avg_cn2(:,1)/a%istep
   a%avg_cn2u53(:,1) = a%avg_cn2u53(:,1)/a%istep
   
   !Compute average values
   call opt_RayQuadratures(a,a%avg_cn2(:,1),a%avg_cn2u53(:,1),cn2_quadrature,cn2u53_quadrature,length_quadrature,r0,vwind,fg)
   
   !Postprocess
   if (a%MPIrank == a%MPIroot) then
      !do ibeam = 1,a%nbeams
      !   write(*,*) cn2_quadrature(ibeam), cn2u53_quadrature, length_quadrature(ibeam)
      !enddo
      
      do ibeam = 1,a%nbeams
         write(a%lun_outres,*)
         write(a%lun_outres,*) '-------------------------------------------'
         write(a%lun_outres,*) 'Averaged Results for beam number: ',ibeam
         write(a%lun_outres,100) r0(ibeam)
         write(a%lun_outres,101) vwind(ibeam)
         write(a%lun_outres,102) fg(ibeam)
         write(a%lun_outres,103) length_quadrature(ibeam)
         write(a%lun_outres,104) cn2_quadrature(ibeam)
         write(a%lun_outres,105) cn2u53_quadrature(ibeam)
         write(a%lun_outres,*) '-------------------------------------------'
      enddo
   endif
   
   
   
   !---------------------------------------------------------------
   !Compute the Cn2 value from averaged values of p, teta etc
   
   !First we average everything
   a%avg_tdiss = a%avg_tdiss/a%istep
   a%avg_vdiss = a%avg_vdiss/a%istep
   a%avg_tempe = a%avg_tempe/a%istep
   a%avg_press = a%avg_press/a%istep
   a%avg_vnorm = a%avg_vnorm/a%istep

   
   if (a%units == 0) then
      !Pressure from Pascals to milibars + 1 atmosphere (1000 milibars)
      a%avg_press = (a%avg_press+1e5)/100.0_rp
      !Temperature to Kelvin
      a%avg_tempe = a%avg_tempe + 273.15
   endif  
   
   call a%Mesh%GetNpoin(npoin)
   do ipoin = 1,npoin
      !First we compute Ct2
      if (abs(a%avg_vdiss(ipoin)) < 1e-14) then
         ct2 = 0.0_rp
      else
         ct2 = (a%a_Obukhov**2)*a%avg_tdiss(ipoin)*(abs(a%avg_vdiss(ipoin))**(-1.0_rp/3.0_rp))
      endif
      !Secondly cn2
      a%avg_cn2(ipoin,2) = (79e-6*a%avg_press(ipoin)/a%avg_tempe(ipoin)**2)**2*ct2
            
      !Then cn2 x u^%/3
      a%avg_cn2u53(ipoin,2) = a%avg_vnorm(ipoin)**(5.0_rp/3.0_rp)*a%avg_cn2(ipoin,2)
   enddo
   
   !Compute average values
   call opt_RayQuadratures(a,a%avg_cn2(:,2),a%avg_cn2u53(:,2),cn2_quadrature,cn2u53_quadrature,length_quadrature,r0,vwind,fg)
   
   
   !Postprocess
   if (a%MPIrank == a%MPIroot) then
      !do ibeam = 1,a%nbeams
      !   write(*,*) cn2_quadrature(ibeam), cn2u53_quadrature, length_quadrature(ibeam)
      !enddo
       
      do ibeam = 1,a%nbeams
         write(a%lun_outres,*)
         write(a%lun_outres,*) '-------------------------------------------'
         write(a%lun_outres,*) 'Results from averaged quantities for beam number: ',ibeam
         write(a%lun_outres,100) r0(ibeam)
         write(a%lun_outres,101) vwind(ibeam)
         write(a%lun_outres,102) fg(ibeam)
         write(a%lun_outres,103) length_quadrature(ibeam)
         write(a%lun_outres,104) cn2_quadrature(ibeam)
         write(a%lun_outres,105) cn2u53_quadrature(ibeam)
         write(a%lun_outres,*) '-------------------------------------------'
      enddo
   endif
  
   100 format('   r0                : ',e12.6)
   101 format('   vwind             : ',e12.6)
   102 format('   fg                : ',e12.6)
   103 format('   length            : ',e12.6)
   104 format('   cn2_quadrature    : ',e12.6)
   105 format('   cn2u53_quadrature : ',e12.6)
   
end subroutine
