subroutine sld_turnofGeneral(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)          :: ndime,npoin,ncomp
   integer(ip)          :: ielem,nelem,pgaus,pnode
   integer(ip)          :: nboun,sz
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNelem(nelem)
   ncomp = a%ncomp

   sz = (ndime*(ndime+1))/2

   
   !Write tail for formatted files
   if (a%MPIrank == a%MPIroot ) then
       if (a%kfl_isSlave .eqv. .false.) then
           write(a%lun_outph,'(//,a,/)') '     * * * END OF ELASTIC SOLID RUN * * *'
           write(a%lun_solve,'(//,a,/)') '     * * * END OF ELASTIC SOLID RUN * * *'
       endif
   endif
   
   call a%Memor%dealloc(ndime,npoin,a%ncomp,a%disp,     'disp','sld_turnof')
   if(a%kfl_NodalStress) call a%Memor%dealloc(sz   ,npoin  ,a%stress,'stress','sld_turnof')
   if(a%kfl_NodalStrain) call a%Memor%dealloc(sz   ,npoin  ,a%strain,'strain','sld_turnof')

   call a%Memor%dealloc(ndime,npoin,a%btraction_nodal,'Bnodaltraction','sld_turnof')

   if(a%kfl_docoupconv) then 
       call a%Memor%dealloc(ndime,npoin,a%disp_cp, 'disp_cp','sld_turnof')
   endif

   if (a%kfl_timei==1) then
       call a%Memor%dealloc(ndime,npoin,2,a%veloc,'veloc','sld_turnof')
       call a%Memor%dealloc(ndime,npoin,2,a%accel,'accel','sld_turnof')
   endif

   call a%Memor%dealloc(npoin,a%kfl_fixrs,'kfl_fixrs','sld_turnof')
   call a%Memor%dealloc(nboun,a%kfl_bours,'kfl_bours','sld_turnof')

   if (a%kfl_sourcesPresent) then
       call a%Memor%dealloc(a%pwSourceDim,a%pwSourceSize,a%pwSource,'pwSource','sld_turnof')
       call a%Memor%dealloc(a%pwSourceSize,a%pwSourceId,'pwSourceId','sld_turnof')
   endif      

   if(a%kfl_isSlave .eqv. .false.) then
      if (a%nptra/=0) then
         if (a%MPIrank == a%MPIroot) then
            call iofile(two,a%lun_trap,a%fil_trap,'SLD TRACKING OF POINTS')
         endif
      endif
   endif

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
       deallocate(a%extForce(ielem)%a)
       if(a%kfl_printStress) deallocate(a%stress_g(ielem)%a)
       if(a%kfl_printStrain) deallocate(a%strain_g(ielem)%a)
       deallocate(a%btraction(ielem)%a)
   end do

   call a%Memor%dealloc(nelem,a%extForce,'extForce','sld_turnof')
   if(a%kfl_printStress) call a%Memor%dealloc(nelem,a%stress_g,'stress_g','sld_turnof')
   if(a%kfl_printStrain) call a%Memor%dealloc(nelem,a%strain_g,'strain_g','sld_turnof')
   call a%Memor%dealloc(nelem,a%btraction,'btraction','sld_turnof')

   call a%SolidSpecificTurnof
  
end subroutine sld_turnofGeneral

subroutine sld_turnofSpecific(a)
   ! DESCRIPTION
   !    This routine closes the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)          :: ndime,ielem,nelem,npoin
   integer(ip)          :: sz,pgaus,pnode
   
   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   
   sz = (ndime*(ndime+1))/2

   call a%Memor%dealloc(sz,sz,a%c_elastic,   'c_elastic','sld_turnof')

   if(a%kfl_printNodalSigma) call a%Memor%dealloc(sz   ,npoin,1,a%sigma,'sigma','sld_memall')
   if(a%kfl_printNodalPress) call a%Memor%dealloc(      npoin,1,a%press,'press','sld_memall')

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
       if (a%kfl_printGaussSigma) then
           deallocate(a%sigma_g(ielem)%a)
       endif

       if (a%kfl_printGaussPress) then
           deallocate(a%press_g(ielem)%a)
       endif

       !Memory error here, beware
       deallocate(a%pStrain(ielem)%a)

       if(a%kfl_printPrincipalStresses) then
           deallocate(a%printPStrain(ielem)%a)
       end if

   end do

   if(a%kfl_printPrincipalStresses) then
       call a%Memor%dealloc(nelem,a%printPStrain,'printPStrain','sld_turnof')
   end if
   call a%Memor%dealloc(nelem,a%pStrain,'pStrain','sld_turnof')
   if (a%kfl_printGaussSigma) call a%Memor%dealloc(nelem,a%sigma_g,'sigma_g','sld_turnof')
   if (a%kfl_printGaussPress) call a%Memor%dealloc(nelem,a%press_g,'press_g','sld_turnof')

end subroutine sld_turnofSpecific
