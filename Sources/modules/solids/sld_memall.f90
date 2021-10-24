subroutine sld_memallGeneral(a)
   !-----------------------------------------------------------------------
   !****f* SOLIDS/sld_memall
   ! NAME 
   !    sld_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_TimeIntegrator
   use Mod_Solids
   implicit none

   class(SolidsProblem)          :: a
   type(TimeIntegratorDt2)       :: Integrator
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: nsteps,ndime,npoin,pgaus,nelem,ncomp,sz,ielem,pnode

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   sz = (ndime*(ndime+1))/2

   call Integrator%Init(a%kfl_tsche_2nd_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)

   a%ncomp = nsteps + 1
   ncomp = a%ncomp

   call a%Memor%alloc(ndime,npoin,ncomp,a%disp, 'disp','sld_memall')

   if(a%kfl_NodalStress) call a%Memor%alloc(sz   ,npoin  ,a%stress,'stress','sldsup_memall')
   if(a%kfl_NodalStrain) call a%Memor%alloc(sz   ,npoin  ,a%strain,'strain','sldsup_memall')

   call a%Memor%alloc(ndime,npoin,a%btraction_nodal, 'Bnodaltraction','sld_memall')

   if(a%kfl_docoupconv) then 
       call a%Memor%alloc(ndime,npoin,a%disp_cp, 'disp_cp','sld_memall')
   endif

   if (a%kfl_timei==1) then
       call a%Memor%alloc(ndime,npoin,2,a%veloc,'veloc','sld_memall')
       call a%Memor%alloc(ndime,npoin,2,a%accel,'accel','sld_memall')
   end if

   call a%Memor%alloc(nelem,a%extForce,'extForce','sld_memall')
   if(a%kfl_printStress) call a%Memor%alloc(nelem,a%stress_g,'stress_g','sld_memall')
   if(a%kfl_printStrain) call a%Memor%alloc(nelem,a%strain_g,'strain_g','sld_memall')
   call a%Memor%alloc(nelem,a%btraction,'btraction','sld_memall')

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_memall')

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

       call a%Mesh%ElementLoad(ielem,e)  
       allocate(a%extForce(ielem)%a(a%ndofn,e%mnode))

       if(a%kfl_printStress) then
           allocate(a%stress_g(ielem)%a(sz,pgaus))
           a%stress_g(ielem)%a  = 0.0_rp
       endif

       if(a%kfl_printStrain) then
           allocate(a%strain_g(ielem)%a(sz,pgaus))
           a%strain_g(ielem)%a  = 0.0_rp
       endif

       allocate(a%btraction(ielem)%a(ndime,pgaus))

       a%extForce(ielem)%a  = 0.0_rp
       a%btraction(ielem)%a = 0.0_rp

   end do

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_memall')

   call a%SolidSpecificMemall

end subroutine sld_memallGeneral


subroutine sld_memallSpecific(a)
   !-----------------------------------------------------------------------
   !****f* SOLIDS/sld_memall
   ! NAME 
   !    sld_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_TimeIntegrator
   use Mod_Solids
   implicit none

   class(SolidsProblem)          :: a
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: nsteps,ndime,npoin,pgaus,nelem,sz,ielem,pnode

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   sz = (ndime*(ndime+1))/2

   if(a%kfl_printNodalSigma) call a%Memor%alloc(sz   ,npoin,1,a%sigma,'sigma','sld_memall')
   if(a%kfl_printNodalPress) call a%Memor%alloc(      npoin,1,a%press,'press','sld_memall')

   if(a%kfl_printGaussSigma) call a%Memor%alloc(nelem,a%sigma_g,'sigma_g','sldup_memall')
   if(a%kfl_printGaussPress) call a%Memor%alloc(nelem,a%press_g,'press_g','sldup_memall')

   call a%Memor%alloc(sz,sz,a%c_elastic,   'c_elastic','sld_memall')

   call a%Memor%alloc(nelem,a%pStrain,'pStrain','sld_memall')

   if(a%kfl_printPrincipalStresses) then
       call a%Memor%alloc(nelem,a%printPStrain,'printPStrain','sld_memall')
   end if

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_memall')

   do ielem=1,nelem
       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

       call a%Mesh%ElementLoad(ielem,e)  


       allocate(a%pStrain(ielem)%a(ndime))
       a%pStrain(ielem)%a   = 0.0_rp

       if(a%kfl_printPrincipalStresses) then
           allocate(a%printPStrain(ielem)%a(sz,pgaus))
           a%printPStrain(ielem)%a = 0.0_rp
       end if

       if (a%kfl_printGaussPress) then
         allocate(a%press_g(ielem)%a(pgaus))
         a%press_g(ielem)%a  = 0.0_rp
       endif

       if (a%kfl_printGaussSigma) then
           allocate(a%sigma_g(ielem)%a(sz,pgaus))
           a%sigma_g(ielem)%a  = 0.0_rp
       endif

   end do

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_memall')

end subroutine sld_memallSpecific

