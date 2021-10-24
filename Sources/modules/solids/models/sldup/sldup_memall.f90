subroutine sldup_memall(a)
   !-----------------------------------------------------------------------
   !****f* SOLIDS/sldup_memall
   ! NAME 
   !    sldup_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_TimeIntegrator
   use Mod_UPSolids
   implicit none

   class(UPSolidsProblem) :: a

   type(TimeIntegratorDt2) :: Integrator
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: nsteps,ndime,npoin,pgaus,nelem,ncomp,ielem,pnode
                  
   !Unknowns      
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   call Integrator%Init(a%kfl_tsche_2nd_datafile)
   call Integrator%GetNumberOfTimeSteps(nsteps)

   a%ncomp = nsteps + 1
   ncomp = a%ncomp

   call a%Memor%alloc(npoin,a%ncomp,        a%press,'press','sldup_memall')

   if(a%kfl_docoupconv) then 
       call a%Memor%alloc(npoin      ,a%press_cp,'press_cp','sldup_memall')
   endif

   !if(a%kfl_printGaussSigma) call a%Memor%alloc(nelem,a%sigma_g,'sigma_g','sldup_memall')

   if(a%kfl_printJ2Stresses)   call a%Memor%alloc(nelem,a%j2_g,'J2','sldup_memall')

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sldup_memall')

   if(a%kfl_printJ2Stresses) then 
       do ielem=1,nelem

       call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
       call a%Mesh%ElementLoad(ielem,e)  

       allocate(a%j2_g(ielem)%a(ndime,pgaus))
       a%j2_g(ielem)%a     = 0.0_rp

       end do
   endif

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sldup_memall')

   call a%MemallExtend

end subroutine sldup_memall

subroutine sldup_memallSpecific(a)
   !-----------------------------------------------------------------------
   !****f* SOLIDS/sldup_memall
   ! NAME 
   !    sldup_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_UPSolids
   implicit none

   class(UPSolidsProblem)          :: a
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ndime,npoin,pgaus,nelem,ncomp,sz,ielem,pnode
   integer(ip) :: usgs_coun,psgs_coun,ncsgs,resid_coun

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   sz = (ndime*(ndime+1))/2

   if(a%kfl_printNodalSigma) call a%Memor%alloc(sz   ,npoin  ,1,a%sigma,'sigma','sldup_memall')
   if(a%kfl_printGaussSigma) call a%Memor%alloc(nelem,a%sigma_g,'sigma_g','sldup_memall')

   !Residual Projection
   if (a%kfl_repro >=1) then
      call a%Memor%alloc(ndime+1,npoin,a%repro,'repro','sldup_memall')
   endif


   if (a%kfl_printGaussSigma) then

       call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_memall')
       do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

         call a%Mesh%ElementLoad(ielem,e)  
         allocate(a%sigma_g(ielem)%a(sz,pgaus))
         a%sigma_g(ielem)%a  = 0.0_rp

       end do

       call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_memall')

   endif

   !Postprocess the residual at the gauss points
   if (a%kfl_printResiduals) then

      resid_coun = 0

      call a%Memor%alloc(nelem,a%residualU,'residual','sldup_memall')
      call a%Memor%alloc(nelem,a%residualP,'residual','sldup_memall')

      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

         allocate(a%residualU(ielem)%a(ndime,pgaus))
         allocate(a%residualP(ielem)%a(pgaus))
         a%residualU(ielem)%a = 0.0_rp
         a%residualP(ielem)%a = 0.0_rp
         resid_coun = resid_coun + (ndime+1)*pgaus

      end do
      call a%Memor%allocObj(0,'residual%a','sldup_memall',resid_coun*rp)
   endif

   !Transient subgrid scales
   if(a%kfl_trasg/=0) then

      usgs_coun = 0
      psgs_coun = 0
      ncsgs     = 2               !Number of components in time for accuracy

      if(a%kfl_tacsg>0) ncsgs=a%ncomp-1

      call a%Memor%alloc(nelem,a%u_sgs,'u_sgs','sldup_memall')
      call a%Memor%alloc(nelem,a%p_sgs,'p_sgs','sldup_memall')

      do ielem=1,nelem

         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%u_sgs(ielem)%a(ndime,ncsgs,pgaus))
         a%u_sgs(ielem)%a = 0.0_rp
         allocate(a%p_sgs(ielem)%a(pgaus))
         a%p_sgs(ielem)%a = 0.0_rp

         usgs_coun = usgs_coun + ndime*ncsgs*pgaus
         psgs_coun = psgs_coun + pgaus

      end do

      call a%Memor%allocObj(0,'u_sgs%a','u_sgs',usgs_coun*rp)
      call a%Memor%allocObj(0,'p_sgs%a','p_sgs',psgs_coun*rp)

   end if

end subroutine sldup_memallSpecific
