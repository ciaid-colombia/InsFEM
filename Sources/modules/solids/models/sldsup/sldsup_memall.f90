subroutine sldsup_memall(a)
   !-----------------------------------------------------------------------
   !****f* SOLIDS/sldsup_memall
   ! NAME 
   !    sldsup_memall
   ! DESCRIPTION
   !    This routine allocates memory for the arrays needed to solve the problem
   !-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_TimeIntegrator
   use Mod_SUPSolids
   implicit none

   class(SUPSolidsProblem) :: a

   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: nsteps,ndime,npoin,pgaus,nelem,sz,ielem,pnode
                  
   !Unknowns      
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   sz = (ndime*(ndime+1))/2

   call a%Memor%alloc(sz   ,npoin,a%ncomp,a%sigma,'sigma','sldsup_memall')

   call a%Memor%alloc(sz   ,npoin  ,a%devstrain,'devstrain','sldsup_memall')

   if(a%kfl_docoupconv) then 
       call a%Memor%alloc(sz   ,npoin,a%sigma_cp,'sigma_cp','sldsup_memall')
   endif

   if(a%kfl_printJ2Stresses) then 
       call a%Memor%alloc(npoin  ,a%j2,'J2','sldsup_memall')
   endif

   call a%MemallExtend2

end subroutine sldsup_memall

subroutine sldsup_memallSpecific(a)
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
   use Mod_SUPSolids
   implicit none

   class(SUPSolidsProblem)          :: a
   class(FiniteElement), pointer :: e => NULL()
   integer(ip) :: ndime,npoin,pgaus,nelem,ncomp,sz,ielem,pnode
   integer(ip) :: usgs_coun,psgs_coun,ssgs_coun,ncsgs,resid_coun

   !Unknowns
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)

   sz = (ndime*(ndime+1))/2

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_memall')

   !Residual Projection
   if (a%kfl_repro >=1) then
      call a%Memor%alloc(ndime+sz+1,npoin,a%repro,'repro','sldsup_memall')
   endif

   !Postprocess the residual at the gauss points
   if (a%kfl_printResiduals) then

      resid_coun = 0

      call a%Memor%alloc(nelem,a%residualU,'residual','sldsup_memall')
      call a%Memor%alloc(nelem,a%residualP,'residual','sldsup_memall')
      call a%Memor%alloc(nelem,a%residualS,'residual','sldsup_memall')

      do ielem=1,nelem
         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)

         allocate(a%residualU(ielem)%a(ndime,pgaus))
         allocate(a%residualP(ielem)%a(pgaus))
         allocate(a%residualS(ielem)%a(sz,pgaus))
         a%residualU(ielem)%a = 0.0_rp
         a%residualP(ielem)%a = 0.0_rp
         a%residualS(ielem)%a = 0.0_rp
         resid_coun = resid_coun + (ndime+sz+1)*pgaus

      end do
      call a%Memor%allocObj(0,'residual%a','sldsup_memall',resid_coun*rp)
   endif

   !Transient subgrid scales
   if(a%kfl_trasg/=0) then

      usgs_coun = 0
      ssgs_coun = 0
      psgs_coun = 0
      ncsgs     = 2               !Number of components in time for accuracy

      if(a%kfl_tacsg>0) ncsgs=a%ncomp-1

      call a%Memor%alloc(nelem,a%u_sgs,'u_sgs','sldsup_memall')
      call a%Memor%alloc(nelem,a%s_sgs,'s_sgs','sldsup_memall')
      call a%Memor%alloc(nelem,a%p_sgs,'p_sgs','sldsup_memall')

      do ielem=1,nelem

         call a%Mesh%GetElemArraySize(ielem,pnode,pgaus)
         allocate(a%u_sgs(ielem)%a(ndime,ncsgs,pgaus))
         a%u_sgs(ielem)%a = 0.0_rp
         allocate(a%s_sgs(ielem)%a(sz,pgaus))
         a%s_sgs(ielem)%a = 0.0_rp
         allocate(a%p_sgs(ielem)%a(pgaus))
         a%p_sgs(ielem)%a = 0.0_rp

         usgs_coun = usgs_coun + ndime*ncsgs*pgaus
         ssgs_coun = ssgs_coun + sz*pgaus
         psgs_coun = psgs_coun + pgaus

      end do

      call a%Memor%allocObj(0,'u_sgs%a','u_sgs',usgs_coun*rp)
      call a%Memor%allocObj(0,'s_sgs%a','s_sgs',ssgs_coun*rp)
      call a%Memor%allocObj(0,'p_sgs%a','p_sgs',psgs_coun*rp)

   end if

   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_memall')

end subroutine sldsup_memallSpecific

