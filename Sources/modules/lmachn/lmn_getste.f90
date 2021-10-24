subroutine lmn_getste(a,dtinv)
!-----------------------------------------------------------------------
! NAME 
!    lmn_getste
! DESCRIPTION
!    This routine computes the time step size for the incompressible NS
!    equation.
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_LowMach
   use Mod_GatherScatterDtcri
   implicit none 
   class(LowMachProblem) :: a
   real(rp) :: dtinv
   real(rp), allocatable :: grvel(:,:)
   real(rp), allocatable :: elvel(:,:), elden(:)
   real(rp), allocatable :: gpvel(:), gpden(:)
   class(FiniteElement), pointer :: e => NULL()
   real(rp) :: acvit,acvis,acden,gpvno,dtmin,dtcri,hclen,rdinv
   
   integer(ip) :: ielem,nelem,ndime,irank
   
   call a%Timer%Total%Tic
   call a%Timer%Getste%Tic

   !Critical time step
   if(a%kfl_timei/=0.and.a%kfl_stead/=1) then
   
      !Dimensions and general variables
      call a%Mesh%GetNdime(ndime) 
      call a%Mesh%GetNelem(nelem)
      rdinv=1.0_rp/real(ndime)
      gpvno = 0.0_rp
      acvit = 0.0_rp
      dtmin = 1e6
      
      !Physical Parameters
      call a%GetPhysicalParameters(acvis=acvis)
      
      
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lmn_getste')
      
      call a%Memor%alloc(ndime,ndime,grvel,'grvel','lmn_getste')
      call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','lmn_getste')
      call a%Memor%alloc(e%mnode,elden,'elden','lmn_getste')
      call a%Memor%alloc(ndime,gpvel,'gpvel','lmn_getste')
      call a%Memor%alloc(1,gpden,'gpden','lmn_getste')
      
      do ielem = 1,nelem
         
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)
         
         !Derivative (and detjm) at the center of gravity
         call e%elmdcg
         
         !Element length
         hclen=(e%weicg*e%detjm)**rdinv
         
         !Critical Time step computation
         dtcri=0.0_rp

         !Interpolate density at the center of gravity
         acden = a%densi

         !Advection exists
         if (a%kfl_advec==1) then
            !Norm of the velocity at the center of gravity
            call e%gather(ndime,elvel,a%veloc(:,:,1))  
            call e%interpc(ndime,elvel,gpvel)
            !Velocity norm
            call vecnor(gpvel,e%ndime,gpvno,2)
            
            !2u/h
            dtcri=dtcri+2.0_rp*gpvno/hclen
         
         endif
         
         !Smagorinsky, modify viscosity
!         if (a%kfl_cotur == -1) then
!            call e%gradient(ndime,elvel,grvel)
!            call nsm_smago(e,grvel,acden,a%turbu,acvit)
!         endif
!         acvis=a%visco+acvit
         
!          4*(mu/rho)/h^2
         dtcri=dtcri+4.0_rp*(acvis/acden)/(hclen*hclen)
         
         dtcri=1.0_rp/dtcri
         dtmin=min(dtmin,dtcri)
      end do
      a%dtcri = dtmin
      
      !Gather Dtcri from all processes, find minumum dtcri value, scatter it
      call php_GatherScatterDtcri(a)
      
      a%dtinv = 1.0_rp/(a%dtcri*a%safet)
      dtinv = a%dtinv
      
      !Memory deallocation
      call a%Memor%dealloc(ndime,ndime,grvel,'grvel','lmn_getste')
      call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','lmn_getste')
      call a%Memor%dealloc(e%mnode,elden,'elden','lmn_getste')
      call a%Memor%dealloc(ndime,gpvel,'gpvel','lmn_getste')
      call a%Memor%dealloc(1,gpden,'gpden','lmn_getste')
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lmn_getste')

   end if
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
end subroutine lmn_getste
