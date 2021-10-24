subroutine nsc_pr_getste(a,dtinv)
!-----------------------------------------------------------------------
! NAME 
!    nsc_getste
! DESCRIPTION
!    This routine computes the time step size for the Compressible NS
!    equation.
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressiblePrimitive
   use Mod_NSCompressibleSubroutines
   use Mod_Mesh
   use Mod_Element
   use Mod_GatherScatterDtcri
   implicit none 
   class(NSCompressiblePrimitiveProblem) :: a
   real(rp) :: dtinv
   
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: eltem(:)
   real(rp), allocatable :: pcvel(:),pctem(:)

   class(FiniteElement), pointer :: e => NULL()
   real(rp) :: pcspd,pcvno,dtmin,dtcri,hclen,rdinv
   integer(ip) :: ielem,nelem,ndime
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: vstar
   
   call a%Timer%Total%Tic
   call a%Timer%Getste%Tic

   !Critical time step
   if(a%kfl_timei/=0.and.a%kfl_stead/=1) then
   
      !Dimensions and general variables
      call a%Mesh%GetNdime(ndime) 
      call a%Mesh%GetNelem(nelem)
      rdinv=1.0_rp/real(ndime)
      pcvno = 0.0_rp
      dtmin = 1e6

      !Physical Parameters
      call a%GetPhysicalParameters(acvis,actco,accph,accvh)

      !Element Initialization
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
     
      call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','nsc_getste')
      call a%Memor%alloc(e%mnode,eltem,'elten','nsc_getste')
      call a%Memor%alloc(ndime,pcvel,'pcvel','nsc_getste')
      call a%Memor%alloc(1_ip,pctem,'pctem','nsc_getste')
      
      do ielem = 1,nelem
         
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)
         
         !Derivative (and detjm) at the center of gravity
         call e%elmdcg
         
         !Element length
         hclen=(e%weicg*e%detjm)**rdinv

         !Critical Time step computation
         dtcri=0.0_rp
         
         !Values at the center of gravity
         call e%gather(ndime,elvel,a%veloc(:,:,1))  
         call e%gather(1_ip,eltem,a%tempe(:,1))  
         call e%interpc(ndime,elvel,pcvel)
         call e%interpc(1_ip,eltem,pctem)

         !Velocity norm
         call vecnor(pcvel,e%ndime,pcvno,2)
           
         !Speed of sound at center of gravity
         call nsc_ComputeSoundSpeed(accph,accvh,(pctem(1)+a%reltem),pcspd)
      
         !Characteristic nscomp velocity
         call a%ComputeNscompNstincVelocity(pcvno,pcspd,vstar)
         
         dtcri=dtcri+vstar/hclen

         if(dtcri>epsilon(0.0_rp)) then
            dtcri=1.0_rp/dtcri
         else
            dtcri=1.0e12_rp
         end if

         dtmin=min(dtmin,dtcri)

      end do
      a%dtcri = dtmin
      
      !Gather Dtcri from all processes, find minumum dtcri value, scatter it
      call php_GatherScatterDtcri(a)
      
      a%dtinv = 1.0_rp/(a%dtcri*a%safet)
      dtinv = a%dtinv

      !Memory deallocation
      call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','nsc_getste')
      call a%Memor%dealloc(e%mnode,eltem,'elten','nsc_getste')
      call a%Memor%dealloc(ndime,pcvel,'pcvel','nsc_getste')
      call a%Memor%dealloc(1_ip,pctem,'pctem','nsc_getste')
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   end if
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
end subroutine nsc_pr_getste
