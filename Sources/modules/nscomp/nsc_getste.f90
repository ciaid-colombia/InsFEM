subroutine nsc_getste(a,dtinv)
!-----------------------------------------------------------------------
! NAME 
!    nsc_getste
! DESCRIPTION
!    This routine computes the time step size for the Compressible NS
!    equation.
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_NSCompressibleSubroutines
   use Mod_Mesh
   use Mod_Element
   use Mod_GatherScatterDtcri
   implicit none 
   class(NSCompressibleProblem) :: a
   real(rp) :: dtinv
   
   real(rp), allocatable :: elmom(:,:)
   real(rp), allocatable :: elden(:),elene(:)
   real(rp), allocatable :: pcmom(:),pcden(:),pcene(:)

   class(FiniteElement), pointer :: e => NULL()
   real(rp) :: pctem,pcspd,pcvno,dtmin,dtcri,hclen,rdinv
   integer(ip) :: ielem,nelem,ndime
   real(rp)    :: acvis,actco,accph,accvh
   
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
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_getste')
     
      call a%Memor%alloc(e%mnode,elden,'elden','nsc_getste')
      call a%Memor%alloc(ndime,e%mnode,elmom,'elmom','nsc_getste')
      call a%Memor%alloc(e%mnode,elene,'elene','nsc_getste')
      call a%Memor%alloc(1_ip,pcden,'pcden','nsc_getste')
      call a%Memor%alloc(ndime,pcmom,'pcmom','nsc_getste')
      call a%Memor%alloc(1_ip,pcene,'pcene','nsc_getste')
      
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
         call e%gather(1_ip,elden,a%densf(:,1))  
         call e%gather(ndime,elmom,a%momen(:,:,1))  
         call e%gather(1_ip,elene,a%energ(:,1))  
         call e%interpc(1_ip,elden,pcden)
         call e%interpc(ndime,elmom,pcmom)
         call e%interpc(1_ip,elene,pcene)

         !Advection velocity at center of gravity
         call nsc_ComputeAdvectionVelocity(e,pcden(1),pcmom,pcvno)
           
         !Temperature at center of gravity
         call nsc_ComputeTemperature(a,e,pcden(1),pcmom,pcene(1),pctem)

         !Speed of sound at center of gravity
         call nsc_ComputeSoundSpeed(accph,accvh,pctem,pcspd)

         dtcri=dtcri+(pcvno+pcspd)/hclen

         dtcri=1.0_rp/dtcri

         dtmin=min(dtmin,dtcri)

      end do
      a%dtcri = dtmin
      
      !Gather Dtcri from all processes, find minumum dtcri value, scatter it
      call php_GatherScatterDtcri(a)
      
      a%dtinv = 1.0_rp/(a%dtcri*a%safet)
      dtinv = a%dtinv

      !Memory deallocation
      call a%Memor%dealloc(e%mnode,elden,'elden','nsc_getste')
      call a%Memor%dealloc(ndime,e%mnode,elmom,'elmom','nsc_getste')
      call a%Memor%dealloc(e%mnode,elene,'elene','nsc_getste')
      call a%Memor%dealloc(1_ip,pcden,'pcden','nsc_getste')
      call a%Memor%dealloc(ndime,pcmom,'pcmom','nsc_getste')
      call a%Memor%dealloc(1_ip,pcene,'pcene','nsc_getste')
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','nsc_getste')
   end if
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
end subroutine nsc_getste
