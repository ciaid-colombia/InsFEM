subroutine tem_getste(a,dtinv)
   use typre
   use Mod_Temperature
   use Mod_Element
   use Mod_GatherScatterDtcri
   use Mod_nsm_Viscosity
   use MPI
   implicit none
   class(TemperatureProblem) :: a
   real(rp) :: dtinv
   
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: gpvel(:),grvel(:,:)
   class(FiniteElement), pointer :: e => NULL()
   
   real(rp) :: acrcp,acden,acsph,actco,dtcri,dtmin,gptem,gpvno,hclen,rdinv,acrea,acsou,acvis,vista
   integer(ip) :: ielem,nelem,irank,ndime
   
   interface
      subroutine tem_coupling_outerr(a)
         use Mod_Temperature
         implicit none
         class(TemperatureProblem) :: a
      end subroutine
   
   end interface
   
   !Check if the coupling between modules is done
   call tem_coupling_outerr(a)
   
   
   call a%Timer%Total%Tic
   call a%Timer%Getste%Tic

   !Critical time step
   if(a%kfl_timei/=0.and.a%kfl_stead/=1) then
   
      !Dimensions and general variables
      call a%Mesh%GetNdime(ndime) 
      call a%Mesh%GetNelem(nelem)
      rdinv=1.0_rp/real(ndime)
      gpvno = 0.0_rp
      dtmin = 1e6
      
      !Element Initialization
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','temgetste')
      
      if (a%kfl_advec >= 1) then
         call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','tem_getste')
         call a%Memor%alloc(ndime,gpvel,'gpvel','tem_getste')
      endif
      
      !Smagorinsky or Wale
      if (a%kfl_cotur == -1 .or. a%kfl_cotur == -2) then
         call a%Memor%alloc(ndime,ndime,grvel,'grvel','tem_getste')
      endif
      
      do ielem = 1,nelem
         call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
         
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)
         
         !Derivative (and detjm) at the center of gravity
         call e%elmdcg
         
         !Element length
         call e%elmlen
         
         !Element length
         hclen=(e%weicg*e%detjm)**rdinv
         
         !Critical Time step computation
         dtcri=0.0_rp
         
         !Advection exists
         if (a%kfl_advec>=1) then
            if (a%kfl_advec == 1) then
               call e%gather(ndime,elvel,a%veloc(:,:))  
            endif
            call e%interpc(ndime,elvel,gpvel)
            call vecnor(gpvel,e%ndime,gpvno,2)
            
            !2u/h
            dtcri=dtcri+2.0_rp*gpvno/hclen
         endif
         
         acvis = actco/acsph
         !Smagorinsky
         if (a%kfl_cotur < 0) then
            call e%gradient(e%ndime,elvel(:,:),grvel)    !Vel. gradient
            if (a%kfl_cotur == -1) then
               call nsm_smago(e,grvel,acden,a%turbu,vista)
            elseif (a%kfl_cotur == -2) then
               call nsm_wale(e,grvel,acden,a%turbu,vista)
            endif
            
            !Prandtl turbulent number
            vista = vista*a%prtur
            acvis = acvis + vista
         endif
         
         !4*(k/rho*cp)/h^2
         dtcri=dtcri+4.0_rp*(acvis/acden)/(hclen*hclen)
         
         dtcri=1.0_rp/dtcri
         dtmin=min(dtmin,dtcri)
      end do
      a%dtcri = dtmin
      
      !Gather Dtcri from all processes, find minimum and scatter
      call php_GatherScatterDtcri(a)
      
      a%dtinv = 1.0_rp/(a%dtcri*a%safet)
      dtinv = a%dtinv
      
      !Memory deallocation
      if (a%kfl_advec >= 1) then
         call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','tem_getste')
         call a%Memor%dealloc(ndime,gpvel,'gpvel','tem_getste')
      endif
      !Smagorinsky
      if (a%kfl_cotur == -1 .or. a%kfl_cotur == -2) then
         call a%Memor%dealloc(ndime,ndime,grvel,'grvel','tem_getste')
      endif
      
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','tem_getste')
   end if
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
   
end subroutine
