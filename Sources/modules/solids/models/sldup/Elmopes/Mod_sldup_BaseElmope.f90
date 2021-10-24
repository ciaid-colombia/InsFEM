submodule (Mod_sld_BaseElmope) sldup_BaseElmope

   contains

   module subroutine dynamicResidualUP
      implicit none
      real(rp) :: massAccel(e%ndime,e%pnode)
      
      call matvecmult(e%ndime,e%pnode,massmat,elaccel(:,:,1),massAccel)

      elrhu = elrhu - massAccel/(2.0_rp*up%mu) 
      
   end subroutine

   module subroutine calculateAU_U_dyn
      implicit none
      integer(ip) :: nd
      integer(ip) :: u1,uf,s1,sf,p1,bc
      real(rp)    :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call getAU_U_dyn(nd,densi,gpaccel,grdpre,divstr,AU(u1:uf))

   end subroutine calculateAU_U_dyn

   module subroutine calculateAU_U
      implicit none
      integer(ip) :: nd
      integer(ip) :: u1,uf,s1,sf,p1,bc
      real(rp)    :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call getAU_U(nd,grdpre,divstr,AU(u1:uf))

   end subroutine calculateAU_U

   module subroutine calculateAU_P
      implicit none
      integer(ip) :: nd,tn
      integer(ip) :: u1,uf,s1,sf,p1,bc
      real(rp)    :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call getAU_P(nd,G,lam,det,B,gppress,AU(p1))

   end subroutine calculateAU_P

   module subroutine pressGather
      implicit none
      integer(ip) :: nd

      call a%Mesh%GetNdime(nd)
      
      call e%gather(1,elpress(1,:),up%press(:,1))
   end subroutine  

   module subroutine calculatePressGradient
       implicit none
       integer(ip) :: nd,tn,inode,idime,jdime
       real(rp)    :: istr(e%ndime,e%ndime,e%mnode)

       call a%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       divdisp = 0.0_rp
       grdpre  = 0.0_rp

       !Pressure gradient
       call e%gradient(1,elpress(1,:),grdpre)

       call e%divergence(eldisp(:,:,1),divdisp)

   end subroutine calculatePressGradient

   module subroutine calculateGPStress
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: b_ii
      real(rp), allocatable :: Bt(:),i(:)

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call a%Memor%alloc(tn,Bt,'Bt','sldup_Elmope')
      call a%Memor%alloc(tn,i,'i','sldup_Elmope')

      call get2ndIITensor_voigt(nd,tn,i)

      call trace(nd,B,b_ii)
      call getVoigtStress(tn,nd,B,Bt)

      gpsigma = 0.0_rp

      !Stress tensor term: (b-(b_ii/nd)*ii)
      gpsigma = (a%mu/det)*(Bt - (b_ii/nd)*i)

      call a%Memor%dealloc(tn,Bt,'Bt','sldup_Elmope')
      call a%Memor%dealloc(tn,i,'i','sldup_Elmope')

   end subroutine calculateGPStress

   module subroutine calculateGPPress
      implicit none
      integer(ip) :: nd,tn
      real(rp)    :: b_ii
      real(rp), allocatable :: Bt(:),i(:)

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call trace(nd,B,b_ii)

      gppress= 0.0_rp

      !Stress tensor term: (b-(b_ii/nd)*ii)
      gppress= (1.0_rp/det)*((a%lambda*log(det)-a%mu) + ((a%mu*b_ii)/nd))

   end subroutine calculateGPPress

   module subroutine calculateStressDivergence2D
      implicit none
      integer(ip) :: nd,inode
      real(rp)    :: i(e%ndime,e%ndime),b_ii
      real(rp)    :: dx,dy,dz

      call a%Mesh%GetNdime(nd)

      call get2ndIITensor(nd,i)
      call trace(nd,B,b_ii)

      divstr = 0.0_rp

      do inode=1,e%pnode

          !Testing against shape function
          dx = e%cartd(1,inode)
          dy = e%cartd(2,inode)

          divstr(1) = divstr(1) + dx*B(1,1) + dy*B(1,2) - dx*(b_ii/nd)
          divstr(2) = divstr(2) + dx*B(2,1) + dy*B(2,2) - dy*(b_ii/nd)

      end do

      divstr=(up%mu/det)*divstr

   end subroutine calculateStressDivergence2D

   module subroutine calculateStressDivergence3D
      implicit none
      integer(ip) :: nd,inode
      real(rp)    :: i(e%ndime,e%ndime),b_ii
      real(rp)    :: dx,dy,dz

      call a%Mesh%GetNdime(nd)

      call get2ndIITensor(nd,i)
      call trace(nd,B,b_ii)

      divstr = 0.0_rp

      do inode=1,e%pnode

          !Testing against shape function
          dx = e%cartd(1,inode)
          dy = e%cartd(2,inode)
          dz = e%cartd(3,inode)

          divstr(1) = divstr(1) + dx*B(1,1) + dy*B(1,2) + dz*B(1,3) - dx*(b_ii/nd)
          divstr(2) = divstr(2) + dx*B(2,1) + dy*B(2,2) + dz*B(2,3) - dy*(b_ii/nd)
          divstr(3) = divstr(3) + dx*B(3,1) + dy*B(3,2) + dz*B(3,3) - dz*(b_ii/nd)

      end do

      divstr=(up%mu/det)*divstr

   end subroutine calculateStressDivergence3D

   module subroutine InterpolateGpPress
      implicit none

      call e%interpg(1_ip,elpress(1,:),gppress(1))
   end subroutine

   module subroutine calculateSigmaPress
       implicit none
       integer(ip) :: kfl_nonlinear

       call ProcHook%PhysicalProp
       call calculateGradientsAndDeter
       call calculateGPStress
       call calculateGPPress

   end subroutine calculateSigmaPress

   module subroutine calculateUPGradients
       implicit none
       integer(ip) :: kfl_nonlinear

       call ProcHook%PhysicalProp
       call calculateGradientsAndDeter
       call calculatePressGradient
       call calculateGPStress

       if (e%ndime.eq.2) then
           call calculateStressDivergence2D
       else
           call calculateStressDivergence3D
       endif

   end subroutine calculateUPGradients

   module subroutine calculateFmatGradient
       implicit none
       integer(ip) :: nd

       call a%Mesh%GetNdime(nd)

       grdFmat = 0.0_rp
       grdFspt = 0.0_rp

       !u = x_sp - X_mat; grad(u) = ii - F(x)^-1
       !d/dX(grad(u)) = d/dx(F(x))
       !d/dx(grad(u)) = d/dx(-F(x)^-1)
       call e%hessian(nd, eldisp(:,:,1),grdFmat)  
       call e%hessian(nd,-eldisp(:,:,1),grdFspt)  

       !calculate div(Fmat) div(Fspat)
       if (nd.eq.2) then

           divFmat(1) = grdFmat(1,1) + grdFmat(1,2)
           divFmat(2) = grdFmat(2,1) + grdFmat(2,2)

           divFspt(1) = grdFspt(1,1) + grdFspt(2,1)
           divFspt(2) = grdFspt(1,2) + grdFspt(2,2)

       elseif (nd.eq.3) then

           divFmat(1) = grdFmat(1,1) + grdFmat(1,2) + grdFmat(1,3)
           divFmat(2) = grdFmat(2,1) + grdFmat(2,2) + grdFmat(2,3)
           divFmat(3) = grdFmat(3,1) + grdFmat(3,2) + grdFmat(3,3)

           divFspt(1) = grdFspt(1,1) + grdFspt(2,1) + grdFspt(3,1)
           divFspt(2) = grdFspt(1,2) + grdFspt(2,2) + grdFspt(3,2)
           divFspt(3) = grdFspt(1,3) + grdFspt(2,3) + grdFspt(3,3)

       endif

   end subroutine calculateFmatGradient

   module subroutine getTauParametersUP(tau_u,tau_s,tau_p)
      implicit none
      real(rp), intent(inout):: tau_u,tau_s,tau_p
      integer(ip)            :: pn,nd,tn,mn,idime
      real(rp)               :: hl,freq,tau_ut,tau_ux

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode
      pn  = e%pnode

      hl=e%hleng(1)
      do idime=2,nd
         hl=hl*e%hleng(idime)
      end do
      hl=hl**(1.0_rp/real(nd))

      freq = (up%tau_u*2.0_rp*up%mu)/(hl*hl)

      if(freq>zeroc) then
         tau_ux= 1.0_rp/freq
      else
         tau_ux= 1.0e12_rp
      end if

      tau_ut  = tau_ux

      if(up%kfl_trasg==1) then
        tau_ut = 1/(1/tau_ux + densi*a%dtinv2*tsch_deriv)
        tautau = (1.0_rp/tau_ux)*tau_ut
      endif

      tau_u  = tau_ut
      tau_p = up%tau_p*up%mu*2.0_rp

   end subroutine getTauParametersUP

  module subroutine AssembleMassMatrixUP
      implicit none
      integer(ip) :: nd,pn
      integer(ip) :: u1,uf,s1,sf,p1,bc

      call a%Mesh%GetNdime(nd)
      pn  = e%pnode

      call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      massmat=massmat/(2.0_rp*up%mu)
      !Assembly of jacobian A(d): mass 
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn) + massmat(1:nd,1:pn,1:nd,1:pn)

  end subroutine AssembleMassMatrixUP

   module subroutine AssembleStiffnessMatrixUP
      implicit none
      integer(ip) :: pn,nd,tn,u1,uf,s1,sf,p1,bc

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      pn  = e%pnode

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      ! Momentum equation
      ! Assembly elmuv to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn)  &
                                   + elmuv(1:nd,1:pn,1:nd,1:pn)
      ! Assembly elmpv to elmat
      elmat(u1:uf,1:pn,p1,1:pn)    = elmat(u1:uf,1:pn,p1,1:pn) &
                                   + elmpv(1:nd,1:pn,1,1:pn)

      ! Continuity Equation
      ! Assembly elmuq to elmat
      elmat(p1,1:pn,u1:uf,1:pn)    = elmat(p1,1:pn,u1:uf,1:pn) &
                                   + elmuq(1,1:pn,1:nd,1:pn)
      ! Assembly elmpq to elmat
      elmat(p1,1:pn,p1,1:pn)       = elmat(p1,1:pn,p1,1:pn) &
                                   + elmpq(1,1:pn,1,1:pn)

   end subroutine AssembleStiffnessMatrixUP

   module subroutine AssembleRhsUP
      implicit none
      integer(ip) :: pn,nd,u1,uf,s1,sf,tn,p1,bc

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      ! Assembly elrhu to elrhs
      elrhs(u1:uf,1:pn) = elrhs(u1:uf,1:pn) + elrhu(1:nd,1:pn)
      elrhs(p1,1:pn)    = elrhs(p1,1:pn)    + elrhp(1,1:pn)   

   end subroutine AssembleRhsUP

   module subroutine AssembleRhsForcesUP
      implicit none
      integer(ip) :: pn,nd,u1,uf,s1,sf,tn,p1,bc

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      ! Assembly elrhu to elrhs
      elrhs(u1:uf,1:pn) = elrhs(u1:uf,1:pn) + a%extForce(ielem)%a(u1:uf,:) 
      elrhs(p1,1:pn)    = elrhs(p1,1:pn)    + a%extForce(ielem)%a(p1,:)    

   end subroutine AssembleRhsForcesUP

   module subroutine AssembleStiffnessMatrixUP_dynsgs
      implicit none
      integer(ip)          :: pn,nd,tn,u1,uf,s1,sf,p1,bc
      real(rp)             :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2
      pn = e%pnode

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

      !-------------TU------------------
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn)  &
                                   - elm_usgsdyn_VU(1:nd,1:pn,1:nd,1:pn)
      ! Assembly elmvp to elmat
      elmat(u1:uf,1:pn,p1,1:pn)    = elmat(u1:uf,1:pn,p1,1:pn) &
                                   - elm_usgsdyn_VP(1:nd,1:pn,1,1:pn)

   end subroutine AssembleStiffnessMatrixUP_dynsgs

   module subroutine AssembleStiffnessMatrixUP_SGS
      implicit none
      integer(ip)          :: pn,nd,tn,u1,uf,s1,sf,p1,bc
      real(rp)             :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2
      pn = e%pnode

      call up%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

      !-------------TU------------------
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn)  &
                                   - tau_u*elm_usgs_VU(1:nd,1:pn,1:nd,1:pn)
      ! Assembly elmvp to elmat
      elmat(u1:uf,1:pn,p1,1:pn)    = elmat(u1:uf,1:pn,p1,1:pn) &
                                   - tau_u*elm_usgs_VP(1:nd,1:pn,1,1:pn)
      ! Assembly elmqu to elmat
      elmat(p1,1:pn,u1:uf,1:pn)    = elmat(p1,1:pn,u1:uf,1:pn)  &
                                   - tau_u*elm_usgs_QU(1,1:pn,1:nd,1:pn)
      ! Assembly elmqp to elmat
      elmat(p1,1:pn,p1,1:pn)       = elmat(p1,1:pn,p1,1:pn) &
                                   - tau_u*elm_usgs_QP(1,1:pn,1,1:pn)

      !-------------TP------------------
      ! Momentum equation
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn) &
                                   - (tau_p/(2.0_rp*up%mu))*elm_psgs_VU(1:nd,1:pn,1:nd,1:pn)
      ! Assembly elmvp to elmat
      elmat(u1:uf,1:pn,p1,1:pn)    = elmat(u1:uf,1:pn,p1,1:pn) &
                                   - (tau_p/(2.0_rp*up%mu))*elm_psgs_VP(1:nd,1:pn,1,1:pn)
      ! Continuity Equation
      ! Assembly elmqu to elmat
      elmat(p1,1:pn,u1:uf,1:pn)    = elmat(p1,1:pn,u1:uf,1:pn) &
                                   - tau_p*elm_psgs_QU(1,1:pn,1:nd,1:pn)
      ! Assembly elmqp to elmat
      elmat(p1,1:pn,p1,1:pn)       = elmat(p1,1:pn,p1,1:pn) &
                                   - tau_p*elm_psgs_QP(1,1:pn,1,1:pn)

      !---------------------------------END SGS----------------------------
   end subroutine AssembleStiffnessMatrixUP_SGS

   module subroutine ResetUPBaseArrays
      implicit none

      !-----------------Assembly-----------------
      elrhs       = 0.0_rp
      elmat       = 0.0_rp

      !-------------------RHS---------------
      elrhu       = 0.0_rp
      elrhp       = 0.0_rp

      !-------------------u---------------
      elmuv       = 0.0_rp
      elmuq       = 0.0_rp
      !-------------------p---------------
      elmpv       = 0.0_rp
      elmpq       = 0.0_rp

   end subroutine ResetUPBaseArrays

   module subroutine ResetUPBaseArrays_dynsgs
      implicit none
      !-------------------Stabilization---------------
      !-------------------Tau_u---------------
      elm_usgsdyn_VU = 0.0_rp
      elm_usgsdyn_VP = 0.0_rp

   end subroutine ResetUPBaseArrays_dynsgs


   module subroutine ResetUPBaseArrays_SGS
      implicit none
      !-------------------Stabilization---------------
      !-------------------Tau_u---------------
      elm_usgs_VU = 0.0_rp
      elm_usgs_VP = 0.0_rp
      elm_usgs_QU = 0.0_rp
      elm_usgs_QP = 0.0_rp
      !-------------------Tau_p---------------
      elm_psgs_VU = 0.0_rp
      elm_psgs_VP = 0.0_rp
      elm_psgs_QU = 0.0_rp
      elm_psgs_QP = 0.0_rp

   end subroutine ResetUPBaseArrays_SGS

   module subroutine AllocateBaseElmopeMatricesUP_dynsgs
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call up%Memor%alloc(nd,mn,nd,mn,elm_usgsdyn_VU,'elm_usgsdyn_VU','sldup_elmope')
      call up%Memor%alloc(nd,mn,1 ,mn,elm_usgsdyn_VP,'elm_usgsdyn_VP','sldup_elmope')

   end subroutine AllocateBaseElmopeMatricesUP_dynsgs

   module subroutine DeallocateBaseElmopeMatricesUP_dynsgs
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call up%Memor%dealloc(nd,mn,nd,mn,elm_usgsdyn_VU,'elm_usgsdyn_VU','sldup_elmope')
      call up%Memor%dealloc(nd,mn,1 ,mn,elm_usgsdyn_VP,'elm_usgsdyn_VP','sldup_elmope')

   end subroutine DeallocateBaseElmopeMatricesUP_dynsgs

   module subroutine AllocateBaseElmopeMatricesUP_SGS
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call up%Memor%alloc(nd,mn,nd,mn,elm_usgs_VU,'elm_usgs_VU','sldup_elmope')
      call up%Memor%alloc(nd,mn,1 ,mn,elm_usgs_VP,'elm_usgs_VP','sldup_elmope')
      !-----------Pressure stabilization
      call up%Memor%alloc(1 ,mn,nd,mn,elm_usgs_QU,'elm_usgs_QS','sldup_elmope')
      call up%Memor%alloc(1 ,mn,1 ,mn,elm_usgs_QP,'elm_usgs_QP','sldup_elmope')

      !-------------------Tau_p matrix calculation----------
      !-----------Displacement stabilization
      call up%Memor%alloc(nd,mn,nd,mn,elm_psgs_VU,'elm_psgs_VU','sldup_elmope')
      call up%Memor%alloc(nd,mn,1 ,mn,elm_psgs_VP,'elm_psgs_VP','sldup_elmope')
      !-----------Pressure stabilization
      call up%Memor%alloc(1 ,mn,nd,mn,elm_psgs_QU,'elm_psgs_QU','sldup_elmope')
      call up%Memor%alloc(1 ,mn,1 ,mn,elm_psgs_QP,'elm_psgs_QP','sldup_elmope')

   end subroutine

   module subroutine DeallocateBaseElmopeMatricesUP_SGS
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call up%Memor%dealloc(nd,mn,nd,mn,elm_usgs_VU,'elm_usgs_ES','sldup_elmope')
      call up%Memor%dealloc(nd,mn,1 ,mn,elm_usgs_VP,'elm_usgs_EP','sldup_elmope')
      !-----------Pressure stabilization
      call up%Memor%dealloc(1 ,mn,nd,mn,elm_usgs_QU,'elm_usgs_QS','sldup_elmope')
      call up%Memor%dealloc(1 ,mn,1 ,mn,elm_usgs_QP,'elm_usgs_QP','sldup_elmope')

      !-------------------Tau_p matrix calculation----------
      !-----------Displacement stabilization
      call up%Memor%dealloc(nd,mn,nd,mn,elm_psgs_VU,'elm_psgs_VU','sldup_elmope')
      call up%Memor%dealloc(nd,mn,1 ,mn,elm_psgs_VP,'elm_psgs_VP','sldup_elmope')
      !-----------Pressure stabilization
      call up%Memor%dealloc(1 ,mn,nd,mn,elm_psgs_QU,'elm_psgs_QU','sldup_elmope')
      call up%Memor%dealloc(1 ,mn,1 ,mn,elm_psgs_QP,'elm_psgs_QP','sldup_elmope')

   end subroutine

   module subroutine AllocateUPAssemblyMatrices
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      call up%Memor%alloc(1,mn,       elrhp,'elrhp','sldup_Elmope')
      !-----------Displacement testing
      call up%Memor%alloc(nd,mn,nd,mn,elmuv,'elmuv','sldup_elmope')
      call up%Memor%alloc(nd,mn,1 ,mn,elmpv,'elmpv','sldup_elmope')
      !-----------Pressure testing
      call up%Memor%alloc(1 ,mn,nd,mn,elmuq,'elmuq','sldup_elmope')
      call up%Memor%alloc(1 ,mn,1 ,mn,elmpq,'elmpq','sldup_elmope')

   end subroutine

   module subroutine DeallocateUPAssemblyMatrices
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      call up%Memor%dealloc(1,mn,       elrhp,'elrhp','sldup_Elmope')
      !-----------Displacement testing
      call up%Memor%dealloc(nd,mn,nd,mn,elmuv,'elmuv','sldup_elmope')
      call up%Memor%dealloc(nd,mn,1 ,mn,elmpv,'elmpv','sldup_elmope')
      !-----------Pressure testing
      call up%Memor%dealloc(1 ,mn,nd,mn,elmuq,'elmuq','sldup_elmope')
      call up%Memor%dealloc(1 ,mn,1 ,mn,elmpq,'elmpq','sldup_elmope')

   end subroutine

   module subroutine AllocateBaseElmopeMatricesUP
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      call up%Memor%alloc(1 ,mn,elpress,'elpress','sldup_Elmope')
      call up%Memor%alloc(nd,   grdpre,'grdpre','sldup_Elmope')

   end subroutine

   module subroutine DeallocateBaseElmopeMatricesUP
      implicit none
      integer(ip) :: mn,nd

      call a%Mesh%GetNdime(nd)
      mn  = e%mnode

      call up%Memor%dealloc(1 ,mn,elpress,'elpress','sldup_Elmope')
      call up%Memor%dealloc(nd,   grdpre,'grdpre','sldup_Elmope')
   end subroutine

   module subroutine AllocateBaseElmopeMatricesSUP
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      mn  = e%mnode

      call a%Memor%alloc(tn,mn,elsigma,'elsigma','sldup_Elmope')
      call a%Memor%alloc(tn,   gpsigma,'gpsigma','sldup_Elmope')
      call a%Memor%alloc(1 ,mn,elpress,'elpress','sldup_Elmope')
            
      call a%Memor%alloc(nd,   divstr,'divstr','sldup_Elmope')
      call a%Memor%alloc(nd,tn,grsig ,'grsig' ,'sldup_Elmope')
      call a%Memor%alloc(nd,   grdpre,'grdpre','sldup_Elmope')

   end subroutine

   module subroutine DeallocateBaseElmopeMatricesSUP
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      call a%Memor%dealloc(tn,mn,elsigma,'elsigma','sldup_Elmope')
      call a%Memor%dealloc(tn,   gpsigma,'gpsigma','sldup_Elmope')
      call a%Memor%dealloc(1 ,mn,elpress,'elpress','sldup_Elmope')

      call a%Memor%dealloc(nd,   divstr,'divstr','sldup_Elmope')
      call a%Memor%dealloc(nd,tn,grsig ,'grsig' ,'sldup_Elmope')
      call a%Memor%dealloc(nd,   grdpre,'grdpre','sldup_Elmope')
   end subroutine

   module subroutine AllocateBaseElmopeMatricesSUPNonLinear
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      mn  = e%mnode

      call up%Memor%alloc(nd,nd,      F_mat   ,'F_mat'   ,'sldup_Elmope')
      call up%Memor%alloc(nd,nd,      F_matt  ,'F_matt'  ,'sldup_Elmope')
      call up%Memor%alloc(nd,nd,      F_spat  ,'F_spat'  ,'sldup_Elmope')
      call up%Memor%alloc(nd,nd,      B       ,'B'       ,'sldup_Elmope')
      call up%Memor%alloc(a%ndofn,    AU      ,'AU'      ,'sldup_Elmope')
   end subroutine

   module subroutine DeallocateBaseElmopeMatricesSUPNonLinear
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      call up%Memor%dealloc(nd,nd,      F_mat   ,'F_mat'   ,'sldup_Elmope')
      call up%Memor%dealloc(nd,nd,      F_matt  ,'F_matt'  ,'sldup_Elmope')
      call up%Memor%dealloc(nd,nd,      F_spat  ,'F_spat'  ,'sldup_Elmope')
      call up%Memor%dealloc(nd,nd,      B       ,'B'       ,'sldup_Elmope')
      call up%Memor%dealloc(a%ndofn,    AU      ,'AU'      ,'sldup_Elmope')
   end subroutine

  module subroutine AllocateHigherOrderDerivativeArrays
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call a%Memor%alloc(nd,divFmat,'divFmat','sldup_Elmope')
      call a%Memor%alloc(nd,divFspt,'divFspt','sldup_Elmope')
      call a%Memor%alloc(nd,tn,grdFmat,'grdFmat','sldup_Elmope')
      call a%Memor%alloc(nd,tn,grdFspt,'grdFspt','sldup_Elmope')
      call a%Memor%alloc(nd,derJFki,'derJFki','sldup_Elmope')

  end subroutine AllocateHigherOrderDerivativeArrays

  module subroutine DeallocateHigherOrderDerivativeArrays
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2

      call a%Memor%dealloc(nd,divFmat,'divFmat','sldup_Elmope')
      call a%Memor%dealloc(nd,divFspt,'divFspt','sldup_Elmope')
      call a%Memor%dealloc(nd,tn,grdFmat,'grdFmat','sldup_Elmope')
      call a%Memor%dealloc(nd,tn,grdFspt,'grdFspt','sldup_Elmope')
      call a%Memor%dealloc(nd,derJFki,'derJFki','sldup_Elmope')

  end subroutine DeallocateHigherOrderDerivativeArrays

end submodule sldup_BaseElmope
