submodule (Mod_sld_BaseElmope) sldsup_BaseElmope

   contains

   !Gathers and interpolates
   module subroutine sigmaGather
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call e%gather(tn,elsigma(:,:),sup%sigma(:,:,1))
   end subroutine  
   
   module subroutine InterpolateGpSigma
      implicit none
      integer(ip) :: nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call e%interpg(tn,elsigma(:,:),gpsigma(:))
   end subroutine

   module subroutine calculatePressGradientsAndStressDiv
       implicit none
       integer(ip) :: nd,tn,inode,idime,jdime
       real(rp)    :: istr(e%ndime,e%ndime,e%mnode)

       call sup_NH%Mesh%GetNdime(nd)
       tn=(nd*(nd+1))/2

       divdisp = 0.0_rp
       grdpre  = 0.0_rp
       divstr  = 0.0_rp
       grsig   = 0.0_rp

       !Pressure gradient
       call e%gradient(1,elpress(1,:),grdpre)

       !Stress divergence 
       istr = 0.0_rp
       do inode=1,e%pnode
           call getStressTensor(tn,nd,elsigma(:,inode),istr(:,:,inode))
       end do
       do idime=1,e%ndime
           call e%divergence(istr(idime,:,:),divstr(idime))
       end do

       call e%divergence(eldisp(:,:,1),divdisp)

       call e%gradient(tn,elsigma,grsig)  


   end subroutine calculatePressGradientsAndStressDiv

   module subroutine calculateGradients
       implicit none
       integer(ip) :: kfl_nonlinear

       call ProcHook%PhysicalProp
       call calculateGradientsAndDeter
       call calculatePressGradientsAndStressDiv

   end subroutine calculateGradients

   module subroutine AssembleStiffnessMatrixSUP
      implicit none
      integer(ip) :: pn,nd,tn,u1,uf,s1,sf,p1,bc

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      pn  = e%pnode

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      ! Momentum equation
      ! Assembly elmsv to elmat
      elmat(u1:uf,1:pn,s1:sf,1:pn) = elmat(u1:uf,1:pn,s1:sf,1:pn)  &
                                   + elmsv(1:nd,1:pn,1:tn,1:pn)
      ! Assembly elmpv to elmat
      elmat(u1:uf,1:pn,p1,1:pn)    = elmat(u1:uf,1:pn,p1,1:pn) &
                                   + elmpv(1:nd,1:pn,1,1:pn)

      !Constitutive equation
      ! Assembly elmut to elmat
      elmat(s1:sf,1:pn,u1:uf,1:pn) = elmat(s1:sf,1:pn,u1:uf,1:pn) &
                                   + elmut(1:tn,1:pn,1:nd,1:pn)
      ! Assembly elmst to elmat
      elmat(s1:sf,1:pn,s1:sf,1:pn) = elmat(s1:sf,1:pn,s1:sf,1:pn) &
                                   + elmst(1:tn,1:pn,1:tn,1:pn)
   
      ! Continuity Equation
      ! Assembly elmuq to elmat
      elmat(p1,1:pn,u1:uf,1:pn)    = elmat(p1,1:pn,u1:uf,1:pn) &
                                   + elmuq(1,1:pn,1:nd,1:pn)
      ! Assembly elmpq to elmat
      elmat(p1,1:pn,p1,1:pn)       = elmat(p1,1:pn,p1,1:pn) &
                                   + elmpq(1,1:pn,1,1:pn)

   end subroutine AssembleStiffnessMatrixSUP

   module subroutine AssembleRhsSUP
      implicit none
      integer(ip) :: pn,nd,u1,uf,s1,sf,tn,p1,bc

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      ! Assembly elrhu to elrhs
      elrhs(u1:uf,1:pn) = elrhs(u1:uf,1:pn) + elrhu(1:nd,1:pn)
      elrhs(s1:sf,1:pn) = elrhs(s1:sf,1:pn) + elrhd(1:tn,1:pn) 
      elrhs(p1,1:pn)    = elrhs(p1,1:pn)    + elrhp(1,1:pn)   

   end subroutine AssembleRhsSUP

   module subroutine AssembleRhsForcesSUP
      implicit none
      integer(ip) :: pn,nd,u1,uf,s1,sf,tn,p1,bc

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      pn  = e%pnode

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      ! Assembly elrhu to elrhs
      elrhs(u1:uf,1:pn) = elrhs(u1:uf,1:pn) + a%extForce(ielem)%a(u1:uf,:) 
      elrhs(s1:sf,1:pn) = elrhs(s1:sf,1:pn) + a%extForce(ielem)%a(s1:sf,:) 
      elrhs(p1,1:pn)    = elrhs(p1,1:pn)    + a%extForce(ielem)%a(p1,:)    

   end subroutine AssembleRhsForcesSUP

   module subroutine getTauParameters_linear(tau_u,tau_s,tau_p)
      implicit none
      real(rp), intent(inout):: tau_u,tau_s,tau_p
      integer(ip)            :: pn,nd,tn,mn
      real(rp)               :: G2,cl,hl

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode
      pn  = e%pnode

      cl  = sup%kfl_PrbCharL
      hl  = e%hleng(nd)

      !Use secant shear modulus 
      !It spoils convergence, so not used generally
      if(sup%kfl_useSecantModulus) then
          if(sup%devnorm .le. zeroc) then
              G2  = 2.0_rp*sup%mu
          else
              G2 = sup%signorm/sup%devnorm
          endif
      else
          G2  = 2.0_rp*sup%mu
      endif

      tau_u = sup%tau_u*((hl*hl)/(G2)) 
      tau_s = sup%tau_s*(hl/cl)
      tau_p = sup%tau_p*(hl/cl) 

   end subroutine getTauParameters_linear

   module subroutine getTauParameters(tau_u,tau_s,tau_p)
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

      freq = (sup%tau_u*sup%mu)/(hl*hl)

      if(freq>zeroc) then
         tau_ux= 1.0_rp/freq
      else
         tau_ux= 1.0e12_rp
      end if

      tau_ut  = tau_ux

      if(sup%kfl_trasg==1) then
        tau_ut = 1/(1/tau_ux + densi*a%dtinv2*tsch_deriv)
        tautau = (1.0_rp/tau_ux)*tau_ut
      endif

      tau_u  = tau_ut

      tau_s = sup%tau_s*sup%mu*2.0_rp
      tau_p = sup%tau_p*sup%mu*2.0_rp

   end subroutine getTauParameters
   
   module subroutine AssembleStiffnessMatrixSUP_SGS_dynsgs
      implicit none
      integer(ip)          :: pn,nd,tn,u1,uf,s1,sf,p1,bc

      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2
      pn = e%pnode

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      !-------------TU------------------
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn)    = elmat(u1:uf,1:pn,u1:uf,1:pn) &
                                   - elm_usgs_VU(1:nd,1:pn,1:nd,1:pn)
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,s1:sf,1:pn)    = elmat(u1:uf,1:pn,s1:sf,1:pn) &
                                   - elm_usgs_VS(1:nd,1:pn,1:tn,1:pn)

      ! Assembly elmqp to elmat
      elmat(u1:uf,1:pn,1,1:pn)       = elmat(u1:uf,1:pn,1,1:pn) &
                                   - elm_usgs_VP(1:nd,1:pn,1,1:pn)

   end subroutine AssembleStiffnessMatrixSUP_SGS_dynsgs

   module subroutine AssembleStiffnessMatrixSUP_SGS_dyn
      implicit none
      integer(ip)          :: pn,nd,tn,u1,uf,s1,sf,p1,bc
      real(rp)             :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2
      pn = e%pnode

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

      !-------------TU------------------
      ! Assembly elmsu to elmat
      elmat(s1:sf,1:pn,u1:uf,1:pn)    = elmat(s1:sf,1:pn,u1:uf,1:pn) &
                                   - tau_u*elm_usgs_EU(1:tn,1:pn,1:nd,1:pn)

      ! Assembly elmqp to elmat
      elmat(p1,1:pn,u1:uf,1:pn)       = elmat(p1,1:pn,u1:uf,1:pn) &
                                   - tau_u*elm_usgs_QU(1,1:pn,1:nd,1:pn)

   end subroutine AssembleStiffnessMatrixSUP_SGS_dyn

   module subroutine AssembleStiffnessMatrixSUP_SGS
      implicit none
      integer(ip)          :: pn,nd,tn,u1,uf,s1,sf,p1,bc
      real(rp)             :: tau_u,tau_s,tau_p

      call a%Mesh%GetNdime(nd)
      tn = (nd*(nd+1))/2
      pn = e%pnode

      call sup%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      call ProcPointer%getTauParameters(tau_u,tau_s,tau_p)

      !-------------TU------------------
      ! Assembly elmes to elmat
      elmat(s1:sf,1:pn,s1:sf,1:pn) = elmat(s1:sf,1:pn,s1:sf,1:pn) &
                                   - tau_u*elm_usgs_ES(1:tn,1:pn,1:tn,1:pn)
      ! Assembly elmep to elmat
      elmat(s1:sf,1:pn,p1,1:pn)    = elmat(s1:sf,1:pn,p1,1:pn) &
                                   - tau_u*elm_usgs_EP(1:tn,1:pn,1,1:pn)
      ! Assembly elmqs to elmat
      elmat(p1,1:pn,s1:sf,1:pn)    = elmat(p1,1:pn,s1:sf,1:pn)  &
                                   - tau_u*elm_usgs_QS(1,1:pn,1:tn,1:pn)
      ! Assembly elmqp to elmat
      elmat(p1,1:pn,p1,1:pn)       = elmat(p1,1:pn,p1,1:pn) &
                                   - tau_u*elm_usgs_QP(1,1:pn,1,1:pn)

      !-------------TS------------------
      ! Momentum equation
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn) &
                                   - tau_s*elm_ssgs_VU(1:nd,1:pn,1:nd,1:pn)
      !Assembly elmvs to elmat
      elmat(u1:uf,1:pn,s1:sf,1:pn) = elmat(u1:uf,1:pn,s1:sf,1:pn)  &
                                   - tau_s*elm_ssgs_VS(1:nd,1:pn,1:tn,1:pn)
      !Constitutive equation
      ! Assembly elmeu to elmat
      elmat(s1:sf,1:pn,u1:uf,1:pn) = elmat(s1:sf,1:pn,u1:uf,1:pn) &
                                   - tau_s*elm_ssgs_EU(1:tn,1:pn,1:nd,1:pn)
      ! Assembly elmes to elmat
      elmat(s1:sf,1:pn,s1:sf,1:pn) = elmat(s1:sf,1:pn,s1:sf,1:pn) &
                                   - tau_s*elm_ssgs_ES(1:tn,1:pn,1:tn,1:pn)
   
      !-------------TP------------------
      ! Momentum equation
      ! Assembly elmvu to elmat
      elmat(u1:uf,1:pn,u1:uf,1:pn) = elmat(u1:uf,1:pn,u1:uf,1:pn) &
                                   - tau_p*elm_psgs_VU(1:nd,1:pn,1:nd,1:pn)
      ! Assembly elmvp to elmat
      elmat(u1:uf,1:pn,p1,1:pn)    = elmat(u1:uf,1:pn,p1,1:pn) &
                                   - tau_p*elm_psgs_VP(1:nd,1:pn,1,1:pn)
      ! Continuity Equation
      ! Assembly elmqu to elmat
      elmat(p1,1:pn,u1:uf,1:pn)    = elmat(p1,1:pn,u1:uf,1:pn) &
                                   - tau_p*elm_psgs_QU(1,1:pn,1:nd,1:pn)
      ! Assembly elmqp to elmat
      elmat(p1,1:pn,p1,1:pn)       = elmat(p1,1:pn,p1,1:pn) &
                                   - tau_p*elm_psgs_QP(1,1:pn,1,1:pn)

      !---------------------------------END SGS----------------------------
   end subroutine AssembleStiffnessMatrixSUP_SGS

   module subroutine ResetSUPBaseArrays
      implicit none

      !-----------------Assembly-----------------
      elrhs       = 0.0_rp
      elmat       = 0.0_rp

      !-------------------RHS---------------
      elrhu       = 0.0_rp
      elrhd       = 0.0_rp
      elrhp       = 0.0_rp

      !-------------------u---------------
      elmst       = 0.0_rp
      elmsv       = 0.0_rp
      !-------------------s---------------
      elmuq       = 0.0_rp
      elmut       = 0.0_rp
      !-------------------p---------------
      elmpv       = 0.0_rp
      elmpq       = 0.0_rp

   end subroutine ResetSUPBaseArrays

   module subroutine ResetSUPBaseArrays_SGS_dynsgs
      implicit none
      !-------------------Stabilization---------------
      !-------------------Tau_u---------------
      elm_usgs_VU = 0.0_rp
      elm_usgs_VS = 0.0_rp
      elm_usgs_VP = 0.0_rp

   end subroutine ResetSUPBaseArrays_SGS_dynsgs

   module subroutine ResetSUPBaseArrays_SGS_dyn
      implicit none
      !-------------------Stabilization---------------
      !-------------------Tau_u---------------
      elm_usgs_EU = 0.0_rp
      elm_usgs_QU = 0.0_rp

   end subroutine ResetSUPBaseArrays_SGS_dyn

   module subroutine ResetSUPBaseArrays_SGS
      implicit none
      !-------------------Stabilization---------------
      !-------------------Tau_u---------------
      elm_usgs_ES = 0.0_rp
      elm_usgs_EP = 0.0_rp
      elm_usgs_QS = 0.0_rp
      elm_usgs_QP = 0.0_rp
      !-------------------Tau_s---------------
      elm_ssgs_VU = 0.0_rp
      elm_ssgs_VS = 0.0_rp
      elm_ssgs_EU = 0.0_rp
      elm_ssgs_ES = 0.0_rp
      !-------------------Tau_p---------------
      elm_psgs_VU = 0.0_rp
      elm_psgs_VP = 0.0_rp
      elm_psgs_QU = 0.0_rp
      elm_psgs_QP = 0.0_rp

   end subroutine ResetSUPBaseArrays_SGS

   module subroutine ResetSUPNonlinearDerTerms
      implicit none
      !-------------------Stabilization---------------
      !-------------------Tau_u---------------
      derJFki = 0.0_rp
      grdFmat = 0.0_rp
      grdFspt = 0.0_rp
      divFmat = 0.0_rp
      divFspt = 0.0_rp

   end subroutine ResetSUPNonlinearDerTerms

   module subroutine AllocateBaseElmopeMatricesSUP_SGS_dynsgs
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Displ stabilization
      call sup%Memor%alloc(nd,mn,nd,mn,elm_usgs_VU,'elm_usgs_VU','sldsup_elmope')
      !-----------Stress stabilization
      call sup%Memor%alloc(nd,mn,tn,mn,elm_usgs_VS,'elm_usgs_VS','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%alloc(nd,mn,1 ,mn,elm_usgs_VP,'elm_usgs_VP','sldsup_elmope')

   end subroutine

   module subroutine DeallocateBaseElmopeMatricesSUP_SGS_dynsgs
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Displ stabilization
      call sup%Memor%dealloc(nd,mn,nd,mn,elm_usgs_VU,'elm_usgs_VU','sldsup_elmope')
      !-----------Stress stabilization
      call sup%Memor%dealloc(nd,mn,tn,mn,elm_usgs_VS,'elm_usgs_VS','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%dealloc(nd,mn,1 ,mn,elm_usgs_VP,'elm_usgs_VP','sldsup_elmope')

   end subroutine

   module subroutine AllocateBaseElmopeMatricesSUP_SGS_dyn
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call sup%Memor%alloc(tn,mn,nd,mn,elm_usgs_EU,'elm_usgs_EU','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%alloc(1 ,mn,nd,mn,elm_usgs_QU,'elm_usgs_QU','sldsup_elmope')

   end subroutine

   module subroutine DeallocateBaseElmopeMatricesSUP_SGS_dyn
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call sup%Memor%dealloc(tn,mn,nd,mn,elm_usgs_EU,'elm_usgs_EU','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%dealloc(1 ,mn,nd,mn,elm_usgs_QU,'elm_usgs_QU','sldsup_elmope')

   end subroutine

   module subroutine AllocateBaseElmopeMatricesSUP_SGS
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call sup%Memor%alloc(tn,mn,tn,mn,elm_usgs_ES,'elm_usgs_ES','sldsup_elmope')
      call sup%Memor%alloc(tn,mn,1 ,mn,elm_usgs_EP,'elm_usgs_EP','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%alloc(1 ,mn,tn,mn,elm_usgs_QS,'elm_usgs_QS','sldsup_elmope')
      call sup%Memor%alloc(1 ,mn,1 ,mn,elm_usgs_QP,'elm_usgs_QP','sldsup_elmope')

      !-------------------Tau_s matrix calculation----------
      !-----------Displacement stabilization
      call sup%Memor%alloc(nd,mn,nd,mn,elm_ssgs_VU,'elm_ssgs_VU','sldsup_elmope')
      call sup%Memor%alloc(nd,mn,tn,mn,elm_ssgs_VS,'elm_ssgs_VS','sldsup_elmope')
      !-----------Stress stabilization
      call sup%Memor%alloc(tn,mn,tn,mn,elm_ssgs_ES,'elm_ssgs_ES','sldsup_elmope')
      call sup%Memor%alloc(tn,mn,nd,mn,elm_ssgs_EU,'elm_ssgs_EU','sldsup_elmope')

      !-------------------Tau_p matrix calculation----------
      !-----------Displacement stabilization
      call sup%Memor%alloc(nd,mn,nd,mn,elm_psgs_VU,'elm_psgs_VU','sldsup_elmope')
      call sup%Memor%alloc(nd,mn,1 ,mn,elm_psgs_VP,'elm_psgs_VP','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%alloc(1 ,mn,nd,mn,elm_psgs_QU,'elm_psgs_QU','sldsup_elmope')
      call sup%Memor%alloc(1 ,mn,1 ,mn,elm_psgs_QP,'elm_psgs_QP','sldsup_elmope')

   end subroutine

   module subroutine DeallocateBaseElmopeMatricesSUP_SGS
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      !-------------------Tau_u matrix calculation----------
      !-----------Stress stabilization
      call sup%Memor%dealloc(tn,mn,tn,mn,elm_usgs_ES,'elm_usgs_ES','sldsup_elmope')
      call sup%Memor%dealloc(tn,mn,1 ,mn,elm_usgs_EP,'elm_usgs_EP','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%dealloc(1 ,mn,tn,mn,elm_usgs_QS,'elm_usgs_QS','sldsup_elmope')
      call sup%Memor%dealloc(1 ,mn,1 ,mn,elm_usgs_QP,'elm_usgs_QP','sldsup_elmope')

      !-------------------Tau_s matrix calculation----------
      !-----------Displacement stabilization
      call sup%Memor%dealloc(nd,mn,nd,mn,elm_ssgs_VU,'elm_ssgs_VU','sldsup_elmope')
      call sup%Memor%dealloc(nd,mn,tn,mn,elm_ssgs_VS,'elm_ssgs_VS','sldsup_elmope')
      !-----------Stress stabilization
      call sup%Memor%dealloc(tn,mn,tn,mn,elm_ssgs_ES,'elm_ssgs_ES','sldsup_elmope')
      call sup%Memor%dealloc(tn,mn,nd,mn,elm_ssgs_EU,'elm_ssgs_EU','sldsup_elmope')

      !-------------------Tau_p matrix calculation----------
      !-----------Displacement stabilization
      call sup%Memor%dealloc(nd,mn,nd,mn,elm_psgs_VU,'elm_psgs_VU','sldsup_elmope')
      call sup%Memor%dealloc(nd,mn,1 ,mn,elm_psgs_VP,'elm_psgs_VP','sldsup_elmope')
      !-----------Pressure stabilization
      call sup%Memor%dealloc(1 ,mn,nd,mn,elm_psgs_QU,'elm_psgs_QU','sldsup_elmope')
      call sup%Memor%dealloc(1 ,mn,1 ,mn,elm_psgs_QP,'elm_psgs_QP','sldsup_elmope')

   end subroutine

   module subroutine AllocateSUPAssemblyMatrices
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      mn  = e%mnode

      call sup%Memor%alloc(tn,mn,      elrhd,'elrhd','sldsup_Elmope')
      call sup%Memor%alloc(1,mn,       elrhp,'elrhp','sldsup_Elmope')
      !-----------Displacement testing
      call sup%Memor%alloc(nd,mn,tn,mn,elmsv,'elmsv','sldsup_elmope')
      call sup%Memor%alloc(nd,mn,1 ,mn,elmpv,'elmpv','sldsup_elmope')
      !-----------Stress testing
      call sup%Memor%alloc(tn,mn,tn,mn,elmst,'elmst','sldsup_elmope')
      call sup%Memor%alloc(tn,mn,nd,mn,elmut,'elmut','sldsup_elmope')
      !-----------Pressure testing
      call sup%Memor%alloc(1 ,mn,nd,mn,elmuq,'elmuq','sldsup_elmope')
      call sup%Memor%alloc(1 ,mn,1 ,mn,elmpq,'elmpq','sldsup_elmope')

   end subroutine

   module subroutine DeallocateSUPAssemblyMatrices
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      mn  = e%mnode

      call sup%Memor%dealloc(tn,mn,      elrhd,'elrhd','sldsup_Elmope')
      call sup%Memor%dealloc(1,mn,       elrhp,'elrhp','sldsup_Elmope')
      !-----------Displacement testing
      call sup%Memor%dealloc(nd,mn,tn,mn,elmsv,'elmsv','sldsup_elmope')
      call sup%Memor%dealloc(nd,mn,1 ,mn,elmpv,'elmpv','sldsup_elmope')
      !-----------Stress testing
      call sup%Memor%dealloc(tn,mn,tn,mn,elmst,'elmst','sldsup_elmope')
      call sup%Memor%dealloc(tn,mn,nd,mn,elmut,'elmut','sldsup_elmope')
      !-----------Pressure testing
      call sup%Memor%dealloc(1 ,mn,nd,mn,elmuq,'elmuq','sldsup_elmope')
      call sup%Memor%dealloc(1 ,mn,1 ,mn,elmpq,'elmpq','sldsup_elmope')

   end subroutine

   module subroutine AllocateBaseElmopeMatricesSUPLinear
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn  =(nd*(nd+1))/2
      mn  = e%mnode

      call sup%Memor%alloc(tn,tn,      P,      'P','sldsup_Elmope')
      call sup%Memor%alloc(tn,tn,      D,      'D','sldsup_Elmope')
      call sup%Memor%alloc(tn,tn,  C_dev,  'C_dev','sldsup_Elmope')
      call sup%Memor%alloc(tn,tn,  D_dev,  'D_dev','sldsup_Elmope')
   end subroutine

   module subroutine DeallocateBaseElmopeMatricesSUPLinear
      implicit none
      integer(ip) :: mn,nd,tn

      call a%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      mn  = e%mnode

      call sup%Memor%dealloc(tn,tn,      P,      'P','sldsup_Elmope')
      call sup%Memor%dealloc(tn,tn,      D,      'D','sldsup_Elmope')
      call sup%Memor%dealloc(tn,tn,  C_dev,  'C_dev','sldsup_Elmope')
      call sup%Memor%dealloc(tn,tn,  D_dev,  'D_dev','sldsup_Elmope')
   end subroutine

end submodule sldsup_BaseElmope
