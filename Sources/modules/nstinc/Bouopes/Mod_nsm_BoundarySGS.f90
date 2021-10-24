module Mod_nsm_BoundarySGS
   use Mod_PointerSetter
   use Mod_nsm_BaseBouope
   implicit none
   private
   public SetPointersBoundarySGS

   type, extends(PointerSetter) :: SPBoundarySGS
contains
      procedure :: SpecificSet => SpecificSetBoundarySGS
   end type
   type(SPBoundarySGS) :: SetPointersBoundarySGS

   real(rp), allocatable, dimension(:)   :: gpetrac,gpevel
   real(rp), allocatable, dimension(:,:) :: eletrac,elevel
   real(rp)                              :: btau,gpetau(1),teS,teN,bbeta,ebeta
   real(rp), allocatable, dimension(:)   :: elvmass,elbtau,eletau
   
contains

   subroutine SpecificSetBoundarySGS(d)
      implicit none
      class(SPBoundarySGS) :: d
  
      !Boundary SGS
      if (a%kfl_bousg /= 0) then
         call ConcatenateProcedures(ProcHook_Initializations,AllocateArraysSGS2)
         call ConcatenateProcedures(ProcHook_Initializations,ToZero)

         call ConcatenateProcedures(ProcHook_ElmatsToZero,ArraysToZeroSGS2)
         call ConcatenateProcedures(ProcHook_Gathers,GatherSGS2)
         call ConcatenateProcedures(ProcHook_Interpolates,InterpolateSGS2)

         call ConcatenateProcedures(ProcHook_InGaussElmats,ElmatsBoundary)
         call ConcatenateProcedures(ProcHook_InGaussElmats,ElrhsBoundary)
         call ConcatenateProcedures(ProcHook_PreDirichlet,AssemblyBaseElmats)
         
         call ConcatenateProcedures(ProcHook_PreGauss,ComputeTauSGS)
         call ConcatenateProcedures(ProcHook_InGauss,ComputeTauMeans)
         call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,BoundaryMassMatrix)
         call ConcatenateProcedures(ProcHook_Assembly,AssemblyTau)
         call ConcatenateProcedures(ProcHook_Finalizations,ExtraTau)
         
         call ConcatenateProcedures(ProcHook_InGauss,ComputeBoundarySubgridScales)

         call ConcatenateProcedures(ProcHook_Finalizations,DeallocateArraysSGS2)
      end if
   end subroutine

   !-----------------------------------------------------------------------
   ! Boundary Tau

   subroutine ToZero
      implicit none

      a%btau = 0.0_rp
      vmassBoundary = 0.0_rp
   end subroutine

   subroutine ComputeTauSGS
      implicit none
      integer(ip) :: inodb,inode

      if (a%kfl_fixbo(iboun) == 51 .or. a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         call ComputeTauBoundary(e,chale,acvis,btau)
         bbeta = 0.0_rp
         if (btau > 0) bbeta = 1/btau
         do inodb = 1,e%pnodb
            inode= e%lboel(inodb)  !inod to ipoin
            elbtau(inode) = btau
         end do
      end if
   end subroutine

   subroutine ComputeTauMeans
      implicit none
      integer(ip) :: inodb,inode

      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         teS = a%bstco(1)*(btau + gpetau(1))
         ebeta = 0.0_rp
         if (gpetau(1) > 0) ebeta = 1/gpetau(1)
         teN = 2*a%bstco(2)*bbeta*ebeta/(ebeta+bbeta)  ! Full harmonic mean Nitsche
      end if
   end subroutine

   subroutine BoundaryMassMatrix
      implicit none
      integer(ip) :: inodb,inode

      do inodb = 1,e%pnodb
         inode= e%lboel(inodb)
         elvmass(inode) = elvmass(inode) + e%shapb(inodb,e%igaub)*dsurf
      end do
   end subroutine
   
   subroutine AssemblyTau
      implicit none

      if (a%kfl_fixbo(iboun) == 51 .or. a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         call a%Mesh%AssemblyToArray(e,1_ip,elvmass,vmassBoundary)
         call a%Mesh%AssemblyToArray(e,1_ip,elbtau,a%btau)
      end if
   end subroutine

   subroutine ExtraTau
      implicit none
      integer(ip) :: ipoin,npoin

      call a%Mesh%GetNpoin(npoin)     
      do ipoin = 1,npoin
         if (vmassBoundary(ipoin) > 0.0_rp) a%btau(ipoin) = a%btau(ipoin)/vmassBoundary(ipoin)
      enddo
   end subroutine

   !-----------------------------------------------------------------------------
   ! Bouope for Neumann condition (iboun = 52)
   ! (F_2,F_2) lhs
   ! (F_1,F_2) + esgs rhs

   subroutine AllocateArraysSGS2
      implicit none
      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)     

      call a%Memor%alloc(e%ndime,e%mnodb,eletrac,'eletrac','nsi_bousgs')
      call a%Memor%alloc(e%ndime,gpetrac,'gpetrac','nsi_bousgs')
      call a%Memor%alloc(e%ndime,e%mnodb,elevel,'elevel','nsi_bousgs')
      call a%Memor%alloc(e%ndime,gpevel,'gpevel','nsi_bousgs')

      call a%Memor%alloc(e%mnode,elbtau,'elbtau','nsi_bousgs')
      call a%Memor%alloc(e%mnode,eletau,'eletau','nsi_bousgs')

      call a%Memor%alloc(e%mnode,elvmass,'elvmass','nsi_bousgs')
      call a%Memor%alloc(npoin,vmassBoundary,'vmassBoundary','nsi_bousgs')

      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','nsi_bousgs')
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','nsi_bousgs')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','nsi_bousgs')
      call a%Memor%alloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','nsi_bousgs')
      call a%Memor%alloc(e%ndime,e%mnode,elrhu,'elrhu','nsi_bousgs')
      call a%Memor%alloc(1,e%mnode,elrhp,'elrhp','nsi_bousgs')
   end subroutine

   subroutine DeallocateArraysSGS2
      implicit none
      integer(ip) :: npoin

      call a%Mesh%GetNpoin(npoin)     

      call a%Memor%dealloc(e%ndime,e%mnodb,eletrac,'eletrac','nsi_bousgs')
      call a%Memor%dealloc(e%ndime,gpetrac,'gpetrac','nsi_bousgs')
      call a%Memor%dealloc(e%ndime,e%mnodb,elevel,'elevel','nsi_bousgs')
      call a%Memor%dealloc(e%ndime,gpevel,'gpevel','nsi_bousgs')

      call a%Memor%dealloc(e%mnode,elbtau,'elbtau','nsi_bousgs')
      call a%Memor%dealloc(e%mnode,eletau,'eletau','nsi_bousgs')

      call a%Memor%dealloc(e%mnode,elvmass,'elvmass','nsi_bousgs')
      call a%Memor%dealloc(npoin,vmassBoundary,'vmassBoundary','nsi_bousgs')

      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuv,'elmuv','nsi_bousgs')
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuq,'elmuq','nsi_bousgs')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmpv,'elmpv','nsi_bousgs')
      call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmpq,'elmpq','nsi_bousgs')
      call a%Memor%dealloc(e%ndime,e%mnode,elrhu,'elrhu','nsi_bousgs')
      call a%Memor%dealloc(1,e%mnode,elrhp,'elrhp','nsi_bousgs')
   end subroutine

   subroutine ArraysToZeroSGS2
      implicit none

      elmuv = 0.0_rp
      elmpv = 0.0_rp
      elmuq = 0.0_rp
      elmpq = 0.0_rp
      elrhu = 0.0_rp
      elrhp = 0.0_rp
   end subroutine

   subroutine GatherSGS2
      implicit none
      real(rp) :: aux_tau(1)
      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         call e%gatherb(e%ndime,eletrac,a%etraction)
         call e%gatherb(e%ndime,elevel,a%eveloc)
         call e%gatherb(1,eletau,a%etau)
      end if
   end subroutine

   subroutine InterpolateSGS2
      implicit none
      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         call e%interpb(e%ndime,eletrac,gpetrac)
         call e%interpb(e%ndime,elevel,gpevel)
         call e%interpb(1,eletau,gpetau)
      end if
   end subroutine

   subroutine ElmatsBoundary
      implicit none
      !Compute contributions to LHS:
      if (a%kfl_fixbo(iboun) == 52) then
         !Stabilization
         call nsm_elmbuv_tr(e,a%fvins,dsurf,acvis,teS,elmuv)
         call nsm_elmbpv_tr(e,a%fvins,dsurf,acvis,teS,elmpv)
         call nsm_elmbuq_tr(e,a%fvins,dsurf,acvis,teS,elmuq)
         call nsm_elmbpq_tr(e,dsurf,teS,elmpq)
      elseif (a%kfl_fixbo(iboun) == 53) then
         !Stabilization
         call nsm_elmbuv_tr(e,a%fvins,dsurf,acvis,teS,elmuv)
         call nsm_elmbpv_tr(e,a%fvins,dsurf,acvis,teS,elmpv)
         call nsm_elmbuq_tr(e,a%fvins,dsurf,acvis,teS,elmuq)
         call nsm_elmbpq_tr(e,dsurf,teS,elmpq)
         !Nitsche
         call nsm_elmbuv_ni(e,dsurf,teN,elmuv)
         !Consistency
         call nsm_elmbuv_cn(e,a%fvins,dsurf,acvis,elmuv)
         call nsm_elmbpv_cn(e,dsurf,elmpv)
         call nsm_elmbuv_acn(e,a%fvins,dsurf,acvis,elmuv)
         call nsm_elmbuq_acn(e,dsurf,elmuq)
      end if
   end subroutine

   subroutine ElrhsBoundary
      implicit none
      !Compute contributions to RHS:
      if (a%kfl_fixbo(iboun) == 52) then
         !Stabilization
         call nsm_elmrhu_tr(e,a%fvins,dsurf,acvis,teS,gpetrac,elrhu)
         call nsm_elmrhp_tr(e,dsurf,teS,gpetrac,elrhp)
      elseif (a%kfl_fixbo(iboun) == 53) then
         !Stabilization
         call nsm_elmrhu_tr(e,a%fvins,dsurf,acvis,teS,gpetrac,elrhu)
         call nsm_elmrhp_tr(e,dsurf,teS,gpetrac,elrhp)
         !Nitsche
         call nsm_elmrhu_ni(e,dsurf,teN,gpevel,elrhu)
         !Consistency
         call nsm_elmrhu_cn(e,dsurf,gpetrac,elrhu)
         call nsm_elmrhu_acn(e,a%fvins,dsurf,acvis,gpevel,elrhu)
         call nsm_elmrhp_acn(e,dsurf,gpevel,elrhp)
      end if
   end subroutine

   subroutine AssemblyBaseElmats
      implicit none
      integer(ip) :: idime

      ! Assemble rest of elmat and elrhs
      if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         elmat(1:e%ndime,:,1:e%ndime,:) = elmat(1:e%ndime,:,1:e%ndime,:) + elmuv(1:e%ndime,:,1:e%ndime,:)
         elmat(1:e%ndime,:,e%ndime+1,:) = elmat(1:e%ndime,:,1+e%ndime,:) + elmpv(1:e%ndime,:,1,:)
         elmat(e%ndime+1,:,1:e%ndime,:) = elmat(e%ndime+1,:,1:e%ndime,:) + elmuq(1,:,1:e%ndime,:)
         elmat(e%ndime+1,:,e%ndime+1,:) = elmat(e%ndime+1,:,e%ndime+1,:) + elmpq(1,:,1,:)

         elrhs(1:e%ndime,:) = elrhs(1:e%ndime,:) + elrhu(1:e%ndime,:)
         elrhs(e%ndime+1,:) = elrhs(e%ndime+1,:) + elrhp(1,:)
      end if
   end subroutine

!---------------------------------------------------------------------------------------------------------

   subroutine ComputeBoundarySubgridScales
      implicit none
      ! memory leak convergence problems
      !if (a%kfl_fixbo(iboun) == 52 .or. a%kfl_fixbo(iboun) == 53) then
         !a%bvesgs(iboun)%a(1:e%ndime,e%igaub) = teS*(gpetrac+tract)
      !end if
   end subroutine

end module
