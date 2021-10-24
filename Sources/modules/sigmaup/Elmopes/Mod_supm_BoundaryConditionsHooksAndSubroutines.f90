module Mod_supm_BoundaryConditionsHooksAndSubroutines
   use typre
   use Mod_supm_BaseElmope
   use Mod_nsm_Viscosity
   implicit none
   !private
   public SetPointersBoundaryConditions
   integer(ip)                                     :: bcstart      
   integer(ip),   allocatable                      :: kfl_IsSet
   real(rp)                                        :: alfa_inv,vista,updbcn,delta
   real(rp),   allocatable,   dimension(:)         :: tract,tractexact
   real(rp),   allocatable,   dimension(:,:)       :: etrac
   !Temperature coupling
   real(rp), allocatable, dimension(:,:)           :: exsigr
   real(rp), allocatable, dimension(:)             :: exsig, eltem
   !Isentropic
   real(rp),   allocatable,   dimension(:,:)       :: ElemVelocMean
   real(rp),   allocatable,   dimension(:,:)       :: BoundVelocMean
   real(rp),   allocatable,   dimension(:,:)       :: BoundPressMean
   
   real(rp),   allocatable,   dimension(:)         :: GpVelocMean
   real(rp),   allocatable,   dimension(:,:)       :: GpPressMean
   real(rp),   allocatable,   dimension(:,:)       :: GradVelocMean
   real(rp)                                        :: DivVelocMean
   real(rp), external :: funcre
   real(rp),   allocatable,   dimension(:)         :: vel_wall_vec
   real(rp),   allocatable,   dimension(:,:,:,:)   :: BoundElmatVU
   real(rp),   allocatable,   dimension(:,:)       :: BoundElrhsU
   real(rp),   allocatable,   dimension(:,:,:,:)   :: BoundElmatVP
   real(rp),   allocatable,   dimension(:,:,:,:)   :: BoundElmatQU
   real(rp),   allocatable,   dimension(:,:)       :: BoundElrhsP
   real(rp),   allocatable,   dimension(:,:,:,:)   :: wmatr
   real(rp),   allocatable,   dimension(:,:)       :: wrhsi
   real(rp),   allocatable,   dimension(:)         :: gpcod
contains
   
   !----------------------------------------------------------------------------------------------!
   !                                 SET POINTERS                                                 !
   !----------------------------------------------------------------------------------------------!
   
   subroutine SetPointersBoundaryConditions(itask)
      implicit none
      integer(ip)    :: itask
      character(6)   :: task
      
      logical        :: aux_logic
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            !Allocates
            call ConcatenateProcedures(ProcHook_PreGauss,AllocateVariablesForBoundaryOperations)
            !ElmatsToZero
            call ConcatenateProcedures(ProcHook_ElmatsToZero,BoundElmatsToZero)
            !FSI
            if(associated(a%etraction) .and. a%doRobin .eqv. .true.) then
               call ConcatenateProcedures(ProcHook_Gathers,GatherFSI)
            endif
            !Smagorinsky
!                if (a%kfl_cotur< 0.or.a%MatProp(imat)%lawvi/=0.or.a%kfl_ExchangeLocationWallLaw==1) &
!                call ConcatenateProcedures(ProcHook_Gathers,GatherSmagorinsky)
!                if (a%kfl_cotur<0) call ConcatenateProcedures(ProcHook_InGauss,InGaussSmagorinsky)
            !Non-newtonian viscosity
            if(a%MatProp(imat)%lawvi /= 0) call ConcatenateProcedures(ProcHook_InGauss,InGaussNonNewtonianViscosity)
            !WLF temperature model 
            if (a%kfl_cotem_WLF==1) call ConcatenateProcedures(ProcHook_Gathers,GathersWLFmodel)
            if (a%kfl_cotem_Arrhenius==1) call ConcatenateProcedures(ProcHook_Gathers,GathersArrheniusmodel)
            
            !Assembly
            call ConcatenateProcedures(ProcHook_PostLoop,AssembleBoundElmats)
            !Deallocates
            call ConcatenateProcedures(ProcHook_Finalizations,DeallocateVariablesForBoundaryOperations)
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine 

   !--------------------------------------------------------------------------------------------!
   !                                      HOOK SUBROUTINES                                      !
   !--------------------------------------------------------------------------------------------!
   subroutine BoundElmatsToZero
         wmatr=0.0_rp
         wrhsi=0.0_rp
         tract = 0.0_rp
         BoundElmatVU = 0.0_rp
         BoundElmatVP = 0.0_rp
         BoundElmatQU = 0.0_rp
         BoundElrhsU  = 0.0_rp
         BoundElrhsP  = 0.0_rp
   end subroutine
   
   subroutine GatherSmagorinsky
      implicit none
      call e%gather(e%ndime,elvel,a%veloc(:,:,1))
   end subroutine

   subroutine GatherFSI
      implicit none
      alfa_inv=0.0_rp
      if(a%alfa_robin < 0.0_rp) then
         if(a%MPIrank == a%MPIroot) write(*,*) &
         'Warning alpha robin = 0, imposing standard Dirichlet B.C on Fluid'
              
         !Reset flag so it does not add traction terms
         a%doRobin = .false.
      else
         alfa_inv=1.0_rp/a%alfa_robin
         etrac=0.0_rp
         call e%gatherb(e%ndime,etrac,a%etraction(:,:))
      endif
   end subroutine
   
   subroutine GathersWLFmodel
      implicit none
      call e%gather(1,eltem,a%tempe)
      call sup_templaw_WLF(e, eltem,  a%ReferenceTemp, a%c1_WLF, a%c2_WLF, a%MatProp(imat)%LawViParam, acvis, lambda)
   end subroutine   
   
   subroutine GathersArrheniusmodel
      implicit none
      call e%gather(1,eltem,a%tempe)
      call sup_templaw_Arrhenius(e, eltem, a%ReferenceTemp, a%alpha_Arrhenius, a%MatProp(imat)%LawViParam, acvis, lambda)
   end subroutine   
   
   subroutine InGaussSmagorinsky 
      implicit none
      call e%elmlen
      call e%gradientb(e%ndime,elvel,grvel)
      if (a%kfl_cotur == -1) then   
         call nsm_smago(e,grvel,acden,a%turbu(1),vista)
      elseif (a%kfl_cotur == -2) then
         call nsm_wale(e,grvel,acden,a%turbu(1),vista)
      endif
      acvis = acvis + vista
   end subroutine

   subroutine InGaussNonNewtonianViscosity
      implicit none

      call e%gradientb(e%ndime,elvel,grvel)
      call nsi_vislaw(e,grvel,a%MatProp(imat)%lawvi,a%MatProp(imat)%LawViParam,acvis)
   end subroutine

   subroutine AssembleBoundElmats
      implicit none
      integer(ip)    :: idime,inode,inodb
      !Assemble wmatr
      elmat((bcstart+1):(bcstart+e%ndime+1),:,(bcstart+1):(bcstart+e%ndime+1),:)= dsurf*wmatr + &
            elmat((bcstart+1):(bcstart+e%ndime+1),:,(bcstart+1):(bcstart+e%ndime+1),:)            
      
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         do idime=1,e%ndime
            elrhs(bcstart+idime,inode)=elrhs(bcstart+idime,inode)&
            +dsurf*(wrhsi(idime,inode)+tract(idime)*e%shapb(inodb,e%igaub)&
            -alfa_inv*etrac(idime,inodb)*e%shapb(inodb,e%igaub))
         end do    
      end do

      !Assemble BoundElmat and BoundElrhs
      elmat(1:e%ndime,:,1:e%ndime,:) = elmat(1:e%ndime,:,1:e%ndime,:)+BoundElmatVU
      elmat(1:e%ndime,:,e%ndime+1,:) = elmat(1:e%ndime,:,1+e%ndime,:)+BoundElmatVP(1:e%ndime,:,1,:)
      elmat(e%ndime+1,:,1:e%ndime,:) = elmat(e%ndime+1,:,1:e%ndime,:)+BoundElmatQU(1,:,1:e%ndime,:)

      elrhs(1:e%ndime,:) = elrhs(1:e%ndime,:) + BoundElrhsU
      elrhs(e%ndime+1,:) = elrhs(e%ndime+1,:) + BoundElrhsP(1,:)
   end subroutine


   subroutine AllocateVariablesForBoundaryOperations
      implicit none
      call a%Memor%alloc(e%ndime,tract,'tract','nsi_bouope')   
      call a%Memor%alloc(e%ndime,e%pnodb,etrac,'etrac','nsi_bouope')   
      call a%Memor%alloc(e%ndime,e%mnode,ElemVelocMean,'ElemVelocMean','bouope')
      call a%Memor%alloc(e%ndime,e%mnodb,BoundVelocMean,'ElemVelocMean','bouope')
      call a%Memor%alloc(e%mnodb,1,BoundPressMean,'ElemVelocMean','bouope')
      call a%Memor%alloc(e%ndime,GpVelocMean,'GpVelocMean','bouope')
      call a%Memor%alloc(1,1,GpPressMean,'GpVelocMean','bouope')
      call a%Memor%alloc(e%ndime,e%ndime,GradVelocMean,'GpVelocMean','bouope')
      call a%Memor%alloc(e%ndime,vel_wall_vec,'vel_wall_vec','bouope')
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,BoundElmatVU,'BoundElmatVU','bouope')
      call a%Memor%alloc(e%ndime,e%mnode,BoundElrhsU,'BoundElrhsU','bouope')
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,BoundElmatQU,'BoundElmaQU','bouope')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,BoundElmatVP,'BoundElmaVP','bouope')
      call a%Memor%alloc(1,e%mnode,BoundElrhsP,'BoundElrhsP','bouope')
      call a%Memor%alloc(e%ndime+1,e%mnode,e%ndime+1,e%mnode,wmatr,'wmatr','bouope')
      call a%Memor%alloc(e%ndime+1,e%mnode,wrhsi,'wrhsi','bouope')
      call a%Memor%alloc(e%ndime,gpcod,'gpcod','bouope')
      call a%Memor%alloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supm_EnditeElmope')
      
      !SigmaUP
      call a%Memor%alloc(bcstart,tractexact,'tractexact','supm_bouop0')
      !Temperature
      call a%Memor%alloc(e%mnode,eltem,'eltem','supm_bouop0')
      !ExactSol 
      call a%Memor%alloc(e%mnode,exsig,'exsig','supm_bouop0')
      call a%Memor%alloc(bcstart,e%ndime,exsigr,'exsigr','supm_bouop0')
      
   endsubroutine

   subroutine DeallocateVariablesForBoundaryOperations
      implicit none
      call a%Memor%dealloc(e%ndime,tract,'tract','nsi_bouope')   
      call a%Memor%dealloc(e%ndime,e%pnodb,etrac,'etrac','nsi_bouope')   
      call a%Memor%dealloc(e%ndime,e%mnode,ElemVelocMean,'ElemVelocMean','bouope')
      call a%Memor%dealloc(e%ndime,e%mnodb,BoundVelocMean,'ElemVelocMean','bouope')
      call a%Memor%dealloc(e%mnodb,1,BoundPressMean,'ElemVelocMean','bouope')
      call a%Memor%dealloc(e%ndime,GpVelocMean,'GpVelocMean','bouope')
      call a%Memor%dealloc(1,1,GpPressMean,'GpVelocMean','bouope')
      call a%Memor%dealloc(e%ndime,e%ndime,GradVelocMean,'GpVelocMean','bouope')
      call a%Memor%dealloc(e%ndime,vel_wall_vec,'vel_wall_vec','bouope')
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,BoundElmatVU,'BoundElmatVU','bouope')
      call a%Memor%dealloc(e%ndime,e%mnode,BoundElrhsU,'BoundElrhsU','bouope')
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,BoundElmatQU,'BoundElmaQU','bouope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,BoundElmatVP,'BoundElmaVP','bouope')
      call a%Memor%dealloc(1,e%mnode,BoundElrhsP,'BoundElrhsP','bouope')
      call a%Memor%dealloc(e%ndime+1,e%mnode,e%ndime+1,e%mnode,wmatr,'wmatr','bouope')
      call a%Memor%dealloc(e%ndime+1,e%mnode,wrhsi,'wrhsi','bouope')
      call a%Memor%dealloc(e%ndime,gpcod,'gpcod','bouope')
      call a%Memor%dealloc(e%ndime,a%ncomp-1,gpvel,'gpvel','supm_EnditeElmope')
      
      !SigmaUP
      call a%Memor%dealloc(bcstart,tractexact,'tractexact','nsm_bouop0')
      !Temperature
      call a%Memor%dealloc(e%mnode,eltem,'eltem','supm_bouop0')
      !ExactSol 
      call a%Memor%dealloc(e%mnode,exsig,'exsig','supm_bouop0')
      call a%Memor%dealloc(bcstart,e%ndime,exsigr,'exsigr','supm_bouop0')
   endsubroutine

   subroutine nsm_bouope_pressure
      implicit none
      if(a%kfl_funbo(iboun)/=0) then
         updbcn=funcre(a%funpa(a%kfl_funbo(iboun))%a,&
         a%kfl_funty(a%kfl_funbo(iboun),2),&
         a%kfl_funty(a%kfl_funbo(iboun),1),a%bctime)
      else
         updbcn=1.0_rp
      end if
      tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*updbcn*e%baloc(1:e%ndime,e%ndime)
   end subroutine

   subroutine nsm_bouope_WallLaw
      implicit none
      integer(ip)       :: inode
      real(rp)          :: ivel_wall,vel_wall

      !nwall=nwall+1
      if (a%kfl_ExchangeLocationWallLaw == 0) then
         if (a%kfl_FixWallLawDeltaForAll == 1) then
            delta = a%WallLawDeltaValue
         else
            delta = a%bvnat(iboun)%a(1)
         endif
         vel_wall_vec = gpvel(:,1)
      elseif (a%kfl_ExchangeLocationWallLaw == 1) then
         !We set as delta the maximum element length
            call e%elmdcg
            call e%elmlen
            delta = e%hleng(1)
         !We set as velocity the maximum velocity of the element, which will correspond to interior nodes
         vel_wall = 0
         vel_wall_vec = 0.0_rp
         do inode = 1,e%pnode
            ivel_wall = sqrt(dot_product(elvel(:,inode,1),elvel(:,inode,1)))
            if (ivel_wall > vel_wall) then
               vel_wall = ivel_wall
               vel_wall_vec = elvel(:,inode,1)
            endif
         enddo   
      endif   
   end subroutine
end module
