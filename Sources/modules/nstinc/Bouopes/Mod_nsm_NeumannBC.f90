module Mod_nsm_NeumannBC
   use Mod_PointerSetter
   use Mod_nsm_BaseBouope
   implicit none
   private
   public SetPointersNeumannBC, wmatr

   type, extends(PointerSetter) :: SPNeumannBC
contains
      procedure :: SpecificSet => SpecificSetNeumannBC
   end type
   type(SPNeumannBC) :: SetPointersNeumannBC

   real(rp), external :: funcre
   real(rp)           :: updbcn=0.0_rp,delta=0.0_rp
   real(rp)           :: density_wall,viscosity_wall
   real(rp)           :: dynpr,coorg,yplus

   real(rp),    allocatable, dimension(:)       :: vel_wall_vec
   real(rp),    allocatable, dimension(:,:,:,:) :: wmatr

contains

   subroutine SpecificSetNeumannBC(d)
      implicit none
      class(SPNeumannBC) :: d
   
      call ConcatenateProcedures(ProcHook_Initializations,AllocateArrays)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocateArrays)
      call ConcatenateProcedures(ProcHook_InGaussElmats,BoundaryInGaussElmats)
      call ConcatenateProcedures(ProcHook_InGaussElmatsAssembly,AssemblyNeumann)
   end subroutine

   subroutine AllocateArrays
      implicit none
      call a%Memor%alloc(e%ndime,vel_wall_vec,'vel_wall_vec','nsm_bouope')
      call a%Memor%alloc(e%ndime+1,e%mnode,e%ndime+1,e%mnode,wmatr,'wmatr','nsm_bouope')
   end subroutine

   subroutine DeallocateArrays
      implicit none
      call a%Memor%dealloc(e%ndime,vel_wall_vec,'vel_wall_vec','nsm_bouope')
      call a%Memor%dealloc(e%ndime+1,e%mnode,e%ndime+1,e%mnode,wmatr,'wmatr','nsm_bouope')
   end subroutine

   subroutine AssemblyNeumann
      implicit none
      integer(ip) :: idime,inode,inodb

      ! Assemble wmatr
      elmat(1:e%ndime+1,:,1:e%ndime+1,:)= dsurf*wmatr + elmat(1:e%ndime+1,:,1:e%ndime+1,:)            
        
      ! Assemble wrhsi
      do inodb=1,e%pnodb
         inode = e%lboel(inodb)
         elrhs(1:e%ndime,inode)=elrhs(1:e%ndime,inode) + e%shapb(inodb,e%igaub)*dsurf*tract(1:e%ndime)
      end do
   end subroutine

   subroutine BoundaryInGaussElmats
      select case(a%kfl_fixbo(iboun))

      case(2)  !Pressure
         call nsm_bouope_pressure

      case(3)  !Wall law
         !Wall law: sig.n=-rho*(U*^2)*u/|u|
         call a%GetPhysicalParameters(imat,density_wall,viscosity_wall)
         call nsm_bouope_WallLaw
         call nsm_bouwal(e,vel_wall_vec,viscosity_wall,density_wall,delta,wmatr,tract,yplus)
         a%ypmin=min(a%ypmin,yplus)
         a%ypmax=max(a%ypmax,yplus)
         a%ypmean=a%ypmean+yplus
         a%nmean_walllaw = a%nmean_walllaw + 1
         if(abs(a%grnor) /= 0.0_rp) call nsm_bopres(e,wmatr)

      case(5)  !Dynamic pressure
         !Dynamic pressure: sig.n=-(-1/2*rho*u^2) n
         dynpr=-0.5_rp*acden*gpvno*gpvno 
         tract(1:e%ndime)=-dynpr*e%baloc(1:e%ndime,e%ndime)

      case(6)  !Open flow
         !Open boundary: assemble 2*mu*Sym(grad(u).n + n·v p
         call nsm_bouopb(e,wmatr,acvis,a%fvins)

      case(8)  !Atmospheric stress
         !Atmospheric stress:  sig.n = - (p_0 + p_atm) n where p_0 is given
         tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*e%baloc(1:e%ndime,e%ndime)
         coorg = dot_product(gpcod(1:e%ndime),a%gravi(1:e%ndime))
         tract(1:e%ndime) = tract(1:e%ndime) - acden*a%grnor*coorg*e%baloc(1:e%ndime,e%ndime)

      case(9)  !Direct traction prescription
         !Prescribed traction: sig.n = tract_0 (tract_0 given)
         tract(1:e%ndime) = a%bvnat(iboun)%a(1:e%ndime)

      case(10) !Pressure prescribed in a weak form
         !Prescribed Pressure in weak: 
         call nsm_bouopb_p(e,wmatr,acvis,a%fvins)
      end select
   end subroutine

   !********************************************************************************************!
   !                                     BC SUBROUTINES                                         !
   !********************************************************************************************!

   ! Pressure: simga·n = -p*n
   subroutine nsm_bouope_pressure
      if(a%kfl_funbo(iboun)/=0) then
         updbcn=funcre(a%funpa(a%kfl_funbo(iboun))%a,&
               a%kfl_funty(a%kfl_funbo(iboun),2),&
               a%kfl_funty(a%kfl_funbo(iboun),1),a%bctime)
      else
         updbcn=1.0_rp
      end if
      tract(1:e%ndime)=-a%bvnat(iboun)%a(1)*updbcn*e%baloc(1:e%ndime,e%ndime)
   end subroutine

   ! Wall law: sig.n=-rho*(U*^2)*u/|u|
   subroutine nsm_bouope_WallLaw
      implicit none
      integer(ip)       :: inode
      real(rp)          :: ivel_wall,vel_wall

      ! nwall=nwall+1
            
      if (a%kfl_ExchangeLocationWallLaw == 0) then
         if (a%kfl_FixWallLawDeltaForAll == 1) then
            delta = a%WallLawDeltaValue
         else
            delta = a%bvnat(iboun)%a(1)
         endif
            vel_wall_vec = gpbve(:,1)
      elseif (a%kfl_ExchangeLocationWallLaw == 1) then
         call e%elmdcg
         call e%elmlen
         delta = e%hleng(1)     ! We set as delta the maximum element length
         vel_wall = 0
         vel_wall_vec = 0.0_rp
         ! We set as velocity the maximum velocity of the element, which corresponds to interior nodes
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
