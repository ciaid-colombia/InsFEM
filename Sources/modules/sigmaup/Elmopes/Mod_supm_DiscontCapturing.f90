module Mod_supm_DiscontCapturing
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateGradients
   implicit none
   private
   public SetPointersDiscontCapturing
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersDiscontCapturing(itask,task)
      integer(ip) :: itask
      character(6) :: task
      !procedure() :: NULL()
      
      select case (itask)   
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            if(a%MatProp(imat)%lawvi<0 .and. a%kfl_shock == 1)then
               if (task .eq. 'Elmope') then             
                  call ConcatenateProcedures(ProcHook_Initializations,AllocRepGrad)
                  call ConcatenateProcedures(ProcHook_Gathers,GatherRepGrad)
                  call ConcatenateProcedures(ProcHook_Interpolates,InterpolateRepGrad)
                  call ConcatenateProcedures(ProcHook_Interpolates,GradRespro)         
                  call ConcatenateProcedures(ProcHook_InGaussElmats,InGaussElmatsRepGrad)
                  call ConcatenateProcedures(ProcHook_Finalizations,DeallocRepGrad)
               else if (task .eq. 'Endite') then
                  call ConcatenateProcedures(ProcHook_PreLoop,PreLoopDiscCapt)
                  call ConcatenateProcedures(ProcHook_Initializations,InitDiscontinuity)
                  call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroRepGrad)
                  call SetPointersInterpolateGradients(1)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,GradientToVector)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,GaussPointAssemblyRepGrad)         
                  call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyGradient)
                  call ConcatenateProcedures(ProcHook_Finalizations,FinDiscontinuity)
                  call ConcatenateProcedures(ProcHook_PostLoop,SmoothGradient)
               end if
            end if
         end if
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select   
      
   end subroutine   
   
   !SUBROUTINES ---------------------------------------
   
   subroutine AllocRepGrad
      implicit none
      cshock(1)=a%shock(1)
      cshock(2)=a%shock(2)
      auxGrad=e%ndime*auxtens
      call a%Memor%alloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supm_elmope')
      call a%Memor%alloc(auxGrad,gpgrad,'gpgrad','supm_elmope')
      call a%Memor%alloc(auxtens,e%ndime,grsigRP,'grsigRP','supm_elmope')
      call a%Memor%alloc(auxtens,e%ndime,grsigRPO,'grsigRPO','supm_elmope')    
   end subroutine

   subroutine DeallocRepGrad
      implicit none
      call a%Memor%dealloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supm_elmope')
      call a%Memor%dealloc(auxGrad,gpgrad,'gpgrad','supm_elmope')
      call a%Memor%dealloc(auxtens,e%ndime,grsigRP,'grsigRP','supm_elmope')      
      call a%Memor%dealloc(auxtens,e%ndime,grsigRPO,'grsigRPO','supm_elmope')       
   end subroutine
   
   subroutine GatherRepGrad
      implicit none
      call e%gather(auxGrad,elrepGrad,a%reproGrad)    
   end subroutine

   subroutine InterpolateRepGrad
      implicit none
      call e%interpg(auxGrad,elrepGrad,gpgrad)
   end subroutine
   
   subroutine GradRespro
      implicit none
      real(rp) :: gradsig(auxtens,e%ndime)
      call sup_grsigRP(e%ndime,auxtens,auxGrad,gpgrad,grsigRP)      
      grsigRPO=0.0_rp
      !Orthogonal Projection
      if (a%LogFormulation==0) gradsig=grsig
      if (a%LogFormulation==1) gradsig=GrExpPsi*auxL

      do kdime=1,auxtens
         do ldime=1,e%ndime
            grsigRPO(kdime,ldime)= grsigRPO(kdime,ldime) + (gradsig(kdime,ldime)-grsigRP(kdime,ldime))
         end do
      enddo

   end subroutine
   
   subroutine InGaussElmatsRepGrad
      use Mod_sup_Kdiscont
      implicit none
      real(rp) :: gradsig(auxtens,e%ndime)
      if (a%LogFormulation==0) gradsig=grsig
      if (a%LogFormulation==1) gradsig=GrExpPsi*auxL
      
      call sup_Kdiscont(e,cshock,auxtens,acvis,lambda,chale,gradsig,grsigRPO,grsigRP,grvel,gpvno,facdisc,kdisc)    
      if(a%kfl_colev==0)then
         a%viscarray(ielem)%a(e%igaus) = Kdisc
      end if
      call  ProcPointer%DiscontinuityCapturing
   end subroutine
   
   subroutine discontinuity2d
      implicit none   
      if (a%LogFormulation==0) then
         call supm_elmbstDC(e,auxtens,kdisc,dvol,elmst)   
      else if (a%LogFormulation==1) then
         call supm_elmbstDC_LCR(e,auxtens,kdisc,auxL,dvol,ExpGpPsi_Matrix,GrExpMatrix,elmst)
         call supm_elmrhcDC_LCR(e,kdisc,auxL,dvol,gpsig,auxtens,ExpGpPsi_Matrix,GrPsiMatrix,GrExpMatrix,elrhc)
      end if
   end subroutine
   
   subroutine PreLoopDiscCapt
      implicit none
      a%reproGrad = 0.0_rp
   end subroutine   
   
   subroutine InitDiscontinuity
      implicit none   
      auxGrad=e%ndime*auxtens
      !Matrices alloc
      call a%Memor%alloc(auxGrad,gpgrad,'gpgrad','supm_EnditeElmope')
      call a%Memor%alloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supm_EnditeElmope')         
   end subroutine
   
   subroutine ElmatsToZeroRepGrad
      implicit none
      elrepGrad = 0.0_rp
   end subroutine   

   subroutine FinDiscontinuity
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(auxGrad,gpgrad,'gpgrad','supm_EnditeElmope')
      call a%Memor%dealloc(auxGrad,e%mnode,elrepGrad,'elrepGrad','supm_EnditeElmope')               
   end subroutine
      
   subroutine GradientToVector
      implicit none
      real(rp) :: gradsig(auxtens,e%ndime)
      
      if (a%LogFormulation==0) gradsig=grsig
      if (a%LogFormulation==1) gradsig=GrExpPsi*auxL
      
      call sup_GradStress(e%ndime,auxtens,auxGrad,gradsig,gpgrad)

   end subroutine
   
   subroutine GaussPointAssemblyRepGrad
      implicit none   
      call supm_elmrep(e,dvol,auxGrad,gpgrad,elrepGrad)
   end subroutine   
   
   subroutine AssemblyGradient
      implicit none
      a%reproGrad(:,e%lnods(1:e%pnode)) = a%reproGrad(:,e%lnods(1:e%pnode)) + elrepGrad(:,1:e%pnode)
   end subroutine
   
   subroutine SmoothGradient
      implicit none
      call a%Project(auxGrad,a%reproGrad) 
   end subroutine 
   
end module
   
   