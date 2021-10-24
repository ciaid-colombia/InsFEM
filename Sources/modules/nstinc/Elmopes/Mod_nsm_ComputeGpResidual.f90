  module Mod_nsm_ComputeGpResidual
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_InterpolateGradients
   use Mod_CutMesh   
   implicit none
   private
   public SetPointersComputeGpResidual, gpres
   
   type, extends(PointerSetter) :: SPComputeGpResidual
contains
      procedure :: SpecificSet => SpecificSetComputeGpResidual
   end type
   type(SPComputeGpResidual) :: SetPointersComputeGpResidual
   
   real(rp), allocatable    :: gpres(:)

contains
   
   subroutine SpecificSetComputeGpResidual(d)
      implicit none
      class(SPComputeGpResidual) :: d
      integer(ip) :: kfl_nonlinear
      
      !We need the gradients in the gauss point
      call SetPointersInterpolateGradients%Set
      call ConcatenateProcedures(ProcHook_Initializations,AllocGpRes)
      call ConcatenateProcedures(ProcHook_Finalizations,DeallocGpRes)

      if (a%kfl_repro == 0 .or. a%kfl_repro == 1) then   
         if (a%kfl_repro == 0 .or. a%kfl_repro_SkipFE == 0) then
            call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpRes)
         elseif (a%kfl_repro == 1 .and. a%kfl_repro_SkipFE == 1) then
            call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpRes_SkipFE)
         elseif(a%kfl_repro == 1 .and. a%kfl_repro_SkipFE == 2) then
            call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResSkipt)            
         else
            call runend('nsm_EndsteElmope: Dissipation option not implemented?')
         endif
         call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
         if (a%kfl_repro == 1 .and. kfl_nonlinear == 1) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpNonLinearRes)
         !For enriched two phase flows
         if (a%kfl_colev==1 .and. a%kfl_EnrichElem==1 .and. a%kfl_repro==1 .and. a%kfl_repro_SkipFE==0) then 
            call PrependProcedure(ProcHook_InGaussElmatsAssembly,Newgpres)               
         end if
         !For FOR Rotating-Axis (Coriolis term)
         if (a%kfl_FORAxesRotation == 1) call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpCoriolisRes)
      elseif (a%kfl_repro == 2) then
         !Split OSS
         call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResSplit)        
         if(a%kfl_colev==1 .and. a%kfl_EnrichElem==1) then 
            call PrependProcedure(ProcHook_InGaussElmatsAssembly,Newgpres)               
         end if
      elseif (a%kfl_repro == 3) then
         !Split OSS, with gravity (free surface case or Boussinesq case)
         call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResSplit2)
         
         if(a%kfl_colev==1 .and. a%kfl_EnrichElem==1)then 
            call PrependProcedure(ProcHook_InGaussElmatsAssembly,Newgpres)               
         end if
      endif
      !Store residual for postprocess
      if (a%npp_stepi(18) /= 0) then
         call ConcatenateProcedures(ProcHook_InGaussElmats,StoreGpRes)
      endif
      
   end subroutine
 
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !-------------------------------------------------------------------

   subroutine AllocGpRes
      implicit none
      
      call a%Memor%alloc(a%ResidualSize,gpres,'gpres','nsm_EnditeElmope')
   end subroutine
   
   subroutine DeAllocGpRes
      implicit none
      
      call a%Memor%dealloc(a%ResidualSize,gpres,'gpres','nsm_EnditeElmope')
   end subroutine
   
   subroutine ComputeGpRes
      implicit none

      call nsm_elmrfe_oto(e,LHSdtinv,acden,gpadv,gpvel,grpre,grvel,elext,eltemp,gpres)
   end subroutine
   
   subroutine ComputeGpResSkipt
      implicit none
      real(rp) :: auxtemp(e%ndime) 
      
      auxtemp = 0.0_rp
      call nsm_elmrfe_oto(e,0.0_rp,acden,gpadv,gpvel,grpre,grvel,elext,auxtemp,gpres)      
   end subroutine
   
   subroutine ComputeGpResSplit
      implicit none

      call nsm_elmrfe_split(e,acden,gpadv,grpre,grvel,gpres)
   end subroutine
   
   subroutine ComputeGpResSplit2
      implicit none
      integer(ip) :: idime,jdime

      call nsm_elmrfe_split(e,acden,gpadv,grpre,grvel,gpres)
      gpres(e%ndime+2:2*e%ndime+1) = gpres(e%ndime+2:2*e%ndime+1) - elext(1:e%ndime)
   end subroutine
   
   subroutine ComputeGpNonLinearRes
      implicit none
      
      call nsm_elmrfe_oto_nonlinear(e,acvis,elvel,gpres)
   end subroutine
   
   subroutine ComputeGpRes_SkipFE
      implicit none

      call nsm_elmrfe_trm(e,acden,gpadv,grpre,grvel,gpres)
   end subroutine
   
   subroutine ComputeGpCoriolisRes
      implicit none
      
      call nsm_elmrfe_oto_coriolis(e,acden,a%FORAxesAngularVeloc,gpvel,gpres)
   end subroutine

   subroutine Newgpres
      implicit none
      integer(ip) :: elemStatus
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      if(elemStatus==0)then      
         gpres=0.0_rp
      end if
   end subroutine
   
   subroutine Newgpres2
      implicit none
      integer(ip) :: elemStatus
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      if(elemStatus<=0)then      
         gpres=0.0_rp
      end if
   end subroutine
   
   subroutine StoreGpRes
      implicit none
      
      !Protect if cut elements and more gaus points, representation is inaccurate
      if (e%igaus <= size(a%ResidualU(ielem)%a,2)) then      
         a%residualU(ielem)%a(:,e%igaus) = gpres(1:e%ndime)
         a%residualP(ielem)%a(e%igaus)   = gpres(e%ndime+1)
         if (a%kfl_repro == 2) a%ResidualGraP2(ielem)%a(:,e%igaus) = gpres(e%ndime+2:2*e%ndime+1)
      endif
   end subroutine
   
end module
