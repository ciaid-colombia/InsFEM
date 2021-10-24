module Mod_nsm_ComputeSubscales
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   use Mod_nsm_SubgridSpaceResidual
   use Mod_nsm_ComputeTaus
   implicit none
   private
   public SetPointersComputeSubscales, GPSGS, SetPointersGetSubscales
   
   type, extends(PointerSetter) :: SPComputeSubscales
contains
      procedure :: SpecificSet => SpecificSetComputeSubscales
   end type
   type(SPComputeSubscales) :: SetPointersComputeSubscales

   type, extends(PointerSetter) :: SPGetSubscales
contains
      procedure :: SpecificSet => SpecificSetGetSubscales
   end type
   type(SPGetSubscales) :: SetPointersGetSubscales

   real(rp) :: GpSGS(7)
   
contains
   
   subroutine SpecificSetComputeSubscales(d)
      implicit none
      class(SPComputeSubscales) :: d
      integer(ip) :: kfl_nonlinear
      
     call SetPointersComputeSubgridSpaceResidual%Set
     if (a%kfl_tacsg == 0) then
        call SetPointersComputeTaus%Set
        if(a%kfl_repro == 0 .or. a%kfl_repro == 1)then
           call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleQSS)
        elseif(a%kfl_repro == 2 .or. a%kfl_repro == 3)then
           call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleQSS_Split)
        else
           call runend('dynsubs not ready for kfl_repro not 0, 1, 2 or 3')
        end if
        
     elseif ( a%kfl_tacsg == 1) then
        call SetPointersComputeTaus%Set
        if(a%kfl_repro == 0 .or. a%kfl_repro == 1)then
           call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS)
        elseif(a%kfl_repro == 2 .or. a%kfl_repro == 3)then
           call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeSubgridScaleDSS_Split) 
        else
           call runend('dynsubs not ready for kfl_repro not 0, 1, 2 or 3')
        end if
     endif
   end subroutine   
   
   subroutine SpecificSetGetSubscales(d)
      implicit none
      class(SPGetSubscales) :: d
            
      if (a%kfl_trasg == 0) then
         call SetPointersComputeSubscales%Set
      else 
         call SetPointersComputeSubgridSpaceResidual%Set
         call SetPointersComputeTaus%Set
         if (a%kfl_repro == 0 .or. a%kfl_repro == 1) then
            call ConcatenateProcedures(ProcHook_InGaussElmats,GetSubgridScaleDSS)
         elseif (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
            call ConcatenateProcedures(ProcHook_InGaussElmats,GetSubgridScaleDSS_Split)
         else
            call runend('dynsubs not ready for kfl_repro not 0, 1, 2 or 3')
         end if
      endif
   end subroutine   
   
   !-------------------------------------------------------------------
   !Tracking of Subscales
   !Dynamic subscales
   subroutine ComputeSubgridScaleDSS
      implicit none
      
      !Uses a Backward Euler scheme
      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = timom*(-gpSubscaleSpaceResidual(1:e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      a%prsgs(ielem)%a(e%igaus) = tidiv*gpSubscaleSpaceResidual(e%ndime+1)
   end subroutine
   
   !Dynamic subscales for Split OSS
   subroutine ComputeSubgridScaleDSS_Split
      implicit none
      
      !Uses a Backward Euler scheme
      !note that the dimension of gpSubscaleSpaceResidual is diferent to the dimension of vesgs
      !First subscale (convective component)
      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = timom*(-gpSubscaleSpaceResidual(1:e%ndime) + acden*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
      !Second subscale (pressure gradient component)
      a%vesgs2(ielem)%a(1:e%ndime,1,e%igaus) = timom*(-gpSubscaleSpaceResidual(e%ndime+2:2*e%ndime+1) + acden*a%vesgs2(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
   end subroutine   

   !Get Subscales
   subroutine GetSubgridScaleDSS
      implicit none
      
      GpSGS(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(e%ndime+1) = tidiv*gpSubscaleSpaceResidual(e%ndime+1)
   end subroutine
   
   subroutine GetSubgridScaleDSS_Split
      implicit none
      
      GpSGS(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(e%ndime+1) = tidiv*gpSubscaleSpaceResidual(e%ndime+1)
      GpSGS(e%ndime+2:2*e%ndime+1) = a%vesgs2(ielem)%a(1:e%ndime,1,e%igaus)
   end subroutine   
   
   !Static subscales  
   subroutine ComputeSubgridScaleQSS
      implicit none

      GpSGS(1:e%ndime) = -timom*gpSubscaleSpaceResidual(1:e%ndime)
      GpSGS(e%ndime+1) = tidiv*gpSubscaleSpaceResidual(e%ndime+1)
      if (a%kfl_trasg /= 0) then
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(1:e%ndime)
      endif
   end subroutine
   
   subroutine ComputeSubgridScaleQSS_Split
      implicit none

      GpSGS(1:e%ndime) = -timom*gpSubscaleSpaceResidual(1:e%ndime)
      GpSGS(e%ndime+2:2*e%ndime+1) = -timom*gpSubscaleSpaceResidual(e%ndime+2:2*e%ndime+1)
      GpSGS(e%ndime+1) = tidiv*gpSubscaleSpaceResidual(e%ndime+1)
      if (a%kfl_trasg /= 0) then
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(1:e%ndime)
         a%vesgs2(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(e%ndime+2:2*e%ndime+1)
      endif
   end subroutine  
   
end module 
