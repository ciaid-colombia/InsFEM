module Mod_lmn_DynSubsElmope
   use typre
   use Mod_lmn_BaseElmope
   implicit none
   private
   public SetPointersDynSubsElmope
   
   integer(ip), allocatable :: kfl_IsSet
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersDynSubsElmope(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            !Tracking of the subscales
            !See the documentation folder
            !Dynamic subscales
            if (a%kfl_tacsg == 1) then
            
               if (a%kfl_repro == 1) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsOSS)
               
               elseif (a%kfl_repro == 0) then
                  call ConcatenateProcedures(ProcHook%InGaussElmats,InGaussElmatsDSS)
               endif
            endif
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   

   !-------------------------------------------------------------------
   !Computations subroutines
   !Dynamic subscales
   
   subroutine InGaussElmatsOSS
      use typre
      implicit none
      !We add the contribution of the subscales at the previous time step 
     
      !Compute contributions to RHS : Block U
      call lmn_elmrhu_dss(e,dvol,ReferenceDtinv,gpden(2),testf_mom,a%vesgs(ielem)%a(:,2,e%igaus),elrhu)
      
      !Compute contributions to RHS : Block P
      call lmn_elmrhp_dss(e,dvol,ReferenceDtinv,acden,gpden(2),timom,a%vesgs(ielem)%a(:,2,e%igaus),elrhp)

      !Compute contributions to RHS : Block T
      call lmn_elmrht_dss(e,dvol,ReferenceDtinv,gpden(2),testf_ene,a%tesgs(ielem)%a(2,e%igaus),elrht)
   end subroutine
   
   subroutine InGaussElmatsDSS
      use typre
      implicit none
      !We add the contribution of the subscales at the previous time step 
     
      !Compute contributions to RHS : Block U
      call lmn_elmrhu_dss(e,dvol,ReferenceDtinv,gpden(2),testf_mom_DSS,a%vesgs(ielem)%a(:,2,e%igaus),elrhu)
      
      !Compute contributions to RHS : Block P
      call lmn_elmrhp_dss(e,dvol,ReferenceDtinv,acden,gpden(2),timom,a%vesgs(ielem)%a(:,2,e%igaus),elrhp)

      !Compute contributions to RHS : Block T
      call lmn_elmrht_dss(e,dvol,ReferenceDtinv,gpden(2),testf_ene_DSS,a%tesgs(ielem)%a(2,e%igaus),elrht)
   end subroutine
      
end module
