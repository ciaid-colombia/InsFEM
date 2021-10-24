module Mod_sldsup_EndElmope_NH
   use typre
   use Mod_sldup_EndElmope
   use Mod_sldsup_ComputeGpResidual
   use Mod_sldsup_ComputeResidualProjection
   use Mod_sldsup_ComputeSubscales
   use Mod_sldsup_SubgridSpaceResidual
   use Mod_sldsup_InterpolateResidualProjection
   implicit none

contains       
    
   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointersSUPEndElmope_NH
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty
      
      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateNonLinearSolidArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateNonLinearSolidArrays)

      ProcPointer%getTauParameters => getTauParameters

      call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp_NH)

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChangeSUPEndElmope_NH
   end subroutine SetPointersSUPEndElmope_NH

   subroutine SetPointersSUPEndElmope_NH_REPRO

      call SetPointersComputeGpResidual(0)
      call SetPointersComputeResidualProjection(0)
      call SetPointersComputeSubscales(0)              
      call SetPointersComputeSubgridSpaceResidual(0)
      call SetPointersInterpolateResidualProjection(0)

      if (sup_NH%kfl_trasg == 1 .or. sup_NH%kfl_tacsg == 1) then
         call SetPointersComputeSubscales(1)
      endif

      if (sup_NH%kfl_repro == 1) call SetPointersComputeResidualProjection(1)

      call SetPointersInterpolateResidualProjection(100)
      call SetPointersComputeGpResidual(100)
      call SetPointersComputeResidualProjection(100)
      call SetPointersComputeSubscales(100)              
      call SetPointersComputeSubgridSpaceResidual(100)

   end subroutine SetPointersSUPEndElmope_NH_REPRO

  subroutine PhysicalProp_NH
      implicit none
      integer(ip) :: nd,tn

      call sup_NH%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call sup_NH%SUPGetPhysicalParameters(nd,tn,densi,lam,G)
  end subroutine

   subroutine SetPointersCalculateStrain_NH
       implicit none

      call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopeCalculateStrain_NH)

   end subroutine SetPointersCalculateStrain_NH

   subroutine EndElmopeCalculateStrain_NH
       implicit none
       integer(ip)::tn,nd,np,idime,inode

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       !Get Gauss Voigt strain
       do idime =1,nd
           stress(idime)  = gpsigma(idime) + gppress(1_ip)
       enddo
       stress(nd+1:tn) = gpsigma(nd+1:tn)

       !----------------------------------------------------------------
       !Strain value
       !  e=D*S  S=depends on constitutive model, D=inv(C)
       !do idime =1,tn
       !    strain(idime)  = dot_product(D(:,idime),stress(:))
       !end do

       !Now we fill strain vector
       do inode = 1,e%pnode

           do idime=1,tn
               estrain(idime,inode)   = estrain(idime,inode)   + strain(idime)*e%shape(inode,e%igaus)*dvol
               estress(idime,inode)   = estress(idime,inode)   + stress(idime)*e%shape(inode,e%igaus)*dvol
           end do

       end do

   end subroutine EndElmopeCalculateStrain_NH

   !--------------------------------------------------------------------
   !Multiple type of elements
   subroutine OnIeltyChangeSUPEndElmope_NH
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointersSUPEndElmope_NH
         ielty0 = e%ielty
      endif

   end subroutine
   
end module

subroutine sldsup_EndElmope_NH(SldProblem,task)
   use Mod_sldsup_EndElmope_NH
   use Mod_SUPSolids_NH
   implicit none
   class(SUPSolidsProblem_NH), target :: SldProblem
   character(6) :: task
   integer               :: ipoin,tn,idime,ndime,inode,npoin
   real(rp)              :: strmodu,sigmodu

   itask = task

   a      => SldProblem
   up    => SldProblem
   sup    => SldProblem
   sup_NH => SldProblem

   call SetPointersGeneral
   call SetPointersSUPEndElmope
   call SetPointersSUPEndElmope_NH

   !Things to be done if endite
   if (itask .eq. 'Endite') then

   !Things to be done if endste
   elseif (itask .eq. 'Endste') then

       call SetPointersSUPEndElmope_NH_REPRO

       if(sup_NH%kfl_printStrain .or. sup_NH%kfl_printStress) then

           call SetPointersCalculateStrain_NH

           if(sup_NH%kfl_printStress) then
               if(a%kfl_GaussStress) then
                   call SetPointersAssembleGaussStress
               endif
               if(a%kfl_NodalStress) then
                   sup%stress = 0.0_rp
                   call SetPointersAssembleStress
               endif

           endif


           if(sup_NH%kfl_printStrain) then
               if(a%kfl_GaussStrain) then
                   call SetPointersAssembleGaussStrain

               endif
               if(a%kfl_NodalStrain) then
                   sup%strain = 0.0_rp
                   call SetPointersAssembleStrain
               endif

           endif 

           if(a%kfl_NodalStrain .or. a%kfl_NodalStress) then
               call SetPointersCalculateVmass
               call SetPointersAssembleVmass
           endif
       end if

       if(sup_NH%kfl_printJ2Stresses) call SetPointersJ2Stresses

   end if

   call sld_elemLoop


end subroutine
