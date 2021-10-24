module Mod_sldsup_EndElmope_lin
   use typre
   use Mod_sldsup_EndElmope
   implicit none
   
contains       
    
   subroutine SetPointersSUPEndElmope_lin
      implicit none
      
      integer(ip) ::  nelty
      
      call ConcatenateProcedures(ProcHook%AllocateArrays,AllocateBaseElmopeMatricesSUPLinear)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUPLinear)

      call ConcatenateProcedures(ProcHook%PhysicalProp,PhysicalProp_lin)

      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChangeSUPEndElmope_lin

   end subroutine

  subroutine PhysicalProp_lin
      implicit none
      integer(ip) :: nd,tn

      call sup_lin%Mesh%GetNdime(nd)
      tn=(nd*(nd+1))/2
      
      call sup_lin%SUPGetPhysicalParameters(nd,tn,densi,K,G,C_dev,D,D_dev,Dv_scalar,P)
  end subroutine

   subroutine EndElmopeCalculateDevStrain_lin
       implicit none
       integer(ip)::tn,nd,idime,inode

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       !----------------------------------------------------------------
       !Strain value
       !  e=D*S  S=depends on constitutive model, D=inv(C)
       do idime =1,tn
            dstrain(idime) = dot_product(D_dev(:,idime),gpsigma(:))
       end do

       !Now we fill strain vector
       do inode = 1,e%pnode

            do idime=1,tn
               devstrain(idime,inode) = devstrain(idime,inode) + dstrain(idime)*e%shape(inode,e%igaus)*dvol
            end do

       end do


   end subroutine EndElmopeCalculateDevStrain_lin

   subroutine EndElmopeCalculateStrain_lin
       implicit none
       integer(ip)::tn,nd,idime,inode

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       !Get Gauss Voigt strain
       do idime =1,nd
           stress(idime)  = gpsigma(idime) + gppress(1_ip)
       enddo
       stress(nd+1:tn) = gpsigma(nd+1:tn)

       !----------------------------------------------------------------
       !Strain value
       !  e=D*S  S=depends on constitutive model, D=inv(C)
       do idime =1,tn
           strain(idime)  = dot_product(D(:,idime),stress(:))
       end do

       !Now we fill strain vector
       do inode = 1,e%pnode

           do idime=1,tn
               estrain(idime,inode)   = estrain(idime,inode)   + strain(idime)*e%shape(inode,e%igaus)*dvol
               estress(idime,inode)   = estress(idime,inode)   + stress(idime)*e%shape(inode,e%igaus)*dvol
           end do

       end do

   end subroutine EndElmopeCalculateStrain_lin

   subroutine SetPointersCalculateStrain_lin
       implicit none

      call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopeCalculateStrain_lin)

   end subroutine SetPointersCalculateStrain_lin

   subroutine SetPointersCalculateDevStrain_lin
       implicit none

      call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopeCalculateDevStrain_lin)

   end subroutine SetPointersCalculateDevStrain_lin

   !--------------------------------------------------------------------
   !Multiple type of elements
   subroutine OnIeltyChangeSUPEndElmope_lin
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointersSUPEndElmope_lin
         ielty0 = e%ielty
      endif
   end subroutine
   
end module

subroutine sldsup_EndElmope_lin(SldProblem,task)
   use Mod_sldsup_EndElmope_lin
   use Mod_SUPSolids_lin
   implicit none
   class(SUPSolidsProblem_lin), target :: SldProblem
   character(6) :: task
   integer(ip)  :: tn,npoin,inode
   real(rp)     :: pstr(3)

   itask = task

   a       => SldProblem
   sup     => SldProblem
   sup_lin => SldProblem

   call SetPointersGeneral
   call SetPointersSUPEndElmope
   call SetPointersSUPEndElmope_lin

   !Things to be done if endite
   if (itask .eq. 'Endite') then

       sup%devstrain = 0.0_rp

       !To be able to calculate devnorm for tau_u
       call SetPointersCalculateDevStrain_lin
       call SetPointersAssembleDevStrain

       call SetPointersCalculateVmass
       call SetPointersAssembleVmass

       !Things to be done if endste
   elseif (itask .eq. 'Endste') then

       if(sup_lin%kfl_printStrain .or. sup_lin%kfl_printStress) then

           call SetPointersCalculateStrain_lin

           if(sup_lin%kfl_printStress) then
               if(a%kfl_GaussStress) then
                   call SetPointersAssembleGaussStress
               endif
               if(a%kfl_NodalStress) then
                   sup%stress = 0.0_rp
                   call SetPointersAssembleStress
               endif

           endif


           if(sup_lin%kfl_printStrain) then
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

       if(sup_lin%kfl_printJ2Stresses) call SetPointersJ2Stresses

   end if

   call sld_elemLoop

   if(sup%kfl_useSecantModulus) then
       call a%Mesh%GetNpoin(npoin)
       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call vecnorALLMPI(sup%sigma(:,:,1),tn,npoin,sup%signorm,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)
       call vecnorALLMPI(sup%devstrain(:,:),tn,npoin,sup%devnorm,2,a%MPIcomm,a%MPIrank,a%MPIroot,a%MPIsize)
   endif

end subroutine
