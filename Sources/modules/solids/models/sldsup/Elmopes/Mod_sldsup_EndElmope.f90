module Mod_sldsup_EndElmope
   use typre
   use Mod_CauchyElement
   use Mod_sld_BaseElmope
   use Mod_SUPSolids
   use Mod_sld_BaseElmopeRoutines
   implicit none

   real(rp), allocatable :: estrain(:,:),dstrain(:)
   real(rp), allocatable :: estress(:,:),esigma(:,:)
   real(rp), allocatable :: devstrain(:,:)
   
contains       
    
   !-------------------------------------------------------------------
   !SetPointers
   subroutine SetPointersGeneralUPEndElmope
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty
      
      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeArrays)

      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateBaseElmopeMatricesSUP)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateBaseElmopeMatricesSUP)

      call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateEndElmopeArrays)
      call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateEndElmopeArrays)

      call ConcatenateProcedures(ProcHook%ResetArrays,ResetEndElmopeArrays)

      call ConcatenateProcedures(ProcHook%PreGauss,GaussResetEndElmopeArrays)

      call SetPointersExternalForces(0)
      call SetPointersExternalForces(1)
      call SetPointersExternalForces(100)
      call ConcatenateProcedures(ProcHook%Initializations,calculateExternalForces)

      call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%ComputeResidual)
      call ConcatenateProcedures(ProcHook%PostInterpolates,ProcHook%InGaussElmats)

   end subroutine SetPointersGeneralUPEndElmope

   subroutine SetPointersUPEndElmope
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty

      call SetPointersGeneralUPEndElmope
      
      call ConcatenateProcedures(ProcHook%Gathers,displacementGather)
      call ConcatenateProcedures(ProcHook%Gathers,pressGather)

      call ConcatenateProcedures(ProcHook%PrePostInterpolates,calculateUPGradients)
      call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpPress)

   end subroutine SetPointersUPEndElmope


   subroutine SetPointersSUPEndElmope
      implicit none
      
      integer(ip) :: kfl_nonlinear, nelty
      
      call SetPointersGeneralUPEndElmope

      call ConcatenateProcedures(ProcHook%Gathers,displacementGather)
      call ConcatenateProcedures(ProcHook%Gathers,sigmaGather)
      call ConcatenateProcedures(ProcHook%Gathers,pressGather)

      call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpSigma)
      call ConcatenateProcedures(ProcHook%Interpolates,InterpolateGpPress)

   end subroutine SetPointersSUPEndElmope

  subroutine calculateExternalForces
      implicit none

      !Compute contributions to RHS :
      elext = 0.0_rp
      call ProcPointer%ExternalForces   

  end subroutine

   subroutine AssembleDevStrainEndElmope
       implicit none
       integer(ip)::ipoin,tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do ipoin = 1,np
           if (vmass(ipoin) > 0.0_rp) then 
               sup%devstrain(:,ipoin) = sup%devstrain(:,ipoin)/vmass(ipoin)
           endif
       enddo
       call sup%Mesh%ArrayCommunicator%GhostCommunicate(tn,sup%devstrain)    

       !Hanging nodes
       if (a%Mesh%kfl_HangingNodes .eqv. .true.) call sup%Mesh%InterpolateHangingValues(tn,sup%devstrain)
       if (a%Mesh%kfl_perio == 1) call sup%Mesh%MasterToSlave(tn,sup%devstrain)

   end subroutine AssembleDevStrainEndElmope

   subroutine AssembleStressEndElmope
       implicit none
       integer(ip)::ipoin,tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do ipoin = 1,np
           if (vmass(ipoin) > 0.0_rp) then 
               sup%stress(:,ipoin)    = sup%stress(:,ipoin)/vmass(ipoin)
           endif
       enddo
       call sup%Mesh%ArrayCommunicator%GhostCommunicate(tn,sup%stress)    

       !Hanging nodes
       if (a%Mesh%kfl_HangingNodes .eqv. .true.) call sup%Mesh%InterpolateHangingValues(tn,sup%stress)
       if (a%Mesh%kfl_perio == 1) call sup%Mesh%MasterToSlave(tn,sup%stress)

   end subroutine AssembleStressEndElmope

   subroutine AssembleSigmaEndElmope
       implicit none
       integer(ip)::ipoin,tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do ipoin = 1,np
           if (vmass(ipoin) > 0.0_rp) then 
               sup%sigma(:,ipoin,1)    = sup%sigma(:,ipoin,1)/vmass(ipoin)
           endif
       enddo
       call sup%Mesh%ArrayCommunicator%GhostCommunicate(tn,sup%sigma(:,:,1))    

       !Hanging nodes
       if (a%Mesh%kfl_HangingNodes .eqv. .true.) call sup%Mesh%InterpolateHangingValues(tn,sup%sigma(:,:,1))
       if (a%Mesh%kfl_perio == 1) call sup%Mesh%MasterToSlave(tn,sup%sigma(:,:,1))

   end subroutine AssembleSigmaEndElmope

   subroutine AssembleStrainEndElmope
       implicit none
       integer(ip)::ipoin,tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do ipoin = 1,np
           if (vmass(ipoin) > 0.0_rp) then 
               sup%strain(:,ipoin)    = sup%strain(:,ipoin)/vmass(ipoin)
           endif
       enddo
       call sup%Mesh%ArrayCommunicator%GhostCommunicate(tn,sup%strain)    

       !Hanging nodes
       if (a%Mesh%kfl_HangingNodes .eqv. .true.) call sup%Mesh%InterpolateHangingValues(tn,sup%strain)
       if (a%Mesh%kfl_perio == 1) call sup%Mesh%MasterToSlave(tn,sup%strain)

   end subroutine AssembleStrainEndElmope

   subroutine ResetEndElmopeArrays
       implicit none

       devstrain = 0.0_rp
       estrain   = 0.0_rp
       estress   = 0.0_rp
       elvmass   = 0.0_rp
       elsigma   = 0.0_rp
       elpress   = 0.0_rp

   end subroutine ResetEndElmopeArrays

   subroutine GaussResetEndElmopeArrays
       implicit none

       dstrain = 0.0_rp
       strain  = 0.0_rp
       stress  = 0.0_rp
       gpsigma = 0.0_rp
       gppress = 0.0_rp

   end subroutine GaussResetEndElmopeArrays

   subroutine SetPointersAssembleVmass
       implicit none

      call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysVmass)

   end subroutine SetPointersAssembleVmass

   subroutine SetPointersAssembleGaussStrain
       implicit none

       call ConcatenateProcedures(ProcHook%EndGauss,EndElmopeAssemblyGaussStrain)

   end subroutine SetPointersAssembleGaussStrain

   subroutine SetPointersAssembleGaussStress
       implicit none

       call ConcatenateProcedures(ProcHook%EndGauss,EndElmopeAssemblyGaussStress)

   end subroutine SetPointersAssembleGaussStress

   subroutine SetPointersAssembleGaussSigma
       implicit none

       call ConcatenateProcedures(ProcHook%EndGauss,EndElmopeAssemblyGaussSigma)

   end subroutine SetPointersAssembleGaussSigma

   subroutine SetPointersAssembleStrain
       implicit none

      call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysStrain)
      call ConcatenateProcedures(ProcHook%Finalizations,AssembleStrainEndElmope)

   end subroutine SetPointersAssembleStrain

   subroutine SetPointersAssembleStress
       implicit none

      call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysStress)
      call ConcatenateProcedures(ProcHook%Finalizations,AssembleStressEndElmope)

   end subroutine SetPointersAssembleStress

   subroutine SetPointersAssembleSigma
       implicit none

      call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysSigma)
      call ConcatenateProcedures(ProcHook%Finalizations,AssembleSigmaEndElmope)

   end subroutine SetPointersAssembleSigma

   subroutine SetPointersAssembleDevStrain
       implicit none

       call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysDevStrain)
       call ConcatenateProcedures(ProcHook%Finalizations,AssembleDevStrainEndElmope)

   end subroutine SetPointersAssembleDevStrain

   subroutine EndElmopeAssemblyToArraysVmass
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call sup%Mesh%AssemblyToArray(e,1_ip,elvmass,vmass)

   end subroutine EndElmopeAssemblyToArraysVmass

   subroutine EndElmopeAssemblyToArraysDevStrain
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call sup%Mesh%AssemblyToArray(e,tn,devstrain,sup%devstrain)

   end subroutine EndElmopeAssemblyToArraysDevStrain

   subroutine EndElmopeAssemblyGaussStrain
       implicit none
       integer(ip)::tn,nd

       sup%strain_g(ielem)%a(:,igaus) = strain

   end subroutine EndElmopeAssemblyGaussStrain

   subroutine EndElmopeAssemblyGaussSigma
       implicit none
       integer(ip)::tn,nd

       sup%sigma_g(ielem)%a(:,igaus) = gpsigma

   end subroutine EndElmopeAssemblyGaussSigma

   subroutine EndElmopeAssemblyGaussStress
       implicit none
       integer(ip)::tn,nd

       sup%stress_g(ielem)%a(:,igaus) = stress

   end subroutine EndElmopeAssemblyGaussStress

   subroutine EndElmopeAssemblyToArraysStress
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call sup%Mesh%AssemblyToArray(e,tn,estress  ,sup%stress   )

   end subroutine EndElmopeAssemblyToArraysStress

   subroutine EndElmopeAssemblyToArraysSigma
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call sup%Mesh%AssemblyToArray(e,tn,esigma,sup%sigma)

   end subroutine EndElmopeAssemblyToArraysSigma

   subroutine EndElmopeAssemblyToArraysStrain
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call sup%Mesh%AssemblyToArray(e,tn,estrain  ,sup%strain   )

   end subroutine EndElmopeAssemblyToArraysStrain

   subroutine EndElmopeCalculateVmass
       implicit none
       integer(ip)::inode

       do inode = 1,e%pnode

            elvmass(inode) = elvmass(inode) + e%shape(inode,e%igaus)*dvol

       end do

   end subroutine EndElmopeCalculateVmass

   subroutine SetPointersCalculateVmass
       implicit none
       integer(ip)::inode

       call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopeCalculateVmass)

   end subroutine SetPointersCalculateVmass

   subroutine calculateJ2Stresses
       implicit none
       integer(ip)::inode,tn,nd,np
       real(rp)   ::pstr(3)

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       sup%J2 = 0.0_rp
       do inode=1,np
          pstr = 0.0_rp
          pstr(1:nd) = getPrincipalStrains(nd,tn,sup%sigma(:,inode,1))
          sup%J2(inode) = 0.5_rp*(pstr(1)*pstr(1) + pstr(2)*pstr(2) + pstr(3)*pstr(3))
       enddo

   end subroutine calculateJ2Stresses

   subroutine SetPointersJ2Stresses
       implicit none
       integer(ip)::inode

       call ConcatenateProcedures(ProcHook%Finalizations,calculateJ2Stresses)

   end subroutine SetPointersJ2Stresses

   subroutine AllocateEndElmopeArrays
       implicit none
       integer(ip)::tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       !Allocation
       call a%Memor%alloc(tn     ,strain           ,'strain'   ,'sldsup_EndElmope')   
       call a%Memor%alloc(tn     ,dstrain          ,'dstrain'  ,'sldsup_EndElmope')   
       call a%Memor%alloc(tn     ,stress           ,'stress'   ,'sldsup_EndElmope')
       call a%Memor%alloc(np     ,vmass            ,'vmass'    ,'sldsup_EndElmope')
       call a%Memor%alloc(e%mnode,elvmass          ,'elvmass'  ,'sldsup_EndElmope')
       call a%Memor%alloc(tn     ,e%mnode,estrain  ,'estrain'  ,'sldsup_EndElmope')   
       call a%Memor%alloc(tn     ,e%mnode,estress  ,'estress'  ,'sldsup_EndElmope')   
       call a%Memor%alloc(tn     ,e%mnode,devstrain  ,'devstrain'  ,'sldsup_EndElmope')   

   end subroutine
   
   subroutine DeallocateEndElmopeArrays
       implicit none
       integer(ip)::tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       !Allocation
       call a%Memor%dealloc(tn     ,strain           ,'strain'   ,'sldsup_EndElmope')   
       call a%Memor%dealloc(tn     ,dstrain          ,'dstrain'  ,'sldsup_EndElmope')   
       call a%Memor%dealloc(tn     ,stress           ,'stress'   ,'sldsup_EndElmope')
       call a%Memor%dealloc(np     ,vmass            ,'vmass'    ,'sldsup_EndElmope')
       call a%Memor%dealloc(e%mnode,elvmass          ,'elvmass'  ,'sldsup_EndElmope')
       call a%Memor%dealloc(tn     ,e%mnode,estrain  ,'estrain'  ,'sldsup_EndElmope')   
       call a%Memor%dealloc(tn     ,e%mnode,estress  ,'estress'  ,'sldsup_EndElmope')   
       call a%Memor%dealloc(tn     ,e%mnode,devstrain  ,'devstrain'  ,'sldsup_EndElmope')   

   end subroutine

end module
