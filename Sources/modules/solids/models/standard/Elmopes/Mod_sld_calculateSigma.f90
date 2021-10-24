module Mod_sld_calculateSigma
   use typre
   use Mod_sld_BaseElmope
   implicit none
   private
   public SetPointersCalculateSigma
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersCalculateSigma(itask)
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

              !-------------------------------------------------------
              call ConcatenateProcedures(ProcHook%ResetArrays,ResetElemSigma)
              call ConcatenateProcedures(ProcHook%PreGauss,GaussResetSigma)

              if(a%kfl_printSigma) call ConcatenateProcedures(ProcHook%PrePostInterpolates,calculateGPStress)

              if(a%kfl_printGaussSigma) call ConcatenateProcedures(ProcHook%EndGauss,EndElmopeAssemblyGaussSigma_IRR)
              if(a%kfl_printNodalSigma) then
                  call ConcatenateProcedures(ProcHook%Initializations,ResetSigma)
                  call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysSigma)
                  call ConcatenateProcedures(ProcHook%Finalizations,AssembleSigmaEndElmope)
                  call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopeSigmaGp2ElemArray)
              endif

          endif  
      case(100)
          deallocate(kfl_IsSet)
          call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine   
   
   !For actual computations
   !------------------------------------------------------
   subroutine ResetSigma
       implicit none

       a%sigma= 0.0_rp
   end subroutine ResetSigma

   subroutine ResetElemSigma
       implicit none

       elsigma   = 0.0_rp
   end subroutine ResetElemSigma

   subroutine GaussResetSigma
       implicit none

       gpsigma = 0.0_rp
   end subroutine GaussResetSigma

   subroutine EndElmopeAssemblyToArraysSigma
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call a%Mesh%AssemblyToArray(e,tn,elsigma,a%sigma(:,:,1))

   end subroutine EndElmopeAssemblyToArraysSigma

   subroutine EndElmopeAssemblyGaussSigma_IRR
       implicit none
       integer(ip)::tn,nd

       a%sigma_g(ielem)%a(:,igaus) = gpsigma

   end subroutine EndElmopeAssemblyGaussSigma_IRR

   subroutine AssembleSigmaEndElmope
       implicit none
       integer(ip)::ipoin,tn,nd,np
       real(rp), pointer :: aux(:,:)

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do ipoin = 1,np
           if (vmass(ipoin) > 0.0_rp) then 
               a%sigma(:,ipoin,1)    = a%sigma(:,ipoin,1)/vmass(ipoin)
           endif
       enddo

       aux => a%sigma(:,:,1)

       !call a%Mesh%ArrayCommunicator%GhostCommunicate(tn,a%sigma(:,:,1))    
       call a%Mesh%ArrayCommunicator%GhostCommunicate(tn,aux)    

       !Hanging nodes
       if (a%Mesh%kfl_HangingNodes .eqv. .true.) call a%Mesh%InterpolateHangingValues(tn,a%sigma(:,:,1))
       if (a%Mesh%kfl_perio == 1) call a%Mesh%MasterToSlave(tn,a%sigma(:,:,1))

   end subroutine AssembleSigmaEndElmope

   subroutine EndElmopeSigmaGp2ElemArray
       implicit none
       integer(ip)::tn,nd,np,idime,inode

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do inode = 1,e%pnode

           do idime=1,tn
               elsigma(idime,inode)    = elsigma(idime,inode)  + gpsigma(idime)*e%shape(inode,e%igaus)*dvol
           end do

       end do

   end subroutine EndElmopeSigmaGp2ElemArray

end module
