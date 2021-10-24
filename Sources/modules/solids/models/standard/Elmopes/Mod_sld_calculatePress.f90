module Mod_sld_calculatePress
   use typre
   use Mod_sld_BaseElmope
   implicit none
   private
   public SetPointersCalculatePress
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersCalculatePress(itask)
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
              call ConcatenateProcedures(ProcHook%ResetArrays,ResetElemPress)
              call ConcatenateProcedures(ProcHook%PreGauss,GaussResetPress)

              if(a%kfl_printPress) call ConcatenateProcedures(ProcHook%PrePostInterpolates,calculateGPPress)

              if(a%kfl_printGaussPress) call ConcatenateProcedures(ProcHook%EndGauss,EndElmopeAssemblyGaussPress)
              if(a%kfl_printNodalPress) then
                  call ConcatenateProcedures(ProcHook%Initializations,ResetPress)
                  call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysPress)
                  call ConcatenateProcedures(ProcHook%Finalizations,AssemblePressEndElmope)
                  call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopePressGp2ElemArray)
              endif
          endif  
      case(100)
          deallocate(kfl_IsSet)
          call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine   
   
   !For actual computations
   !------------------------------------------------------
   subroutine ResetPress
       implicit none

       a%press = 0.0_rp
   end subroutine ResetPress

   subroutine ResetElemPress
       implicit none

       elpress   = 0.0_rp
   end subroutine ResetElemPress

   subroutine GaussResetPress
       implicit none

       gppress= 0.0_rp
   end subroutine GaussResetPress

   subroutine EndElmopeAssemblyToArraysPress
       implicit none
       integer(ip) :: inode

       !Has memory fault 
       !call a%Mesh%AssemblyToArray(e,1,elpress(1,:),a%press(:,1))

       do inode = 1,e%pnode
         a%press(e%lnods(inode),1) = a%press(e%lnods(inode),1) + elpress(1,inode)
       enddo

   end subroutine EndElmopeAssemblyToArraysPress

   subroutine EndElmopeAssemblyGaussPress
       implicit none
       integer(ip)::tn,nd

       a%press_g(ielem)%a(igaus) = gppress(1_ip)

   end subroutine EndElmopeAssemblyGaussPress

   subroutine AssemblePressEndElmope
       implicit none
       integer(ip)::ipoin,tn,nd,np
       real(rp), pointer :: aux(:,:)

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do ipoin = 1,np
           if (vmass(ipoin) > 0.0_rp) then 
               a%press(ipoin,1)    = a%press(ipoin,1)/vmass(ipoin)
           endif
       enddo

       call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%press(:,1))    

       !Hanging nodes
       if (a%Mesh%kfl_HangingNodes .eqv. .true.) call a%Mesh%InterpolateHangingValues(1,a%press(:,1))
       if (a%Mesh%kfl_perio == 1) call a%Mesh%MasterToSlave(1,a%press(:,1))

   end subroutine AssemblePressEndElmope

   subroutine EndElmopePressGp2ElemArray
       implicit none
       integer(ip)::tn,nd,np,idime,inode

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       do inode = 1,e%pnode

           elpress(1,inode)    = elpress(1,inode)  + gppress(1)*e%shape(inode,e%igaus)*dvol

       end do

   end subroutine EndElmopePressGp2ElemArray

end module
