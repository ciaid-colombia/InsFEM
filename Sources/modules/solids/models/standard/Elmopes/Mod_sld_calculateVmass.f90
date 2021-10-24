module Mod_sld_calculateVmass
   use typre
   use Mod_sld_BaseElmope
   implicit none
   private
   public SetPointersCalculateVmass
   integer(ip), allocatable :: kfl_IsSet
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersCalculateVmass(itask)
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
              if(a%kfl_printNodalSigma .or. a%kfl_printNodalPress) then
                  call ConcatenateProcedures(ProcHook%ResetArrays,ResetVmass)

                  call ConcatenateProcedures(ProcHook%AllocateArrays  ,AllocateVmass)
                  call ConcatenateProcedures(ProcHook%DeallocateArrays,DeallocateVmass)

                  call ConcatenateProcedures(ProcHook%PostInterpolates,EndElmopeCalculateVmass)
                  call ConcatenateProcedures(ProcHook%ToLinearSystem,EndElmopeAssemblyToArraysVmass)
              endif
          endif  
      case(100)
          deallocate(kfl_IsSet)
          call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine   
   
   !For actual computations
   !------------------------------------------------------
   subroutine EndElmopeAssemblyToArraysVmass
       implicit none
       integer(ip)::tn,nd

       call a%Mesh%GetNdime(nd)
       tn  = (nd*(nd+1))/2

       call a%Mesh%AssemblyToArray(e,1_ip,elvmass,vmass)

   end subroutine EndElmopeAssemblyToArraysVmass

   subroutine EndElmopeCalculateVmass
       implicit none
       integer(ip)::inode

       do inode = 1,e%pnode

            elvmass(inode) = elvmass(inode) + e%shape(inode,e%igaus)*dvol

       end do

   end subroutine EndElmopeCalculateVmass

   subroutine ResetVmass
       implicit none

       elvmass   = 0.0_rp
   end subroutine ResetVmass

   subroutine AllocateVmass
       implicit none
       integer(ip)::tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       !Allocation
       call a%Memor%alloc(np     ,vmass            ,'vmass'    ,'sldup_EndElmope')
       call a%Memor%alloc(e%mnode,elvmass          ,'elvmass'  ,'sldup_EndElmope')

   end subroutine
   
   subroutine DeallocateVmass
       implicit none
       integer(ip)::tn,nd,np

       call a%Mesh%GetNdime(nd)
       call a%Mesh%GetNpoin(np)     
       tn  = (nd*(nd+1))/2

       !Allocation
       call a%Memor%dealloc(np     ,vmass            ,'vmass'    ,'sldup_EndElmope')
       call a%Memor%dealloc(e%mnode,elvmass          ,'elvmass'  ,'sldup_EndElmope')

   end subroutine
   
end module
