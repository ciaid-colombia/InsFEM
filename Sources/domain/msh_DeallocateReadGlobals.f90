subroutine msh_DeallocateReadGlobals(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   integer(ip) :: auxcount,ielem
   
   
      !Inverse numbering 
      if (allocated(a%igPointGlobNumber)) then
         call a%Memor%dealloc(a%gnpoin,a%igPointGlobNumber,'igPointGlobNumber','ScatterData')
         call a%Memor%dealloc(a%gnpoin,a%gPointGlobNumber,'gPointGlobNumber','Partitionate')
         call a%Memor%dealloc(a%gnpoin,a%gPointProcNumber,'gPointProcNumber','Partitionate')
      endif
      
   if (a%MPIrank == a%MPIroot) then
      if (a%kfl_perio == 1) then
         call a%Memor%dealloc(a%gnpoin,a%gMasterSlave,'gMasterSlave','DeallocateLocals')
      endif   
   endif   
      
   !Parallel reading
   if (a%kfl_MeshInputType == 1) then
         if (a%kfl_ReadType == 1) then
            !ParallelOrdering
            call a%RP_LocalOrdering%Deallocate(a%Memor)
            call a%ParallelLibrary%DeallocateOrdering(a%RP_LocalOrdering,a%Memor)
            call a%RP_InitialOrdering%Deallocate(a%Memor)
            call a%ParallelLibrary%DeallocateOrdering(a%RP_InitialOrdering,a%Memor)
            
            call a%Memor%dealloc(a%MPIsize,a%RP_NpoinOthersNeed,'RP_NpoinOthersNeed','DeallocateReadGlobals')
            call a%Memor%dealloc(a%MPIsize,a%RP_NpoinINeed,'RP_NpoinINeed','DeallocateReadGlobals')
         endif
   endif
   
   if (a%kfl_MeshInputType == 0) then
   !Mesh is read from disk
   !This part is not done for the manufactured mesh
      !Initial data for element reading
      call a%ElementLocal2Initial%Deallocate(a%Memor)
      call a%ElementNaive2Initial%Deallocate(a%Memor)
      call a%ParallelLibrary%DeallocateOrdering(a%ElementLocal2Initial, a%Memor)
      call a%ParallelLibrary%DeallocateOrdering(a%ElementNaive2Initial, a%Memor)
      
      !Serial read
      if (a%kfl_ReadType == 0) then
         if (associated(a%gElementInitialSendToProcessor)) then
            auxcount = 0
            do ielem = 1,size(a%gElementInitialSendToProcessor)
               if (associated(a%gElementInitialSendToProcessor(ielem)%l)) then
                  auxcount = auxcount + size(a%gElementInitialSendToProcessor(ielem)%l)
                  deallocate(a%gElementInitialSendToProcessor(ielem)%l)
               endif
            enddo
            call a%Memor%deallocObj(0,'gElementInitialSendToProcessor','DeallocateReadGlobals',auxcount*ip)
            call a%Memor%deallocObj(0,'gElementInitialSendToProcessor','DeallocateReadGlobals',size(a%gElementInitialSendToProcessor)*ip)
            deallocate(a%gElementInitialSendToProcessor)
         endif
      elseif (a%kfl_ReadType == 1) then
         !Naive Partitioned read
         auxcount = 0
         if (allocated(a%ElementNaiveSendToProcessor)) then
            do ielem = 1,size(a%ElementNaiveSendToProcessor)
               if (associated(a%ElementNaiveSendToProcessor(ielem)%l)) then
                  auxcount = auxcount + size(a%ElementNaiveSendToProcessor(ielem)%l)
                  deallocate(a%ElementNaiveSendToProcessor(ielem)%l)
               endif
            enddo
            call a%Memor%deallocObj(0,'gElementInitialSendToProcessor','Partitionate',auxcount*ip)
            call a%Memor%dealloc(size(a%ElementNaiveSendToProcessor),a%ElementNaiveSendToProcessor,'gElementInitialSendToProcessor','Partitionate')
         endif
      
      endif
   endif   
   
end subroutine