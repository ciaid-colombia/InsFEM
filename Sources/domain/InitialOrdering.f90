subroutine InitialOrdering(a)
   use MPI
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   integer(ip), allocatable :: NpoinOthersLocal(:),NpoinOthersDisp(:)
   integer(ip), allocatable :: L2Groot(:), G2Iroot(:)
   integer(ip) :: L2G(a%npoinLocal)
   integer(ip) :: irank,ipoin,ierr,lpoin
   integer(ip) :: LocalArray(a%RP_npoinLocal),LocalArray2(a%RP_npoinGhost)
   integer(ip) :: ProcNumber(a%RP_npoinLocal),ProcNumber2(a%RP_npoinGhost)
   integer(ip) :: GlobalArray(a%RP_npoinLocal),GlobalArray2(a%RP_npoinGhost)
   

   if (a%kfl_MeshInputType == 0) then
      if (a%kfl_ReadType == 0) then
         if (a%npoinLocal > 0) then
            call a%Memor%alloc(a%RP_npoinLocal,a%lPointGlobNumber,'lPointGlobNumber','InitialOrdering')
            call a%Memor%alloc(a%npoinLocal,a%ilPointGlobNumber,'ilPointGlobNumber','InitialOrdering')
            call a%Memor%alloc(a%RP_npoinLocal,a%lPointProcNumber,'lPointProcNumber','InitialOrdering') 
           
            do ipoin = a%RP_poin0+1,a%RP_npoinLocal+a%RP_poin0
               call a%Initial2Local(ipoin,lpoin)
               a%lPointProcNumber(lpoin) = a%gPointProcNumber(ipoin)
               a%lPointGlobNumber(lpoin) = a%gPointGlobNumber(ipoin)
            end do
         end if
         
         !Inverse numbering
            call a%Memor%alloc(a%MPIsize,NpoinOthersLocal,'NpoinOthersLocal','InitialOrdering')
            call a%Memor%alloc(a%MPIsize,NpoinOthersDisp,'NpoinOthersDisp','InitialOrdering')
         if (a%MPIrank == a%MPIroot) then 
            call a%Memor%alloc(a%gnpoin,L2Groot,'L2Groot','InitialOrdering')
            call a%Memor%alloc(a%gnpoin,G2Iroot,'G2Iroot','InitialOrdering')
         else
            call a%Memor%alloc(1,L2Groot,'L2Groot','InitialOrdering')
            call a%Memor%alloc(1,G2Iroot,'G2Iroot','InitialOrdering')
         end if
         call MPI_Gather(a%npoinLocal,1,MPI_INTEGER4,NpoinOthersLocal,1,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
         do ipoin = 1,a%npoinLocal
            L2G(ipoin) = ipoin
         enddo
         call a%LocalOrdering%Local2Global(a%npoinLocal,L2G,L2G)
   
         if (a%MPIrank == a%MPIroot) then 
            NpoinOthersDisp = 0
            do irank = 0,a%MPIsize-2
               NpoinOthersDisp(irank+2) = NpoinOthersLocal(irank+1) + NpoinOthersDisp(irank+1)
            end do
         end if
         call MPI_GatherV(L2G,a%npoinLocal,MPI_INTEGER4,L2Groot,NpoinOthersLocal,NpoinOthersDisp,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
         if (a%MPIrank == a%MPIroot) then 
            do ipoin = 1,a%gnpoin
               G2Iroot(ipoin) = a%igPointGlobNumber(L2Groot(ipoin))
            end do
         end if
         call MPI_ScatterV(G2Iroot,NpoinOthersLocal,NpoinOthersDisp,MPI_INTEGER4,a%ilPointGlobNumber,a%npoinLocal,MPI_INTEGER4,a%MPIroot,a%MPIcomm,ierr)
            call a%Memor%dealloc(a%MPIsize,NpoinOthersLocal,'NpoinOthersLocal','InitialOrdering')
            call a%Memor%dealloc(a%MPIsize,NpoinOthersDisp,'NpoinOthersDisp','InitialOrdering')
         if (a%MPIrank == a%MPIroot) then
            call a%Memor%dealloc(a%gnpoin,L2Groot,'L2Groot','InitialOrdering')
            call a%Memor%dealloc(a%gnpoin,G2Iroot,'G2Iroot','InitialOrdering')
         else
            call a%Memor%dealloc(1,L2Groot,'L2Groot','InitialOrdering')
            call a%Memor%dealloc(1,G2Iroot,'G2Iroot','InitialOrdering')
         end if
    
      endif   
   endif   
      
end subroutine

subroutine DeallocateInitial(a)
   use typre
   use Mod_Mesh
   implicit none
   class(FemMesh) :: a
   
   if (a%kfl_MeshInputType == 0) then
      if (a%kfl_ReadType == 0) then
         if (allocated(a%ilPointGlobNumber)) then
              call a%Memor%dealloc(a%RP_npoinLocal,a%lPointGlobNumber,'lPointGlobNumber','InitialOrdering')
              call a%Memor%dealloc(size(a%ilPointGlobNumber),a%ilPointGlobNumber,'ilPointGlobNumber','InitialOrdering')
              call a%Memor%dealloc(a%RP_npoinLocal,a%lPointProcNumber,'lPointProcNumber','InitialOrdering')
         end if
      end if
      if (a%kfl_ReadType == 1) then
          call a%RP_LocalOrdering%Deallocate(a%Memor)
          call a%ParallelLibrary%DeallocateOrdering(a%RP_LocalOrdering,a%Memor)
          call a%RP_InitialOrdering%Deallocate(a%Memor)
          call a%ParallelLibrary%DeallocateOrdering(a%RP_InitialOrdering,a%Memor)
          
          call a%Memor%dealloc(a%MPIsize,a%RP_NpoinOthersNeed,'RP_NpoinOthersNeed','DeallocateReadGlobals')
          call a%Memor%dealloc(a%MPIsize,a%RP_NpoinINeed,'RP_NpoinINeed','DeallocateReadGlobals')
      endif
   end if
   
end subroutine
