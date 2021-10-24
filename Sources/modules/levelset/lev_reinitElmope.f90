subroutine lev_reinitElmope(LevelProblem,auxlevel)
   use typre
   use Mod_LevelSet
   use Mod_lev_baseElmope
   use Mod_lev_HangingNodes
   use Mod_int2str
   implicit none
   class(LevelSetProblem), target :: LevelProblem
   real(rp) :: auxlevel(:)
   
   real(rp) :: grlev(3)
   integer(ip) :: inode,ndime,idime,ipoin,jnode
   
   integer(ip) :: kfl_HangingNodes, kfl_perio
   integer(ip) :: elemstatus, kpoin,nipoin
   
   real(rp) :: weigp(100),pdsurf
   real(rp), pointer :: xloc(:,:) => NULL()
   
   integer(ip) :: iboun,nboun,igaub,inodb
   real(rp) :: xmuit

   external :: NULLSUB
   
   a=>LevelProblem
   
   !SettingThePointers
   call ResetProcedureComposition
   ProcHook%Initializations => NULLSUB
   ProcHook%Gathers         => NULLSUB
   ProcHook%OnIeltyChange   => NULLSUB
   ProcHook%PreGauss        => NULLSUB
   ProcHook%Interpolates    => NULLSUB
   ProcHook%Elext           => NULLSUB
   ProcHook%InGauss         => NULLSUB
   ProcHook%InGaussElmats   => NULLSUB
   ProcHook%PostGaussElmats => NULLSUB
   ProcHook%Testf           => NULLSUB
   ProcHook%PreDirichlet    => NULLSUB
   ProcHook%Finalizations   => NULLSUB
   ProcHook%PhysicalProp    => NULLSUB
   call SetPointersHangingNodes(0)
   call SetPointersHangingNodes(1)
   call SetPointersHangingNodes(100)
   
      
   !Linear System To Zero
   call a%LinearSystem%ToZero
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_elmope')

   !Matrices Alloc
   call a%Memor%alloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','lev_elmope')
   call a%Memor%alloc(1_ip,e%mnode,elrhs,'elrhs','lev_elmope')
   
   !Other Arrays Alloc
   call a%Memor%alloc(e%mnode,a%ncomp-1,ellev,'ellev','lev_elmope')
   call a%Memor%alloc(a%ncomp-1,gplev,'gplev','lev_elmope')
   
   call ProcHook%Initializations
   
   
   elements : do ielem = 1,nelem
      !Load Element
      call a%Mesh%ElementLoad(ielem,e)   
      
      !Hook 
      call ProcHook%PreGauss
      
      elmat = 0.0_rp
      elrhs = 0.0_rp
      
      call e%gather(1,ellev(:,1),auxlevel)
            
      !Cartesian derivatives and Jacobian at center of gravity
      call e%elmdcg
      dvol = e%detjm
      
      call e%gradient(1_ip,ellev,grlev)
      
      call e%interpc(1_ip,ellev,gplev)
      
      do inode = 1,e%pnode
         elmat(1,inode,1,inode) = elmat(1,inode,1,inode) + dvol
         elrhs(1,inode) = elrhs(1,inode) + dvol*gplev(1)/abs(gplev(1))
      enddo
      
      !System matrices
      do inode = 1,e%pnode
         do jnode = 1,e%pnode
            do idime = 1,ndime
               elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + 0.1*e%cartd(idime,inode)*e%cartd(idime,jnode)*dvol
            enddo
         enddo
         !Source term, 1 in one side, -1 in the other
         elrhs(1,inode) = elrhs(1,inode) + 0.1*e%shacg(inode)*dvol*gplev(1)/abs(gplev(1))
      enddo
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      
      if(elemStatus==0)then
         
         call a%CutMesh%GetLocalIntersectionPoints(ielem,xloc)
         call a%CutMesh%GetSurfaceIntersection(ielem,e,pdsurf)
         call a%CutMesh%GetNinters(ielem,nipoin)
         !Compute shape functions at the intersection points
         do kpoin = 1, nipoin         
            weigp(kpoin) = 1.0_rp/nipoin
         end do      
         !The rutine gives the needed shape functions associated 
         call e%SetParticularGaussPoints(a%Memor,nipoin,xloc,weigp)
         
         do igaus = 1,e%pgaus
            e%igaus = igaus
            do inode = 1,e%pnode
               do jnode = 1,e%pnode
                  elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) + &
                     2*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*e%weigp(e%igaus)*pdsurf/e%detjm*200000
                  
               enddo
            enddo
         enddo
      endif
      
      !PreDirichlet
      call ProcHook%PreDirichlet
      
      !Assembly
      call a%LinearSystem%Assembly(e,elmat,elrhs)
   enddo elements

   !Hook 
   call ProcHook%Finalizations

   !Matrices Alloc  
   call a%Memor%dealloc(1_ip,e%mnode,1_ip,e%mnode,elmat,'elmat','lev_elmope')
   call a%Memor%dealloc(1_ip,e%mnode,elrhs,'elrhs','lev_elmope')
   
   !Other Arrays Alloc
   call a%Memor%dealloc(e%mnode,a%ncomp-1,ellev,'ellev','lev_elmope')
   call a%Memor%dealloc(a%ncomp-1,gplev,'gplev','lev_elmope')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_elmope')
   
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%AssemblyPeriodicBC(a%ndofn,a%LinearSystem,a%Memor)

   !Hanging nodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%AssemblyHangingNodesDiag(a%ndofn,a%LinearSystem,a%Memor)

   call a%LinearSystem%Solve(a%unkno)
   
   !Ghostcommunicate
   call a%Mesh%ArrayCommunicator%GhostCommunicate(a%ndofn,a%unkno)
   
   !HangingNodes
   call a%Mesh%GetHanging(kfl_HangingNodes)
   if (kfl_HangingNodes == 1) call a%Mesh%InterpolateHangingValues(a%ndofn,a%unkno)
   
   !Periodic boundary conditions
   call a%Mesh%GetPerio(kfl_perio)
   if (kfl_perio == 1) call a%Mesh%MasterToSlave(a%ndofn,a%unkno)

   auxlevel = a%unkno(1,:)   
   
   !call a%FilePostpr%postpr(auxlevel,'auxlevel2',a%istep,a%ctime,a%Mesh)
      


end subroutine
