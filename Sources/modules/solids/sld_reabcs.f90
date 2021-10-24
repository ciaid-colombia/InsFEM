subroutine  sld_Reabcs(a,itask,kflag)
   ! NAME
   !    sld_reabcs
   ! DESCRIPTION
   !    This routine reads the specific boundary conditions for SOLID.
   !
   !    For conditions on nodes, bvess(idime,ipoin,1) 
   !    The different codes for kfl_fixno_nsi(ipoin) are:
   !    = 1 ... Dirichlet
   !    = 0 ... Free or initial
   !
   !    For conditions on boundaries, bvnat(iboun) is allocated ONLY if
   !    a condition is imposed on it. Overmore, its length depends on the
   !    condition type. The different codes for kfl_fixbo(iboun) are:
   !    = 1 ... Dirichlet ............ u
   !    At the end of the subroutine conditions on boundaries are
   !    transfered to conditions on nodes.
   !------------------------------------------------------------------------
   use MPI
   use typre
   use Mod_Solids
   implicit none
  
   class(SolidsProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag
   integer(ip) :: iboun,ipoin
   integer(ip) :: nboun,npoin,ndime
   character(150) :: outstr
   integer :: ierr
   
   !Initializations
   if (itask == 0) then
      
      !Output
      outstr = adjustl(trim(a%exmod))//'_SPECIFIC_REABCS'
      
      !Get Mesh data
      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNboun(nboun)
      call a%Mesh%GetNdime(ndime)

      !Arrays for conditions on boundaries
      call a%Memor%alloc(npoin,a%kfl_fixrs,'kfl_fixrs',outstr)
      call a%Memor%alloc(nboun,a%kfl_bours,'kfl_bours', outstr)

   !Header   
   elseif(itask == 1) then
      
    
    
    !Finalization
   elseif(itask == 100) then
      
      !-------------------------------------------------------------------------
      !COMMUNICATE NODAL GHOST DATA WITH NEIGHBOURS
      call a%Mesh%ArrayCommunicator%GhostCommunicate(1,a%kfl_fixrs)
      
      !Final operations
      call a%Mesh%Getnboun(nboun)
      call a%Mesh%Getnpoin(npoin)
      
     
      !Check if there is a local prescription
      nodes: do ipoin=1,npoin
         if(a%kfl_fixrs(ipoin)/=0) then
            a%kfl_local=1
            exit nodes
         end if
      end do nodes

   endif

end subroutine

subroutine sld_ReadOnNodes(a)
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)::ndime,idime

   call a%Mesh%GetNdime(ndime)
      
end subroutine

subroutine sld_ReadSource(a)
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)::ndime,idime,npoin
   integer(ip) :: u1,uf,s1,sf,p1,bc
   logical,save :: doSet=.true.

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

   if(doSet) then
       doSet = .false.
       call a%Memor%alloc(a%ndofn+1,npoin,a%pwSource,'pwSource','sld_reabcs')
       call a%Memor%alloc(npoin,a%pwSourceId,'pwSourceId','sld_reabcs')
   endif


   a%kfl_constPointForce = a%Listener%param(2)
   a%kfl_PointForceTime = a%Listener%param(3)
   a%kfl_sourcesPresent=.true.

   a%pwSource(1,a%gipoin) = a%Listener%param(1)
   do idime=1,ndime
       a%pwSource(bc+idime+1,a%gipoin) = a%Listener%param(3+idime)
   enddo

end subroutine

subroutine sld_ReadOnBoundaries(a)
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)                :: nbody,iboun,nboun,pnodb,ipsta
   
   !Dirichlet
   if(a%kfl_fixbo(a%giboun) == 1) then
      if(a%kfl_conbc==1) then
         a%kfl_bours(a%giboun) = int(a%Listener%param(a%gipsta+1))
      else
         a%kfl_bours(a%giboun) = int(a%Listener%param(a%gipsta+2))
      endif
      
   !Prescribed traction: sig.n = tract_0 (tract_0 given)
   else if(a%kfl_fixbo(a%giboun) == 2) then
      allocate(a%bvnat(a%giboun)%a(3))
      a%bvnat(a%giboun)%a(1:3) = a%Listener%param(a%gipsta+1:a%gipsta+3)
      a%kfl_imposedTraction=1_ip
   endif

end subroutine


