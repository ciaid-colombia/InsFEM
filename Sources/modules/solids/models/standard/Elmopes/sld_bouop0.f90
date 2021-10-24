subroutine sld_bouop0(a)
!-----------------------------------------------------------------------
!****f* SOLIDS/sld_bouope
! NAME 
!    sld_bouop0
! DESCRIPTION
!    Solids boundary elemental operations:
!    1. Compute elemental matrix and RHS 
!    2. Impose Dirichlet boundary conditions
!    3. Assemble them
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_Solids
   use Mod_sld_elmdir
   implicit none
   class(SolidsProblem) :: a
   
   class(FiniteElement) , pointer  :: e => NULL()
   real(rp), allocatable :: elmat(:,:,:,:)
   real(rp),target, allocatable :: elrhs(:,:)
   real(rp), allocatable :: welmat(:,:,:,:)
   real(rp), allocatable :: welrhs(:,:)
   real(rp), pointer    :: exnor(:,:) => NULL()
   real(rp) :: vnor,prod,part_prod,w
   
   integer(ip) :: iboun,ibody,nboun,idime,ndime,sz,ierr
   integer              :: ibopo,ipoin
   integer(ip), pointer :: lbody => NULL()
   integer  icount,npoin,npoinLocal,inodb,inod
   real(rp), pointer:: pointRhs(:) => NULL()
   logical  :: cycleflag
   integer(ip) :: kfl_HangingNodes,u1,uf,s1,sf,p1,bc
   
    interface
       subroutine sld_bouope(a,ndime,iboun,e,sz,bc,elmat,elrhs,ndofn,npoinLocal)
          !This subroutine computes the boundary contribution
          !-----------------------------------------------------------------------
             use typre
             use Mod_Element
             use Mod_Solids
             implicit none
             class(SolidsProblem) :: a
             class(FiniteElement)        :: e
             integer(ip), intent(in)    :: ndofn,npoinLocal,bc
             real(rp),intent(inout)     :: elmat(ndofn,e%mnode,ndofn,e%mnode)              ! Element matrices
             real(rp),intent(inout)     :: elrhs(ndofn,e%mnode)
             integer(ip)                :: iboun,sz,ndime
       end subroutine
       
    end interface
   
   !Initializations
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNpoin(npoin)     
   call a%Mesh%GetNpoinLocal(npoinLocal)     
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
   call a%Mesh%GetHanging(kfl_HangingNodes)
   
   call a%Memor%alloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','sld_bouop0')
   call a%Memor%alloc(a%ndofn,e%mnode,elrhs,'elrhs','sld_bouop0')
   
   
   sz = (ndime*(ndime+1))/2
   call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

   ibody=0

   ! Loop over boundaries
   boundaries: do iboun=1,nboun
       !Load Element
       call a%Mesh%BoundaryLoad(iboun,e)

       cycleflag = .false.
       do inodb = 1,e%pnodb
           ipoin = e%lnodb(inodb)
           call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
           if (ibopo == 0) then
               cycleflag =  .true.
           else
               call vecnor(exnor(:,1),e%ndime,vnor,2)
               if (vnor == 0.0_rp) cycleflag =  .true. 
           end if
       end do

       if (cycleflag) cycle

       call a%Mesh%GetLbody(lbody,iboun)

       ibody=lbody

       !Initialize
       elmat=0.0_rp
       elrhs=0.0_rp

       !to delete repeated entities 
       icount = 0
       do inodb = 1,e%pnodb
          if (e%lnodb(inodb) <= npoinLocal) then
              icount = icount +1
          endif
      enddo
      
      !Compute boundary terms
      if(ibody>0 ) then   

          if (associated(a%trac) .or. a%kfl_fixbo(iboun) == 2) then
              call sld_bouope(a,ndime,iboun,e,sz,bc,elmat,elrhs,a%ndofn,npoinLocal)
          end if
      end if

      call a%ModifyBouopeRHS(a%ndofn,e%mnode,elrhs)

      if (kfl_HangingNodes == 1) call a%Mesh%PrepareHangingMatrices(e,elmat,elrhs)
      
      !Boundary conditions
      call sld_elmdir(a,e,elmat,elrhs)

      ! Assembly
      elrhs = a%force_factor*elrhs

      call a%LinearSystem%Assembly(e,elmat,elrhs)

  end do boundaries

  call a%Memor%dealloc(a%ndofn,e%mnode,a%ndofn,e%mnode,elmat,'elmat','sld_bouop0')
  call a%Memor%dealloc(a%ndofn,e%mnode,elrhs,'elrhs','sld_bouop0')


  call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')

end subroutine sld_bouop0
