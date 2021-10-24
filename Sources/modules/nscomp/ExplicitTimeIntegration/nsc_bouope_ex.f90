module Mod_nsc_bouope_ex
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_NSCompressibleExplicit
   use Mod_Element

   implicit none
   class(NSCompressibleExplicitProblem), pointer :: a => NULL()

   class(FiniteElement), pointer :: e => NULL()
   real(rp), allocatable :: boden(:),gpden(:)
   real(rp), allocatable :: bomom(:,:),gpmom(:)
   real(rp), allocatable :: botem(:),gptem(:)
   real(rp), allocatable :: tract(:)

   integer(ip) :: iboun,nboun,nelem,ndime
   integer(ip) :: igaub,inodb,inode,idofn
   real(rp)    :: acvis,actco,accph,accvh
   real(rp)    :: dsurf,gptno

end module

subroutine nsc_bouope_ex(NSCompExplicitProblem)
   
!-----------------------------------------------------------------------
!****f* Nscomp/nsc_bouope_ex
! NAME 
!    nsc_bouope_ex
! DESCRIPTION
!    Compressible Navier-Stokes boundary elemental operations:
!    1. Compute elemental RHS 
!    1. Assemble them
!***
!-----------------------------------------------------------------------
   use Mod_nsc_bouwal
   use Mod_nsc_bouope_ex
   implicit none
   class(NSCompressibleExplicitProblem), target :: NSCompExplicitProblem
   a=>NSCompExplicitProblem
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   
   
   !Memory allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','allocs')
  
   !Variables!
   call a%Memor%alloc(e%pnodb,boden,'boden','nsc_bouope_ex')
   call a%Memor%alloc(1,gpden,'gpden','nsc_bouope_ex')
   call a%Memor%alloc(e%ndime,e%pnodb,bomom,'bomom','nsc_bouope_ex')
   call a%Memor%alloc(e%ndime,gpmom,'gpmom','nsc_bouope_ex')
   call a%Memor%alloc(a%ndofn,tract,'tract','nsc_bouope_ex')

   !Primitive! 
   call a%Memor%alloc(e%pnodb,botem,'botem','nsc_bouope_ex')
   call a%Memor%alloc(1,gptem,'gptem','nsc_bouope_ex')
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun

      !Load Element
      call a%Mesh%BoundaryLoad(iboun,e)

      if(    a%kfl_fixbo(iboun)==3.or.&   !Mom=Neu  Ene=Free
         & a%kfl_fixbo(iboun)==4.or.&     !Mom=Neu  Ene=Neu
         & a%kfl_fixbo(iboun)==5.or.&     !Mom=Neu  Ene=Rob
         & a%kfl_fixbo(iboun)==6.or.&     !Mom=Free  Ene=Neu
         & a%kfl_fixbo(iboun)==7.or.&     !Mom=Free  Ene=Rob
         & a%kfl_fixbo(iboun)==8.or.&     !Mom=Wall_law Ene=Free
         & a%kfl_fixbo(iboun)==9 &        !Mom= backlow penalty
         ) then  

         !Initialize
         tract=0.0_rp
         
         !Physical Parameters
         call a%GetPhysicalParameters(acvis,actco,accph,accvh)
         
         call e%elmdel
         
         !Gather operations
         call e%gatherb(1_ip,boden,a%densf(:,1))
         call e%gatherb(e%ndime,bomom,a%momen(:,:,1))
         call e%gatherb(1_ip,botem,a%tempe(:,1))
         
         dsurf = 0.0_rp
         
         !Gauss-Point Loop
         gauss_points : do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Initialize
            tract=0.0_rp
      
            !Calculate exterior Normal
            call e%bounor
            
            dsurf=e%weigb(e%igaub)*e%eucta
            
            !Derivatives at the boundary
            call e%elmderb
   
            !Mom=Neu Ene=Free
     	    if(a%kfl_fixbo(iboun)==3) then
               tract(2:e%ndime+1)=-a%bvnat(iboun)%a(1:e%ndime)*e%baloc(1:e%ndime,e%ndime)
            
            !Mom=Neu Ene=Neu
            else if(a%kfl_fixbo(iboun)==4) then
               tract(2:e%ndime+1)=-a%bvnat(iboun)%a(1:e%ndime)*e%baloc(1:e%ndime,e%ndime)
               tract(e%ndime+2)=-a%bvnat(iboun)%a(4)

            !Mom=Neu Ene=Rob
            else if(a%kfl_fixbo(iboun)==5) then
               tract(2:e%ndime+1)=-a%bvnat(iboun)%a(1:e%ndime)*e%baloc(1:e%ndime,e%ndime)
   	       !Temperature value
               call e%interpb(1_ip,botem,gptem(1))
               call vecnor(gptem,1_ip,gptno,2)
               tract(e%ndime+2)=-a%bvnat(iboun)%a(6)*(a%bvnat(iboun)%a(5)-gptno)

            !Mom=Fre Ene=Neu
            else if(a%kfl_fixbo(iboun)==6) then
               tract(e%ndime+2)=-a%bvnat(iboun)%a(1)

            !Mom=Fre Ene=Rob
            else if(a%kfl_fixbo(iboun)==7) then
   	       !Temperature value
               call e%interpb(1_ip,botem,gptem(1))
               call vecnor(gptem,1_ip,gptno,2)
               tract(e%ndime+2)=-a%bvnat(iboun)%a(3)*(a%bvnat(iboun)%a(2)-gptno)

            !Mom=Wlaw Ene=Free
            !Wall law: sig.n=-rho*(U*^2)*u/|u|
            else if(a%kfl_fixbo(iboun)==8) then
               call e%interpb(1_ip,boden,gpden(1))
               call e%interpb(e%ndime,bomom,gpmom)
               call nsc_bouwal(e,gpmom/gpden(1),acvis,gpden(1),a%delta,tract)

            !Backflow penalty
            else if(a%kfl_fixbo(iboun)==9) then
               call e%interpb(e%ndime,bomom,gpmom)
               call nsc_bouback_ex(e,gpmom,a%bvnat(iboun)%a(1),tract)

            end if

            !Traction to vrhs
            do inodb=1,e%pnodb
               inode = e%lnodb(inodb)
               do idofn=1,a%ndofn
                     a%femti(idofn,inode)=a%femti(idofn,inode)+dsurf*tract(idofn)*e%shapb(inodb,e%igaub)
               end do
            end do
   
         end do gauss_points
	
!  call nsc_elmdir_ex(a)

      end if 

   end do boundaries

   call a%Memor%dealloc(e%pnodb,boden,'boden','nsc_bouope_ex')
   call a%Memor%dealloc(1,gpden,'gpden','nsc_bouope_ex')
   call a%Memor%dealloc(e%ndime,e%pnodb,bomom,'bomom','nsc_bouope_ex')
   call a%Memor%dealloc(e%ndime,gpmom,'gpmom','nsc_bouope_ex')
   call a%Memor%dealloc(e%pnodb,botem,'botem','nsc_bouope_ex')
   call a%Memor%dealloc(1,gptem,'gptem','nsc_bouope_ex')
   call a%Memor%dealloc(a%ndofn,tract,'tract','nsc_bouope_ex')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','allocs')
   
end subroutine nsc_bouope_ex
