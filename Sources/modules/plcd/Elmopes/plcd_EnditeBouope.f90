subroutine plcd_Enditebouope(a)
!-----------------------------------------------------------------------
!****f* PLCD/plcd_Enditebouope
! NAME 
!    plcd_Enditebouope
! DESCRIPTION
!    PLCD boundary elemental operations at the end of the iteration:
!    1. Compute elemental tractions
!    2. Assemble them
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_PLCD
   use Mod_Element
   implicit none
   class(PLCDProblem) :: a
   class(FiniteElement) , pointer     :: e => NULL()
   
   integer(ip) :: iboun, nboun, igaub, inodb, inode, ndime, ipoin, ibopo, icount, npoinLocal
   real(rp), allocatable   :: traction(:), elTractionForces(:,:), fluidtraction(:,:)
   real(rp) :: dsurf, tractionnorm, vnor
   integer(ip), pointer :: lbody
   logical  :: cycleflag
   real(rp), pointer    :: exnor(:,:)
   
   if (any(a%Stages(a%CurrentStage)%kfl_fixbo(:) == 2) .or. associated(a%fluidtraction)) then
      
      !Initializations
      call a%Mesh%GetNboun(nboun)
      call a%Mesh%GetNdime(ndime)
      call a%Mesh%GetNpoinLocal(npoinLocal) 
   
      !Memory allocation
      
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','plcd_Enditebouope')
      call a%Memor%alloc(ndime,e%mnode,elTractionForces,'elTractionForces','plcd_Enditebouope')
      call a%Memor%alloc(ndime,traction,'traction','plcd_Enditebouope')
      call a%Memor%alloc(ndime,e%pnodb,fluidtraction,'fluidtraction','plcd_Enditebouope')
   
      boundaries: do iboun=1,nboun
      
         call a%Mesh%GetLbody(lbody,iboun)
      
         !Tractions
         if (a%Stages(a%CurrentStage)%kfl_fixbo(iboun) == 2 .or. associated(a%fluidtraction)) then
      
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
            
            !To delete repeated entities 
            icount = 0
            do inodb = 1,e%pnodb
               if (e%lnodb(inodb) <= npoinLocal) then
                  icount = icount +1
               endif
            enddo
            
            !Tractions
            traction=0.0_rp
            elTractionForces=0.0_rp
            fluidtraction=0.0_rp
            
            if (a%Stages(a%CurrentStage)%kfl_fixbo(iboun) == 2) then 
               traction(1:e%ndime) = a%Stages(a%CurrentStage)%bvnat(iboun)%a(1:e%ndime)
            elseif (lbody .GT. 0) then
               call e%gatherb(ndime,fluidtraction,a%fluidtraction(:,:))
            end if
            
            dsurf = 0.0_rp
               
            !Gauss-Point Loop
            do igaub=1,e%pgaub
               e%igaub = igaub
         
               !Derivatives at the boundary
               call e%elmderb
   
               !Calculate exterior Normal
               call e%bounor
               
!                if (a%kfl_LargeStrains == 1) then
!                   call vecnor(traction,e%ndime,tractionnorm,2)
!                   traction(1:e%ndime) = tractionnorm*(-e%baloc(:,e%ndime))
!                endif
         
               dsurf=e%weigb(e%igaub)*e%eucta

               do inodb=1,e%pnodb
                  inode = e%lboel(inodb)
                     elTractionForces(:,inode) = elTractionForces(:,inode) + e%shapb(inodb,e%igaub)*traction(:)*dsurf - e%shapb(inodb,e%igaub)*fluidtraction(:,inodb)*dsurf
               end do
            end do
         
         call a%Mesh%AssemblyToArray(e,a%ndofn,elTractionForces,a%ExternalForcesVector)
         endif
      end do boundaries
   
      call a%Memor%dealloc(ndime,e%mnode,elTractionForces,'elTractionForces','plcd_Enditebouope')
      call a%Memor%dealloc(ndime,traction,'traction','plcd_Enditebouope')
      call a%Memor%dealloc(ndime,e%pnodb,fluidtraction,'fluidtraction','plcd_Enditebouope')
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','plcd_Enditebouope')
      
   endif
   
end subroutine plcd_Enditebouope
