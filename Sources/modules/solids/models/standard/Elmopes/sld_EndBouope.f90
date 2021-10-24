subroutine sld_EndBouope(sldProblem)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!  Calculates tractions in the boundary to check convergence for the 
!  FSI loops
!
!-----------------------------------------------------------------------
   use typre 
   use Mod_Mesh
   use Mod_Element
   use Mod_memor
   use Mod_iofile
   use Mod_int2str
   use Mod_Solids
   use Mod_sld_BaseElmope
   
   implicit none
   class(SolidsProblem), target :: SldProblem
   type(MemoryMan)       :: Memor
   type(FemMesh)         :: Mesh  
   integer(ip)           :: ndime,idime,iboun,igaub,nboun,ibody
   integer(ip)           :: inode,npoin,inodb,npoinLocal,sz
   integer(ip), pointer  :: lbody => NULL()
   integer               :: ibopo,ipoin
   real(rp), pointer     :: exnor(:,:) => NULL()
   real(rp)              :: dsurf,dsurf0,vnor
   real(rp), allocatable :: tract(:), etraction(:,:),stress_t(:,:)
   real(rp), allocatable :: vmassBoundary(:)
   logical  :: cycleflag

   a=>SldProblem

   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)   
   call a%Mesh%GetNpoinLocal(npoinLocal)     
   call a%Mesh%GetNpoin(npoin)     

   sz  = (ndime*(ndime+1))/2
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_EndBouope')      

   call AllocateSolidBase
   call AllocateBaseElmopeArrays
   if(a%sld_type== 'NONLI' ) then
       call AllocateNonLinearSolidArrays
   endif
   
   !Allocation
   call a%Memor%alloc(ndime,ndime,stress_t,'stress_t','sld_EndBouope')
   call a%Memor%alloc(ndime,tract,'tract','sld_EndBouope')   
   call a%Memor%alloc(ndime,e%mnode,etraction,'etraction','sld_EndBouope')   
   call a%Memor%alloc(e%mnode,elvmass,'elvmass','sld_EndBouope')
   call a%Memor%alloc(npoin,vmassBoundary,'vmassBoundary','sld_EndBouope')

   !Inializations
   a%btraction_nodal=0.0_rp

   ibody=0
   elvmass = 0.0_rp
   vmassBoundary = 0.0_rp
   
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
            call vecnor(exnor(:,1),ndime,vnor,2)
            if (vnor == 0.0_rp) cycleflag =  .true. 
         end if
      end do

      if (cycleflag) cycle
      
      call a%Mesh%GetLbody(lbody,iboun)
       
      ibody=lbody
      

      call e%elmdel
      eldisp=0.0_rp

      call e%gather(ndime,eldisp(:,:,1),a%disp(:,:,1))


      dsurf0 = 0.0_rp           
      etraction=0.0_rp
      elvmass = 0.0_rp

      !Gauss-Point Loop
      do igaub=1,e%pgaub
          e%igaub = igaub

          !Traction at gauss point
          tract  = 0.0_rp
          strain = 0.0_rp
          stress = 0.0_rp
          stress_t=0.0_rp

          !Calculate exterior Normal
          call e%bounor

          dsurf  = e%weigb(e%igaub)*e%eucta
          dsurf0 = dsurf0 + dsurf

          !Derivatives at the boundary
          call e%elmderb         

          !Stress calculation
          call calculateBoundaryGradientsAndDeter
          call calculateStress(ndime)
          call getStressTensor(sz,ndime,stress,stress_t)

          !----------------------------------------------------------------
          !Traction value
          !  t=n*S  S=depends on constitutive model
          do idime =1,ndime
              tract(idime) = dot_product(e%baloc(:,ndime),stress_t(:,idime))
          end do

          !Now we fill traction vector
          do inodb = 1,e%pnodb
             inode = e%lboel(inodb)

              do idime=1,ndime
                  etraction(idime,inode) = etraction(idime,inode) + tract(idime)*e%shapb(inodb,e%igaub)*dsurf
              end do
              elvmass(inode) = elvmass(inode) + e%shapb(inodb,e%igaub)*dsurf

          end do


      end do 

      call a%Mesh%AssemblyToArray(e,1_ip,elvmass,vmassBoundary)
      call a%Mesh%AssemblyToArray(e,e%ndime,etraction,a%btraction_nodal)

  end do boundaries 

   do ipoin = 1,npoin
      if (vmassBoundary(ipoin) > 0.0_rp) a%btraction_nodal(:,ipoin) = a%btraction_nodal(:,ipoin)/vmassBoundary(ipoin)
   enddo
   call a%Mesh%ArrayCommunicator%GhostCommunicate(e%ndime,a%btraction_nodal)    
   !Hanging nodes
   if (a%Mesh%kfl_HangingNodes .eqv. .true.) call a%Mesh%InterpolateHangingValues(e%ndime,a%btraction_nodal)
   if (a%Mesh%kfl_perio == 1) call a%Mesh%MasterToSlave(e%ndime,a%btraction_nodal)

  if(a%sld_type== 'NONLI' ) then
      call DeallocateNonLinearSolidArrays
  endif
  call DeallocateSolidBase
  call DeallocateBaseElmopeArrays

  !Dellocation
  call a%Memor%dealloc(ndime,tract,'tract','sld_EndBouope')   
  call a%Memor%dealloc(ndime,e%mnode,etraction,'etraction','sld_EndBouope')   

  !if FSI now we take the values from the boundaries we need
  call a%Memor%dealloc(e%mnode,elvmass,'elvmass','sld_EndBouope')
  call a%Memor%dealloc(npoin,vmassBoundary,'vmassBoundary','sld_EndBouope')
  call a%Memor%dealloc(ndime,ndime,stress_t,'stress_t','sld_EndBouope')

  !Element Deallocation
  call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','sld_EndBouope')
end subroutine

