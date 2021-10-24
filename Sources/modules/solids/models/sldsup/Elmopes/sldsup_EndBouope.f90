subroutine sldsup_EndBouope(sldProblem)
!-----------------------------------------------------------------------
!  
! DESCRIPTION
!  Calculates tractions in the boundary to check convergence for the 
!  FSI loops
!
!-----------------------------------------------------------------------
   use typre 
   use Mod_SUPSolids
   use Mod_sld_BaseElmope
   
   implicit none
   class(SUPSolidsProblem), target :: sldProblem
   integer(ip)           :: ndime,idime,iboun,igaub,nboun,ibody
   integer(ip)           :: inode,npoin,inodb
   integer(ip), pointer  :: lbody => NULL()
   integer               :: ibopo,ipoin,tn
   real(rp), pointer     :: exnor(:,:) => NULL()
   real(rp)              :: dsurf,dsurf0,vnor
   real(rp), allocatable :: tract(:), etraction(:,:),stress_t(:,:)
   real(rp), allocatable :: sigma_t(:,:),bosigma(:,:),bopre(:)
   real(rp), allocatable :: vmassBoundary(:)
   logical  :: cycleflag

   a  =>sldProblem
   sup=>sldProblem

   !Initializations
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)   
   call a%Mesh%GetNpoin(npoin)     
   tn  = (ndime*(ndime+1))/2

   a%btraction_nodal=0.0_rp
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sldsup_EndBouope')      

   call AllocateBaseElmopeArrays
   
   !Allocation
   call a%Memor%alloc(ndime  ,tract  ,          'tract'    ,'sldsup_EndBouope')   
   call a%Memor%alloc(e%mnode,elvmass,          'elvmass'  ,'sldsup_EndBouope')
   call a%Memor%alloc(ndime  ,ndime  ,sigma_t  ,'sigma_t'  ,'sldsup_EndBouope')
   call a%Memor%alloc(ndime  ,ndime  ,stress_t ,'stress_t' ,'sldsup_EndBouope')
   call a%Memor%alloc(        e%mnodb,bopre    ,'bopre'    ,'sldsup_EndBouope')
   call a%Memor%alloc(tn     ,e%mnodb,bosigma  ,'bosigma'  ,'sldsup_EndBouope')
   call a%Memor%alloc(tn,     gpsigma,'gpsigma',            'sldsup_EndBouope')
   call a%Memor%alloc(ndime  ,e%mnode,etraction,'etraction','sldsup_EndBouope')   
   call a%Memor%alloc(npoin  ,vmassBoundary    ,'vmassBoundary','sldsup_EndBouope')

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

      bosigma=0.0_rp
      bopre  =0.0_rp

      call e%gatherb(tn  ,bosigma(:,:),sup%sigma(:,:,1))
      call e%gatherb(1_ip,bopre,sup%press(:,1))

      dsurf0    = 0.0_rp           
      etraction = 0.0_rp
      elvmass   = 0.0_rp

      !Gauss-Point Loop
      do igaub=1,e%pgaub
          e%igaub = igaub

          !Traction at gauss point
          tract    = 0.0_rp
          stress_t = 0.0_rp
          sigma_t  = 0.0_rp
          gppress  = 0.0_rp

          !Calculate exterior Normal
          call e%bounor

          dsurf  = e%weigb(e%igaub)*e%eucta
          dsurf0 = dsurf0 + dsurf

          !Derivatives at the boundary
          call e%elmderb         

          call e%interpb(tn,bosigma(:,:),gpsigma(:))
          call e%interpb(1_ip,bopre,gppress(1_ip))

          !Deviatoric stress to tensor notation
          call getStressTensor(tn,ndime,gpsigma,sigma_t)

          call get2ndIITensor(ndime,ii)
          stress_t = sigma_t + ii*gppress(1_ip)

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

  call DeallocateBaseElmopeArrays

  !Dellocation
  call a%Memor%dealloc(ndime  ,tract  ,          'tract'    ,'sldsup_EndBouope')   
  call a%Memor%dealloc(e%mnode,elvmass,          'elvmass'  ,'sldsup_EndBouope')
  call a%Memor%dealloc(ndime  ,ndime  ,sigma_t  ,'sigma_t'  ,'sldsup_EndBouope')
  call a%Memor%dealloc(ndime  ,ndime  ,stress_t ,'stress_t' ,'sldsup_EndBouope')
  call a%Memor%dealloc(        e%mnodb,bopre    ,'bopre'    ,'sldsup_EndBouope')
  call a%Memor%dealloc(tn     ,e%mnodb,bosigma  ,'bosigma'  ,'sldsup_EndBouope')
  call a%Memor%dealloc(tn,     gpsigma,'gpsigma',            'sldsup_EndBouope')
  call a%Memor%dealloc(ndime  ,e%mnode,etraction,'etraction','sldsup_EndBouope')   
  call a%Memor%dealloc(npoin  ,vmassBoundary    ,'vmassBoundary','sldsup_EndBouope')

  !Element Deallocation
  call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','sld_EndBouope')

end subroutine
