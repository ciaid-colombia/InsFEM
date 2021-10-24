subroutine lmn_ComputeThermPressure(a)
   use MPI
   use typre
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_Element
   use Mod_LowMach

   implicit none
   class(FiniteElement), pointer :: e => NULL()
   type(TimeIntegratorDt1)  :: Integrator
   class(LowMachProblem)    :: a  
   real(rp)                 :: aux1,aux2,acgas,acsph,actco,con1,con2,dsurf,fact,gammagas,LHSDtinv,rhspther(1),ReferenceDtinv,dvol
   integer(ip)              :: nelem,ielem,inode,inodb,igaus,iboun,igaub,nboun,ndime,npoinlocal,nodecount,nsteps
   real(rp), allocatable    :: heatf(:),grate(:),gpvel(:),elvel(:,:),eltem(:),eltem_ini(:),elsou(:),gpsou(:)
   real(rp)                 :: gptem(1),gptem_ini(1),divve,aux1g,aux2g,vnor
   integer                  :: ierr,ibopo,ipoin
   real(rp), pointer :: exnor(:,:) => NULL()
   logical  :: cycleflag

   aux1   = 0.0_rp
   aux2   = 0.0_rp
   aux1g  = 0.0_rp
   aux2g  = 0.0_rp

   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNboun(nboun)
   call a%Mesh%GetNpoinLocal(npoinLocal)      

   !Element Allocation
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lmn_ComputeThermPressure')

   if (a%kfl_confi == 1) then

      call a%Memor%alloc(e%mnode,eltem,'eltem','lmn_ComputeThermPressure')
      call a%Memor%alloc(e%mnode,eltem_ini,'eltem_ini','lmn_ComputeThermPressure') 

      call a%Mesh%ElementSetPointers(e)
      elements : do ielem = 1,nelem   
         call a%Mesh%ElementLoad(ielem,e)  
         call e%gather(1,eltem(:),a%tempe(:,1))
         call e%gather(1,eltem_ini(:),a%itemp(:))

         nodecount = 0
         do inode = 1,e%pnode
            if (e%lnods(inode) <= npoinLocal) then
               nodecount = nodecount+1
            endif
         enddo
         fact = real(nodecount)/real(e%pnode)

         gauss_points : do igaus=1,e%pgaus
            e%igaus = igaus
            !Cartesian derivatives and Jacobian at center of gravity
            call e%elmdcg
            dvol = e%weigp(e%igaus)*e%detjm
            call e%interpg(1,eltem,gptem)
            call e%interpg(1,eltem_ini,gptem_ini)
            aux1 = aux1 + fact*dvol/gptem_ini(1)
            aux2 = aux2 + fact*dvol/gptem(1)
         end do gauss_points
      end do elements

      CALL MPI_ALLREDUCE(aux2,aux2g,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
      CALL MPI_ALLREDUCE(aux1,aux1g,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
      
      a%pther(1) = a%itpre*aux1g/aux2g

      call a%Memor%dealloc(e%mnode,eltem,'eltem','lmn_ComputeThermPressure')
      call a%Memor%dealloc(e%mnode,eltem_ini,'eltem_ini','lmn_ComputeThermPressure') 

   else

      call php_SetTimeIntegrator(a,Integrator,LHSDtinv,nsteps)
      ReferenceDtinv = a%dtinv
      call Integrator%GetRHS(1,a%pther(3),rhspther)

      call a%Memor%alloc(e%mnode,eltem,'eltem','lmn_ComputeThermPressure')
      call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','lmn_ComputeThermPressure')
      call a%Memor%alloc(ndime,gpvel,'gpvel','lmn_ComputeThermPressure')
      call a%Memor%alloc(ndime,grate,'grate','lmn_ComputeThermPressure')
      call a%Memor%alloc(ndime,heatf,'heatf','lmn_ComputeThermPressure')
      call a%Memor%alloc(e%mnode,elsou,'elsou','lmn_ComputeThermPressure')
      call a%Memor%alloc(1,gpsou,'gpsou','lmn_ComputeThermpressure')

      call MPI_ALLREDUCE(MPI_IN_PLACE,a%dvolt,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
      call a%GetPhysicalParameters(acsph=acsph,actco=actco,acgas=acgas)
      gammagas = acsph / (acsph - acgas)
      con1 = a%dvolt/(gammagas-1)
      con2 = gammagas / (gammagas - 1.0_rp)

!      call a%Mesh%ElementSetPointers(e)
!      do ielem = 1,nelem   
!         call a%Mesh%ElementLoad(ielem,e)  
!         call e%gather(ndime,elvel(:,1:e%pnode),a%veloc(:,:,1))
!         call e%elmdcg
!         call e%divergence(elvel,divve)
!
!         nodecount = 0
!         do inode = 1,e%pnode
!            if (e%lnods(inode) <= npoinLocal) then
!               nodecount = nodecount+1
!            endif
!         enddo
!         fact = real(nodecount)/real(e%pnode)
!
!         do igaus=1,e%pgaus
!            e%igaus = igaus
!            call e%elmdcg
!            dvol = e%weigp(e%igaus)*e%detjm
!            aux2 = aux2 + fact*divve*dvol
!         end do
!      end do

      boundaries: do iboun=1,nboun
         !Load Element
         call a%Mesh%BoundaryLoad(iboun,e)
         call e%elmdel

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

         call e%gatherb(1,eltem(1:e%pnodb),a%tempe(:,1))
         call e%gatherb(ndime,elvel(:,1:e%pnodb),a%veloc(:,:,1))

         nodecount = 0
         do inodb = 1,e%pnodb
            if (e%lnodb(inodb) <= npoinLocal) then
               nodecount = nodecount+1
            endif
         enddo
         fact = real(nodecount)/real(e%pnodb)

         !Compute boundary terms
         do igaub=1,e%pgaub
            e%igaub = igaub
            
            !Derivatives at the boundary
            call e%elmderb
            call e%gradientb(1,eltem,grate)
            heatf   = actco*grate
            
         
            !Calculate exterior Normal
            call e%bounor 
            dsurf=e%weigb(e%igaub)*e%eucta
            call e%interpb(ndime,elvel(:,1:e%pnodb),gpvel)
            aux1 = aux1 + fact*dot_product(e%baloc(:,ndime),heatf)*dsurf
            aux2 = aux2 + fact*dot_product(e%baloc(:,ndime),gpvel)*dsurf
         end do 
      end do boundaries

      call MPI_ALLREDUCE(aux2,aux2g,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
      call MPI_ALLREDUCE(aux1,aux1g,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)

      a%pther(1) = (aux1g + con1*a%dtinv*rhspther(1)) / (LHSDtinv*con1 + con2*aux2g)

      call a%Memor%dealloc(e%mnode,eltem,'eltem','lmn_ComputeThermPressure')
      call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','lmn_ComputeThermPressure')
      call a%Memor%dealloc(ndime,gpvel,'gpvel','lmn_ComputeThermPressure')
      call a%Memor%dealloc(ndime,grate,'grate','lmn_ComputeThermPressure')
      call a%Memor%dealloc(ndime,heatf,'heatf','lmn_ComputeThermPressure')
      call a%Memor%dealloc(e%mnode,elsou,'elsou','lmn_ComputeThermPressure')
      call a%Memor%dealloc(1,gpsou,'gpsou','lmn_ComputeThermpressure')
   end if

   !Element Deallocation
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','lmn_ComputeThermPressure')

end subroutine   
