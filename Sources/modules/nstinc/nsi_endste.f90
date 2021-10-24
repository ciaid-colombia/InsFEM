subroutine nsi_endste(a,itask)
   !This routine ends a time step of the incoma%pressible NS equations.
   use typre
   use Mod_TimeIntegrator
   use Mod_php_SetTimeIntegrator
   use Mod_NavierStokes
   use Mod_nsm_ElasticBoundaryDir
   use Mod_nsi_Statistics
   use Mod_NsiExacso
   implicit none
   class(NavierStokesProblem) :: a
   type(TimeIntegratorDt1)    :: Integrator1
   type(NsiExacso)            :: exacso

   integer(ip) :: itask
   integer(ip) :: ifout
   integer(ip) :: ielem, nelem, ndime, npoin, icomp, ipoin, nsteps, npoinlocal,idime,ibopo,iroty
   real(rp)    :: prhs(1), dpdt, LHSDtinv,poinveloc(a%ndofn-1),poindisp(a%ndofn-1)
   real(rp), pointer     :: exnor(:,:) => NULL()
   real(rp), allocatable :: auxsource(:,:)

   !Exact values
   real(rp), allocatable, dimension(:,:)  :: exveg
   real(rp), allocatable, dimension(:)    :: exvel,exprg
   real(rp)                               :: expre
   real(rp), pointer, dimension(:)        :: coord 

   call a%Mesh%GetNpoinLocal(npoinlocal)
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin) 

   !Pre-CrankNicolson
   if (itask == 1) then
      !Also compute the dynamic subgrid scales
      call a%EndElmope('Endste')

   !Post-CrankNicolson
   elseif (itask == 2) then

      !Update a%velocity and a%pressure
      !Higher order components
      if (a%ncomp > 3) then
         do icomp = a%ncomp, 4,-1
            a%veloc(:,:,icomp) = a%veloc(:,:,icomp-1)
            a%press(:,icomp) = a%press(:,icomp-1)
         enddo
      endif
      !Previous time step component
      a%veloc(:,:,3) = a%veloc(:,:,1)  ! Vn = Vn+1
      a%press(:,3) = a%press(:,1)     

      if(a%kfl_exacs==0)then
      
         !Boundary operation at the end
         call a%EndBouope('Endste')
      
         !Update the subgrid a%velocity
         if (a%kfl_tacsg/=0) then
            call a%Mesh%GetNelem(nelem)
            if (a%kfl_tacsg>0) then
            !Crank Nicolson schemes
               if (a%kfl_tsche_1st_current == 'CN   ') then   
                  do ielem=1,nelem
                     a%vesgs(ielem)%a(:,2,:)=2.0_rp*a%vesgs(ielem)%a(:,1,:)-a%vesgs(ielem)%a(:,2,:)
                     if (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
                        a%vesgs2(ielem)%a(:,2,:)=2.0_rp*a%vesgs2(ielem)%a(:,1,:)-a%vesgs2(ielem)%a(:,2,:)
                     endif
                  end do
               else if(a%kfl_tacsg>0) then
                  do ielem=1,nelem
                     a%vesgs(ielem)%a(:,2,:)=a%vesgs(ielem)%a(:,1,:)
                     if (a%kfl_repro == 2 .or. a%kfl_repro == 3) then
                        a%vesgs2(ielem)%a(:,2,:)=a%vesgs2(ielem)%a(:,1,:)
                     endif
                  end do
               end if
            endif
         end if
      
         if (a%kfl_ElasticBoundary == 1) call EndsteElasticBoundaryDir
      
         !For statistics
         call nsi_CalculateStatistics(a)
      
      elseif (a%kfl_exacs/=0 .AND. a%kfl_timei==1)then
         !Allocate exact components
         call a%Memor%alloc(ndime,exvel,'exvel','nsi_endste')  
         call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsi_endste')      
         call a%Memor%alloc(ndime,exprg,'exprg','nsi_endste') 

         do ipoin=1,npoin
            
            call a%Mesh%GetPointCoord(ipoin,coord)
            
            !Exact solution
            call exacso%nsi_ComputeSolution(ndime,coord,a%ctime,a)
            call exacso%nsi_GetPressure(ndime,expre,exprg)
            call exacso%nsi_GetVelocity(ndime,exvel,exveg)
            
            !Store exact solution values, for postprocess
            do idime=1,ndime
               a%AnalyticalVelocity (idime,ipoin,1) = exvel(idime)
            enddo
            a%AnalyticalPressure(ipoin,1) = expre
         enddo
         
         !Deallocate exact component
         call a%Memor%dealloc(ndime,exvel,'exvel','nsi_endste')  
         call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsi_endste')      
         call a%Memor%dealloc(ndime,exprg,'exprg','nsi_endste')
      endif
   endif

   20 format(21(1x,e14.7))

end subroutine nsi_endste
