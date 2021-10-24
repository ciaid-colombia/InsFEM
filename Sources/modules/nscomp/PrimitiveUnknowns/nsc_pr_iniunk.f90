subroutine nsc_pr_iniunk(a)
!DESCRIPTION
!   This routine sets up the initial conditions for the pressure,
!   velocity and temperature fields.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressiblePrimitive
   use Mod_NSCompressibleSubroutines
   use Mod_NscExacso
   implicit none
   class(NSCompressiblePrimitiveProblem) :: a
   type(NscExacso)  :: exacso   
   
   integer(ip) :: icomp,ipoin,idime,ndime,npoin
   real(rp) :: energ,acvis,actco,accph,accvh
   real(rp) :: expre, extem

   !Exact Values
   real(rp), allocatable   :: exprg(:),exvel(:),exveg(:,:),exteg(:)
   !Nodal coordinates
   real(rp), pointer       :: coord(:) => NULL()   
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)
   if(accph-accvh<zensi) then
      call runend('Nsc_iniunk: Gas constant is negative')
    end if

   do ipoin=1,npoin
      if((a%kfl_fixno(1,ipoin)==1) .or. (a%kfl_fixno(1,ipoin) == 0)) then
         a%press(ipoin,a%ncomp) = a%bvess(1,ipoin,1)
      end if
      if((a%kfl_fixno(ndime+2,ipoin)==1) .or. (a%kfl_fixno(ndime+2,ipoin)==0)) then
         a%tempe(ipoin,a%ncomp) = a%bvess(ndime+2,ipoin,1)
      end if
      a%densf(ipoin,1) = (a%press(ipoin,a%ncomp)+a%relpre)/((accph-accvh)*(a%tempe(ipoin,a%ncomp)+a%reltem))
      do idime=1,ndime
         if((a%kfl_fixno(idime+1,ipoin)==1) .or. (a%kfl_fixno(idime+1,ipoin) == 0)) then
            a%veloc(idime,ipoin,a%ncomp) = a%bvess(idime+1,ipoin,1)
            a%momen(idime,ipoin,1) = a%densf(ipoin,1)*a%bvess(idime+1,ipoin,1)
         end if
      end do
      call nsc_ComputeEnergy(ndime,accvh,(a%tempe(ipoin,a%ncomp)+a%reltem),a%veloc(:,ipoin,a%ncomp),energ)
      a%energ(ipoin,1) = a%densf(ipoin,1)*energ
   end do
  
   if(minval(a%densf(:,1))<(-zensi)) then
      call runend('Nsc_iniunk: Initial density is negative')
   end if

   if(a%kfl_incnd==1) then
      call runend('Nsc_iniunk: Initial conditions not implemented')
   end if

   if(a%kfl_exacs/=0) then

      ! Allocate exact components
      call a%Memor%alloc(ndime,exprg,'exprg','nsc_iniunk')   
      call a%Memor%alloc(ndime,exvel,'exvel','nsc_iniunk')   
      call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsc_iniunk')     
      call a%Memor%alloc(ndime,exteg,'exteg','nsc_iniunk')   
         do ipoin = 1,npoin 
            call a%Mesh%GetPointCoord(ipoin,coord)
            
            call exacso%nsc_ComputeSolution(ndime,coord,a)
            call exacso%nsc_GetPressure(ndime,expre,exprg)           
            call exacso%nsc_GetVelocity(ndime,exvel,exveg)           
            call exacso%nsc_GetTemperature(ndime,extem,exteg)           

            a%press(ipoin,a%ncomp) = expre
            do idime=1,ndime
               a%veloc(idime,ipoin,a%ncomp) = exvel(idime)
            end do
            a%tempe(ipoin,a%ncomp) = extem
         end do
      ! Allocate exact components
      call a%Memor%dealloc(ndime,exprg,'exprg','nsc_iniunk')   
      call a%Memor%dealloc(ndime,exvel,'exvel','nsc_iniunk')   
      call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsc_iniunk')     
      call a%Memor%dealloc(ndime,exteg,'exteg','nsc_iniunk')   

   end if

   !Assign var(n,i,*) <-- var(n-1,*,*), initial guess after initialization (or reading restart)

   do icomp = 1,a%ncomp-1
      a%press(1:npoin,icomp) = a%press(1:npoin,a%ncomp)
      a%veloc(1:ndime,1:npoin,icomp) = a%veloc(1:ndime,1:npoin,a%ncomp) 
      a%tempe(1:npoin,icomp) = a%tempe(1:npoin,a%ncomp) 
   enddo
   
   call a%Ifconf

end subroutine nsc_pr_iniunk

