subroutine nsc_iniunk(a)
!DESCRIPTION
!   This routine sets up the initial conditions for the density,
!   momentum and energy fields.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_NSCompressibleSubroutines
   use Mod_NscExacso
   implicit none
   class(NSCompressibleProblem) :: a
   type(NscExacso)  :: exacso   
   
   integer(ip) :: icomp,ipoin,idime,ndime,npoin
   real(rp) :: energ,acvis,actco,accph,accvh
   real(rp) :: exden, exene

   !Exact Values
   real(rp), allocatable   :: exdeg(:),exmom(:),exmog(:,:),exeng(:)
   !Nodal coordinates
   real(rp), pointer       :: coord(:) => NULL()   
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)
   if(accph-accvh<zensi) then
      call runend('Nsc_iniunk: Gas constant is negative')
    end if

   if(minval(a%bvess(1,:,1))<(-zensi)) then
      call runend('Nsc_iniunk: Initial density is negative')
   end if

   if(minval(a%bvess(ndime+2,:,1))<(-zensi)) then
      call runend('Nsc_iniunk: Initial absolute temperature is negative')
   end if

   do ipoin=1,npoin
      if((a%kfl_fixno(1,ipoin)==1) .or. (a%kfl_fixno(1,ipoin) == 0)) then
         a%densf(ipoin,a%ncomp) = a%bvess(1,ipoin,1)
      end if
      if((a%kfl_fixno(ndime+2,ipoin)==1) .or. (a%kfl_fixno(ndime+2,ipoin)==0)) then
         a%tempe(ipoin,1) = a%bvess(ndime+2,ipoin,1)
      end if
      do idime=1,ndime
         if((a%kfl_fixno(idime+1,ipoin)==1) .or. (a%kfl_fixno(idime+1,ipoin) == 0)) then
            a%veloc(idime,ipoin,1) = a%bvess(idime+1,ipoin,1)
            if((a%kfl_fixno(1,ipoin)==1) .or. (a%kfl_fixno(1,ipoin) == 0)) then
               a%momen(idime,ipoin,a%ncomp) = a%densf(ipoin,a%ncomp)*a%bvess(idime+1,ipoin,1)
               if((a%kfl_fixno(ndime+2,ipoin)==1) .or. (a%kfl_fixno(ndime+2,ipoin)==0)) then
                 call nsc_ComputeEnergy(ndime,accvh,a%tempe(ipoin,1),&
                           a%veloc(:,ipoin,1),energ)
                 a%energ(ipoin,a%ncomp) = a%densf(ipoin,a%ncomp)*energ
               end if
            end if
         end if
      end do
   end do
  

   if(a%kfl_incnd==1) then
      call runend('Nsc_iniunk: Initial conditions not implemented')
   end if

   if(a%kfl_exacs/=0) then

      ! Allocate exact components
      call a%Memor%alloc(ndime,exdeg,'exdeg','nsc_iniunk')   
      call a%Memor%alloc(ndime,exmom,'exmom','nsc_iniunk')   
      call a%Memor%alloc(ndime,ndime,exmog,'exmog','nsc_iniunk')     
      call a%Memor%alloc(ndime,exeng,'exeng','nsc_iniunk')   
         do ipoin = 1,npoin 
            call a%Mesh%GetPointCoord(ipoin,coord)
            
            call exacso%nsc_ComputeSolution(ndime,coord,a)
            call exacso%nsc_GetDensity(ndime,exden,exdeg)           
            call exacso%nsc_GetMomentum(ndime,exmom,exmog)           
            call exacso%nsc_GetEnergy(ndime,exene,exeng)           

            a%densf(ipoin,a%ncomp) = exden
            do idime=1,ndime
               a%momen(idime,ipoin,a%ncomp) = exmom(idime)
            end do
            a%energ(ipoin,a%ncomp) = exene
         end do
      ! Allocate exact components
      call a%Memor%dealloc(ndime,exdeg,'exdeg','nsc_iniunk')   
      call a%Memor%dealloc(ndime,exmom,'exmom','nsc_iniunk')   
      call a%Memor%dealloc(ndime,ndime,exmog,'exmog','nsc_iniunk')     
      call a%Memor%dealloc(ndime,exeng,'exeng','nsc_iniunk')   

   end if
   !Assign var(n,i,*) <-- var(n-1,*,*), initial guess after initialization (or reading restart)

   do icomp = 1,a%ncomp-1
      a%densf(1:npoin,icomp) = a%densf(1:npoin,a%ncomp)
      a%momen(1:ndime,1:npoin,icomp) = a%momen(1:ndime,1:npoin,a%ncomp) 
      a%energ(1:npoin,icomp) = a%energ(1:npoin,a%ncomp) 
   enddo
   
end subroutine nsc_iniunk

