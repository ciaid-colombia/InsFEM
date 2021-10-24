module Mod_nsc_solite_ex
   use typre
   use Mod_Mesh
   use Mod_Memor
   use Mod_NSCompressibleExplicit
   use Mod_NSCompressibleSubroutines

   implicit none   
   real(rp), pointer :: vmass(:) => NULL()

   real(rp)              :: dtime,incre,increeva
   real(rp)              :: ExplicitCoeff(4,2)
   integer(ip)           :: istage,nstage
   integer(ip)           :: ndime,ipoin,npoin,ielem,nelem


contains

   subroutine InitExplicitTimeIntegrator(a)
      implicit none   
      class(NSCompressibleExplicitProblem) :: a
      if (a%kfl_tsche_ex == 'EULER') then
          nstage = 1
         ExplicitCoeff(1,1) =  1.0_rp  
         ExplicitCoeff(1,2) =  1.0_rp
      else if (a%kfl_tsche_ex == 'RK2HE') then
         nstage = 2
         ExplicitCoeff(1,1) =  1.0_rp  
         ExplicitCoeff(1,2) =  0.5_rp
         ExplicitCoeff(2,1) =  1.0_rp  
         ExplicitCoeff(2,2) =  0.5_rp 
      else if (a%kfl_tsche_ex == 'RK223') then
         nstage = 2
         ExplicitCoeff(1,1) =  1.0_rp  
         ExplicitCoeff(1,2) =  1.0_rp/4.0_rp
         ExplicitCoeff(2,1) =  2.0_rp/3.0_rp  
         ExplicitCoeff(2,2) =  3.0_rp/4.0_rp 
      else if (a%kfl_tsche_ex == 'RK4') then
         nstage = 4
         ExplicitCoeff(1,1) =  1.0_rp  
         ExplicitCoeff(1,2) =  1.0_rp/6.0_rp
         ExplicitCoeff(2,1) =  0.5_rp  
         ExplicitCoeff(2,2) =  2.0_rp/6.0_rp 
         ExplicitCoeff(3,1) =  0.5_rp  
         ExplicitCoeff(3,2) =  2.0_rp/6.0_rp  
         ExplicitCoeff(4,1) =  1.0_rp
         ExplicitCoeff(4,2) =  1.0_rp/6.0_rp
      else 
      call runend('Nsc_solite_ex: Other explicit integration methods not ready')
      endif
   end subroutine

end module   


subroutine nsc_solite_ex(a)
   use Mod_nsc_solite_ex
   
   implicit none   
   class(NSCompressibleExplicitProblem) :: a

   interface
      subroutine nsc_elmdir_ex(a)
         import NSCompressibleExplicitProblem
         implicit none
         class(NSCompressibleExplicitProblem) :: a
      end subroutine
      subroutine nsc_elmope_ex(a,tincr,task)
         use typre
         import NSCompressibleExplicitProblem
         implicit none
         class(NSCompressibleExplicitProblem) :: a
         real(rp) :: tincr
         character(6) :: task
      end subroutine
      subroutine nsc_bouope_ex(a)
         use typre
         import NSCompressibleExplicitProblem
         implicit none
         class(NSCompressibleExplicitProblem) :: a
      end subroutine
      subroutine nsc_GhostCommunicate(a)
         use typre
         import NSCompressibleExplicitProblem
         implicit none
         class(NSCompressibleExplicitProblem) :: a
      end subroutine

   end interface

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetVmass(vmass)

   call InitExplicitTimeIntegrator(a)

   !Matrices to zero
   a%femti = 0.0_rp
   if (a%kfl_tacsg == 1) then
      do ielem = 1,nelem
         a%sgsti(ielem)%a = 0.0_rp      
      end do
   end if
   
   a%unkno(1,:) = a%densf(:,3)
   a%unkno(2:ndime+1,:) = a%momen(1:ndime,:,3)
   a%unkno(a%ndofn,:) = a%energ(:,3)
   
   do istage = 1,nstage

      dtime = 1.0 / a%dtinv
      incre = dtime * ExplicitCoeff(istage,1)
  
      a%densf(:,1) = a%densf(:,3) + incre*a%femti(1,:)
      a%momen(1:ndime,:,1) = a%momen(1:ndime,:,3) + incre*a%femti(2:ndime+1,:)
      a%energ(:,1) = a%energ(:,3) + incre*a%femti(a%ndofn,:)

      !Time Dependent SGS stage update
      if (a%kfl_tacsg == 1) then
         do ielem = 1,nelem
            a%cosgs(ielem)%a(1,:) = a%cosgs(ielem)%a(3,:) + incre*a%sgsti(ielem)%a(1,:)       
            a%mosgs(ielem)%a(1:ndime,1,:)=a%mosgs(ielem)%a(1:ndime,3,:)+incre*a%sgsti(ielem)%a(2:ndime+1,:)
            a%ensgs(ielem)%a(1,:) = a%ensgs(ielem)%a(3,:) + incre*a%sgsti(ielem)%a(a%ndofn,:)  
         end do
      end if

      call nsc_elmdir_ex(a)

      call nsc_GhostCommunicate(a)

      a%femti = 0.0_rp
     
      !Residual and gradient projection.
      if (a%kfl_shock == 2 .or. a%kfl_repro == 1) call nsc_elmope_ex(a,incre,'Elmpro')

      call nsc_elmope_ex(a,incre,'Elmope')

      call nsc_bouope_ex(a)
     
      do ipoin = 1,npoin
        a%femti(:,ipoin) = vmass(ipoin)*a%femti(:,ipoin)
      enddo

      increeva = dtime * ExplicitCoeff(istage,2)
   
      a%unkno = a%unkno + increeva*a%femti

      !Time Dependent SGS unknown update
      if (a%kfl_tacsg == 1) then
         do ielem = 1,nelem
            a%cosgs(ielem)%a(2,:) = a%cosgs(ielem)%a(2,:) + increeva*a%sgsti(ielem)%a(1,:)
            a%mosgs(ielem)%a(1:ndime,2,:)=a%mosgs(ielem)%a(1:ndime,2,:)+increeva*a%sgsti(ielem)%a(2:ndime+1,:)
            a%ensgs(ielem)%a(2,:) = a%ensgs(ielem)%a(2,:) + increeva*a%sgsti(ielem)%a(a%ndofn,:)
         end do
      end if

   enddo

   !Subscales time step update
   if (a%kfl_tacsg == 1) call a%DynamicSGSUpdate

end subroutine
