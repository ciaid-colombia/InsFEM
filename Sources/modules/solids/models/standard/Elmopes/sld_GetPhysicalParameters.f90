subroutine sld_GetPhysicalParameters(a,J,sz,densi,c_elas,ielem)
   use typre
   use Mod_Solids
   implicit none
   class(SolidsProblem), target :: a
   integer(ip) :: ndime,sz
   real(rp)    :: J,densi, lam, mu,young,nu
   real(rp)    :: c_elas(sz,sz), aux, aux_esf, aux_dev
   real(rp), allocatable :: pStrain(:)
   integer(ip), optional :: ielem
   
   call a%Mesh%GetNdime(ndime)

   densi = a%densi

   if(a%sld_model == 'NEOHO') then


       !Elastic tensor for neohookean materials
       lam = a%lambda
       mu  = a%mu -lam*log(J) 

       if (ndime.eq.2) then

           c_elas =  reshape([ lam + 2.0_rp*mu, lam             ,  0.0_rp, &
                                    lam       , lam + 2.0_rp*mu ,  0.0_rp, &
                                   0.0_rp     ,  0.0_rp         , mu    ], &
               [3,3])

       elseif (ndime.eq.3) then

           c_elas =  reshape([ lam + 2.0_rp*mu, lam, lam, 0.0_rp, 0.0_rp , 0.0_rp, &
                               lam,  lam + 2.0_rp*mu, lam, 0.0_rp, 0.0_rp , 0.0_rp, &
                               lam, lam,  lam + 2.0_rp*mu, 0.0_rp, 0.0_rp , 0.0_rp, &
                               0.0_rp,  0.0_rp, 0.0_rp, mu,  0.0_rp, 0.0_rp, &
                               0.0_rp,  0.0_rp, 0.0_rp,  0.0_rp, mu, 0.0_rp, &
                               0.0_rp,  0.0_rp, 0.0_rp,  0.0_rp,  0.0_rp, mu], &
               [6,6])

       endif

       c_elas =  (1.0_rp/J)*c_elas


   else if(a%sld_model == 'STVEN') then

       c_elas = (1.0_rp/J)*a%c_elastic

   else 

       if(a%kfl_modProperties) then

           call a%Memor%alloc(ndime,pStrain,'pStrain','sld_getPhysical')

           pStrain=0.0_rp
           pStrain(:)=a%pStrain(ielem)%a(:)
           nu    = a%poisson
           !From eqn (13), based on square mean value of strains
           !TODO

           !From eqn (15), based on element strain energy density
           !young = ((a%young*(1.0_rp-nu))/(2.0_rp*(1.0_rp-2.0_rp*nu)*(1.0_rp+nu)))*((pStrain(1)**2+pStrain(2)**2) +&
           !    &(2.0_rp*nu/(1.0_rp-nu))*(pStrain(1) + pStrain(2)))
           
           !From eqn (17), based on distortion energy density, default
           if (ndime.eq.2) then
               young = (a%young/(12.0_rp*(1.0_rp+nu)))*((pStrain(1)-pStrain(2))**2 +&
                       &(pStrain(1))**2 + (pStrain(2))**2)
           else
               young = (a%young/(12.0_rp*(1.0_rp+nu)))*((pStrain(1)-pStrain(2))**2 +&
                       &(pStrain(3)-pStrain(1))**2 + (pStrain(2)-pStrain(3))**2)
           endif


           if(young <= zesld) then
               !If the special case of a totally constrained element occurs,
               !lets say all nodes lie on the boundary, the strain would
               !be equal to zero giving a null young coefficient, hence 
               !it's replaced by a small value
               young = 0.1_rp
           end if

          if(a%kfl_printPrincipalStresses) then

              !this is allocated by the ale module
              a%printPStrain(ielem)%a(:,:) = young

          end if

           aux_esf=nu/(1.0_rp-nu)
           aux_dev=(1.0_rp-2.0_rp*nu)/(2.0_rp-2.0_rp*nu)
           aux=young*(1.0_rp-nu)/(1.0_rp-nu-2.0_rp*nu*nu)

           !Elastic tensor for isotropic materials
           if (ndime.eq.2) then

               c_elas =  reshape([ 1.0_rp, aux_esf,  0.0_rp, &
                                  aux_esf,  1.0_rp,  0.0_rp, &
                                   0.0_rp,  0.0_rp, aux_dev], &
                                  [3,3])

           elseif (ndime.eq.3) then

               c_elas =  reshape([ 1.0_rp, aux_esf, aux_esf, 0.0_rp, 0.0_rp , 0.0_rp, &
                                  aux_esf,  1.0_rp, aux_esf, 0.0_rp, 0.0_rp , 0.0_rp, &
                                  aux_esf, aux_esf,  1.0_rp, 0.0_rp, 0.0_rp , 0.0_rp, &
                                  0.0_rp,  0.0_rp, 0.0_rp, aux_dev,  0.0_rp, 0.0_rp, &
                                  0.0_rp,  0.0_rp, 0.0_rp,  0.0_rp, aux_dev, 0.0_rp, &
                                  0.0_rp,  0.0_rp, 0.0_rp,  0.0_rp,  0.0_rp, aux_dev], &
                                  [6,6])

           endif

           c_elas =  aux*c_elas

           call a%Memor%dealloc(ndime,pStrain,'pStrain','sld_getPhysical')

       else

           c_elas = a%c_elastic

       end if
   end if

   end subroutine
