subroutine sld_turnonGeneral(a)
   ! DESCRIPTION
   !    This routine opens the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_Solids
   implicit none
   class(SolidsProblem)  :: a
   integer(ip)           :: ndime,idime,count,icount,npoin,nod,sz
   real(rp), allocatable :: auxpoin(:,:)
   logical(lg)           :: pass

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)

   !Material constants
   a%mu       = a%young/(2.0_rp + 2.0_rp*a%poisson)
   a%lambda   = a%poisson*a%young/((1.0_rp + a%poisson)*(1.0_rp - 2.0_rp*a%poisson))

   a%bulkMod  = a%young/(3.0_rp*(1.0_rp - 2.0_rp*a%poisson))

   !We modify the size of pwSource so that PETSC can recieve only
   !the nodes that will be modified
   if (allocated(a%pwSource)) then

       call a%Memor%alloc(a%ndofn+1,npoin,auxpoin,'auxpoin','sld_turnon')

       !Find non zero entries
       count=0
       do nod = 1,npoin

           pass =.false.
           if(a%pwSource(1,nod)/=0) then 
               pass=.true.
           endif

           if (pass) then
               auxpoin(:,count+1) = a%pwSource(:,nod)
               count = count +1
           endif

       enddo

       a%pwSourceDim = a%ndofn + 1

       if(count > 0) then

           a%pwSourceSize=count

           call a%Memor%dealloc(a%ndofn+1,npoin,a%pwSource,'pwSource','sld_turnon')
           call a%Memor%dealloc(npoin,a%pwSourceId,'pwSourceId','sld_turnon')

           !Realloc with reduced entries
           call a%Memor%alloc(a%ndofn+1,count,a%pwSource,'pwSource','sld_turnon')
           call a%Memor%alloc(count,a%pwSourceId,'pwSourceId','sld_turnon')

           do icount = 1,count
               do idime= 1,a%ndofn+1
                   a%pwSource(idime,icount)=auxpoin(idime,icount)
               enddo
           enddo
       else
           a%pwSourceSize= npoin

       end if

       call a%Memor%dealloc(a%ndofn+1,npoin,auxpoin,'auxpoin','sld_turnon')

       call a%ModifyPointForces(a%ndofn,a%pwSourceSize,a%pwSource(2:a%pwSourceDim,1:a%pwSourceSize))

   end if

   call a%SolidSpecificTurnon

   if(a%sld_model .eq.'LINEA') then
       a%sld_type = 'LINEA'
   else

       a%sld_type = 'NONLI'
       a%kfl_linea= 2 !For now we default nonlinear sld to NR linearization

       !So we can move the mesh
       call a%Mesh%SetALE(1_ip)
       call a%Mesh%SetDisplacements(a%disp)
   endif

end subroutine sld_turnonGeneral

subroutine sld_turnonSpecific(a)
   ! DESCRIPTION
   !    This routine opens the run for the elastic solid prob
   !-----------------------------------------------------------------------
   use typre
   use Mod_iofile
   use Mod_Solids
   implicit none
   class(SolidsProblem) :: a
   integer(ip)          :: ndime,sz
   
   interface 
      subroutine sld_initElasticTe(a,ndime,sz,young,nu,c)
         use typre
         use Mod_Solids
         implicit none
         class(SolidsProblem) :: a
         integer(ip)          :: ndime,sz
         real(rp)             :: young,nu
         real(rp)             :: c(sz,sz)
      end subroutine
   end interface

   call a%Mesh%GetNdime(ndime)
   sz = (ndime*(ndime+1))/2

   call sld_initElasticTe(a,ndime,sz,a%young,a%poisson,a%c_elastic)

end subroutine sld_turnonSpecific

subroutine sld_initElasticTe(a,ndime,sz,young,nu,c)
   use typre
   use Mod_Solids
   implicit none

   class(SolidsProblem) :: a
   integer(ip)          :: ndime,sz
   real(rp) :: c(sz,sz)
   real(rp) :: young  ,&           ! Young's modulus
                nu     ,&          ! Poisson
                lam    ,&
                mu     ,&
                aux    ,&
                aux_esf,&
                aux_dev

        !Auxiliary variables for esferic and deviatoric
        !components
        aux_esf=nu/(1.0_rp-nu)
        aux_dev=(1.0_rp-2.0_rp*nu)/(2.0_rp-2.0_rp*nu)
        aux=young*(1.0_rp-nu)/(1.0_rp-nu-2.0_rp*nu*nu)

        !Elastic tensor for isotropic materials
        if (ndime.eq.2) then

            c =  reshape([ 1.0_rp, aux_esf,  0.0_rp, &
                          aux_esf,  1.0_rp,  0.0_rp, &
                           0.0_rp,  0.0_rp, aux_dev], &
                          [3,3])

        elseif (ndime.eq.3) then

            c =  reshape([ 1.0_rp, aux_esf, aux_esf,  0.0_rp, 0.0_rp , 0.0_rp, &
                          aux_esf,  1.0_rp, aux_esf,  0.0_rp, 0.0_rp , 0.0_rp, &
                          aux_esf, aux_esf,  1.0_rp,  0.0_rp, 0.0_rp , 0.0_rp, &
                           0.0_rp,  0.0_rp,  0.0_rp, aux_dev,  0.0_rp, 0.0_rp, &
                           0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp, aux_dev, 0.0_rp, &
                           0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp,  0.0_rp, aux_dev], &
                          [6,6])

        endif

        c =  aux*c

end subroutine
