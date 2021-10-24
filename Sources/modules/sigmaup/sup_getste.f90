subroutine sup_getste(a,dtinv)
!-----------------------------------------------------------------------
! NAME 
!    nsi_getste
! DESCRIPTION
!    This routine computes the time step size for the incompressible NS
!    equation.
!-----------------------------------------------------------------------
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_Mesh
   use Mod_Element
   use Mod_PhysicalProblem   
   use Mod_NavierStokes, only : NavierStokesProblem   
   use Mod_GatherScatterDtcri
   use Mod_ThreeField
   use MPI   

   implicit none 
   class(ThreeFieldNSProblem), target :: a
   real(rp) :: dtinv
   
   real(rp), allocatable :: grvel(:,:)
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: gpvel(:)
   real(rp), allocatable :: eltem(:)
   class(FiniteElement), pointer :: e => NULL()
   real(rp) :: divel,acvit,acvis,acden,gpvno,dtmin,dtcri,hclen,rdinv,gradnorm
   real(rp) :: lambda,auxUnorm,dtcri2,dtmin2,dtuser,dtfactormin,dtfactormax
   integer(ip) :: ielem,nelem,ndime,irank,idime,jdime
   real(rp)    :: aux_dtinv(a%MPIsize)
   
   !todo multy materials
   integer(ip) :: imat=1
   
   !MPI
   integer, parameter :: mtag1 = 1, mtag2 = 2
   integer status(MPI_STATUS_SIZE)
   integer :: ierr,irequest
   
   call a%Timer%Total%Tic
   call a%Timer%Getste%Tic

   !Critical time step
   if(a%kfl_timei/=0.and.a%kfl_stead/=1) then
   
      !Dimensions and general variables
      call a%Mesh%GetNdime(ndime) 
      call a%Mesh%GetNelem(nelem)
      rdinv=1.0_rp/real(ndime)
      gpvno = 0.0_rp
      acvit = 0.0_rp
      dtmin = 1e6
      dtmin2= 1e6
      
      !Physical Parameters
      call a%GetPhysicalParameters(imat,acden,acvis)
      
      !Element Initialization
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sup_getste')
      
      call a%Memor%alloc(ndime,ndime,grvel,'grvel','sup_getste')
      call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','sup_getste')
      call a%Memor%alloc(ndime,gpvel,'gpvel','sup_getste')
      call a%Memor%alloc(e%mnode,eltem,'eltem','sup_getste')
      do ielem = 1,nelem
         
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)
         
         !Derivative (and detjm) at the center of gravity
         call e%elmdcg
         
         !Element length
         hclen=(e%weicg*e%detjm)**rdinv
         
         !Critical Time step computation
         dtcri=0.0_rp
         dtcri2=0.0_rp
         
         !viscoelastic parameters
         call a%IncrementalLambda(1,a%MatProp(1)%LawViParam(3),lambda)
         
         !lambda = a%MatProp(imat)%LawViParam(3)
         
          if (a%kfl_cotem ==1) then
            call e%gather(1,eltem,a%tempe)
            if (a%kfl_cotem_WLF==1) call sup_templaw_WLF(e, eltem, a%ReferenceTemp, a%c1_WLF, a%c2_WLF, a%MatProp(imat)%LawViParam, acvis, lambda)
            if (a%kfl_cotem_Arrhenius==1) call sup_templaw_Arrhenius(e, eltem, a%ReferenceTemp, a%alpha_Arrhenius, a%MatProp(imat)%LawViParam, acvis, lambda)
          endif 
            
         dtcri2 = dtcri2 + a%staco(1)*e%npol*e%npol/(2.0_rp*acvis)
         
         !Norm of the velocity at the center of gravity
         call e%gather(ndime,elvel,a%veloc(:,:,1))  
         call e%interpc(ndime,elvel,gpvel)
         !Velocity norm
         call vecnor(gpvel,e%ndime,gpvno,2)
         
         !Advection exists
         if (a%kfl_advec==1) then         
            !2u/h
            dtcri=dtcri+2.0_rp*gpvno/hclen         
         endif
         
         dtcri2=dtcri2 + a%staco(3)*(lambda/(2.0_rp*acvis))*(gpvno/hclen)
         
         !Smagorinsky, modify viscosity
         if (a%MatProp(imat)%lawvi < 0) then
            call e%gradient(ndime,elvel,grvel)
                  
            !second invariant of velocity gradient
            gradnorm=0.0_rp
            do idime=1,e%ndime
               do jdime=1,e%ndime
                  gradnorm=gradnorm + grvel(idime,jdime)&
                  *(grvel(idime,jdime))
               end do
            end do     
            gradnorm=0.5_rp*(gradnorm)**(0.5_rp)
            
            dtcri2=dtcri2 + a%staco(4)*(lambda/(acvis))*(gradnorm)            

         endif
         
         acvis=a%MatProp(imat)%visco+acvit
         
         ! 4*(mu/rho)/h^2
         dtcri=dtcri+4.0_rp*(acvis/acden)/(hclen*hclen)
         
         dtcri=1.0_rp/dtcri
         dtcri2=1.0_rp/dtcri2         
         
         dtmin=min(dtmin,dtcri)
         dtmin2=min(dtmin2,dtcri2)
      end do
      
      write(*,*) 'dtminmom'
      write(*,*) dtmin
      write(*,*) 'dtmincons'
      write(*,*) dtmin2  
      
      a%dtcri = min(dtmin,dtmin2)
      
      !Gather Dtcri from all processes, find minumum dtcri value, scatter it
      call php_GatherScatterDtcri(a)

      a%dtinv = 1.0_rp/(a%dtcri*a%safet)
      dtinv = a%dtinv
                     
      !Memory deallocation
      call a%Memor%dealloc(ndime,ndime,grvel,'grvel','sup_getste')
      call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','sup_getste')
      call a%Memor%dealloc(ndime,gpvel,'gpvel','sup_getste')
      call a%Memor%dealloc(e%mnode,eltem,'eltem','sup_getste')
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sup_getste')
   end if
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
end subroutine sup_getste
