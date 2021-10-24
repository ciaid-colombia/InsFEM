module Mod_sld_BaseElmopeRoutines
   use Mod_php_SetTimeIntegrator
   use Mod_CauchyElement
   use Mod_sld_BaseElmope
   use Mod_sld_NonLinearDerivatives
   use Mod_sld_ExternalForces
   use Mod_sld_HangingNodes
   implicit none

   contains

   subroutine SetPointersGeneral
      use typre
      implicit none
      integer(ip) :: kfl_nonlinear
            
      !Set All Pointers To NULLSUB (in Mod_sld_BaseElmope)
      call ResetProcedureComposition
      call SetPointersAndHooksToNULLSUB

      !Non-linear elements
      call SetPointersNonLinearDerivatives(0)
      call SetPointersHangingNodes(0)

      !Non linear elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
          call SetPointersNonLinearDerivatives(1)
          call ConcatenateProcedures(ProcHook%PreGauss,InGaussVolumesNonLinear)
      end if

      call ConcatenateProcedures(ProcHook%InGauss,CalculateVolume)

      call SetPointersHangingNodes(1)

      !Mod deallocs
      call SetPointersNonLinearDerivatives(100)
      call SetPointersHangingNodes(100)

   end subroutine SetPointersGeneral

   subroutine SetPointersPointForces

      !Point forces
      if(a%kfl_sourcesPresent) then
          call setPointForcesFlag

          if(kfl_setPointForce==1) then
              call ConcatenateProcedures(ProcHook%Finalizations, setPointForces)
          endif
      endif

   end subroutine SetPointersPointForces


   subroutine SetRHSPointers

      if(a%kfl_timei==1) then
          call ConcatenateProcedures(ProcHook%preAssemblyLhs,completeMass)
          call ConcatenateProcedures(ProcHook%DynamicForce,dynamicForce)
          if (a%kfl_tsche_2nd_current == 'NEWMA') then
              call ConcatenateProcedures(ProcHook%Interpolates,InterpolateNewmark)
          endif
      endif
      call ConcatenateProcedures(ProcHook%InGaussVec,calculateRHS)
      call ConcatenateProcedures(ProcHook%preAssemblyRhs,addExternalForces)

   end subroutine SetRHSPointers

   subroutine addExternalForces
      implicit none
      integer(ip) :: u1,uf,s1,sf,p1,bc

      call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)

      elrhu = elrhu + a%extForce(ielem)%a(u1:uf,:)

  end subroutine

   subroutine calculateRHS
      implicit none

      !Compute contributions to RHS :
      call ProcPointer%sld_elmrhu(e,dvol,elext,elrhu)

  end subroutine

   subroutine setPointForcesFlag
      implicit none

      kfl_setPointForce=0

      !If constant force
      if (a%kfl_constPointForce==1) then
          kfl_setPointForce=1 

      !If time dependant force
      elseif (a%kfl_PointForceTime>=a%ctime) then
          kfl_setPointForce=1

      endif
   end subroutine

   subroutine setPointForces
      implicit none
      integer(ip) :: npoin,ndime
      real(rp), pointer:: pointRhs(:,:) => NULL()

      call a%Mesh%GetNpoin(npoin)
      call a%Mesh%GetNdime(ndime)

      call a%Mesh%Global2Local(a%pwSourceSize,int(a%pwSource(1,:)),a%pwSourceId)
      pointRhs => a%pwSource(2:a%pwSourceDim,:)
      call a%LinearSystem%AssemblyPointRhs(a%pwSourceSize,a%pwSourceId,pointRhs)

   end subroutine

   subroutine InitTimeIntegration
      implicit none

      call php_SetTimeIntegrator(a,Integrator,LHSdtinv2,nsteps)
      if (a%kfl_tsche_2nd_current == 'NEWMA') then
          call php_SetNewmarkCoefficients(a,Integrator,a%beta,a%gamma)
      endif
      call php_GetTimeSchemeDerivative(a,Integrator,tsch_deriv)

   end subroutine

   subroutine dynamicForce
      implicit none
      
      !not used in Newton Raphson
      call sld_TimeInt2extF(e,nsteps,Integrator,a%dtinv2,gpdisp(:,2:nsteps),elext,densi)
      
   end subroutine

   subroutine staticResidual
      implicit none
      
      elrhu = elrhu + a%extForce(ielem)%a(:,:) - intForce
      
   end subroutine

   subroutine dynamicResidual
      implicit none
      real(rp) :: massAccel(e%ndime,e%pnode)
      
      call matvecmult(e%ndime,e%pnode,massmat,elaccel(:,:,1),massAccel)

      elrhu = elrhu - massAccel 
      
   end subroutine

   subroutine buildMassMat
      implicit none
      integer(ip) :: nd
      
      call a%Mesh%GetNdime(nd)

      call ProcHook%PhysicalProp
      call sld_buildMassMatrix(e,dvol,nd,densi,massmat,1.0_rp)

  end subroutine

   subroutine completeMass
      implicit none

      !Add integration terms u_n-1
      massmat = LHSdtinv2*massmat

   end subroutine

   subroutine InGaussResetExternal
       implicit none

       !reset
       elext=0.0_rp

   end subroutine

   subroutine CalculateVolume
       implicit none

       dvol = 0.0_rp
       dvol = e%weigp(e%igaus)*e%detjm
       dvolt0 = dvol + dvolt0

   end subroutine

   !NonLinear Elements   
   subroutine InGaussVolumesNonLinear
      implicit none
      dvolt0=0.0_rp
   end subroutine

   subroutine AssembleLinearSystem
       implicit none

       !For adaptivity we share hanging nodes
       call ProcHook%PreDirichlet

       !Dirichlet Boundary Conditions, set residual for dirichlet nodes
       call sld_elmdir(a,e,elmat,elrhs)

       call a%LinearSystem%Assembly(e,elmat,elrhs)

   end subroutine

   subroutine CloseEigenProblem
       implicit none

       !If problem is linear, this needs to be done only once
       !Fix set linear flag and stuff
       a%kfl_eigenter= .false.

   end subroutine

   subroutine sld_elemLoop
      implicit none
   
      ielem = 1
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_ElmLoop')
   
      !Allocate Arrays in BaseElmope
      call ProcHook%AllocateArrays
      
      call ProcHook%PhysicalProp
      call ProcHook%Initializations

      call a%Mesh%GetNelem(nelem)

      elements : do ielem = 1,nelem
         !Load Element
         a%kfl_foldedElements=.false.
         call a%Mesh%ElementLoad(ielem,e)  
   
         call ProcHook%OnIeltyChange
   
         !ElmatsToZero
         call ProcHook%ResetArrays
   
         call ProcHook%Gathers
   
         !Cartesian derivatives and Jacobian at center of gravity
         call e%elmdcg
   
         !Element length at center of gravity
         call e%elmlen
   
         dvolt0 = 0.0_rp
   
         call ProcHook%PhysicalProp
   
         !Build Jacobian, A(d) matrix for the latest state of the body
         gauss_points: do igaus=1,e%pgaus
   
             e%igaus = igaus
   
             call ProcHook%PreGauss
   
             call ProcHook%InGauss
   
             !Interpolate         
             call ProcHook%Interpolates
   
             call ProcHook%PrePostInterpolates

             call ProcHook%PostInterpolates

             call ProcHook%PostPostInterpolates

             call ProcHook%EndGauss
   
         enddo gauss_points
   
         call ProcPointer%PostGauss
   
         call ProcHook%PostGauss
   
         call ProcHook%preAssemblyLhs
   
         call ProcHook%AssemblyLhs

         call ProcHook%postAssemblyLhs
   
         call ProcHook%preAssemblyRhs
   
         call ProcHook%AssemblyRhs

         call ProcHook%postAssemblyRhs
   
         call ProcHook%ToLinearSystem
   
         if(a%kfl_foldedElements) then
             exit
         endif
   
      enddo elements

      call ProcHook%Finalizations

      call ProcHook%DeallocateArrays
   
      !DeallocateElement
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_ElmLoop')
   
   end subroutine sld_elemLoop

end module Mod_sld_BaseElmopeRoutines
