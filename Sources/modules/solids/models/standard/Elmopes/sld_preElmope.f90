subroutine sld_preElmope(SldProblem)
   use Mod_sld_BaseElmope
   use Mod_sld_Elmope

   implicit none
   class(SolidsProblem), target :: SldProblem
   integer(ip) :: ndime,sz
   
   a=>SldProblem
   
   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','sld_preElmope')
   !Allocate Arrays in BaseElmope
   call php_SetTimeIntegrator(a,Integrator,LHSdtinv2,nsteps)

   call AllocateBaseElmopeArrays
   call AllocateAssemblyArrays
   call AllocateSolidBase
   call AllocateDynamicArray

   call a%Mesh%GetNdime(ndime)
   sz        = (ndime*(ndime+1))/2
   LHSdtinv2 = 1.0_rp


   !Hook
   call ProcHook%Initializations

   call a%Mesh%GetNelem(nelem)


   elements : do ielem = 1,nelem
       !Load Element
       call a%Mesh%ElementLoad(ielem,e)  

       !ElmatsToZero
       elmat=0.0_rp
       massmat=0.0_rp
       elrhs=0.0_rp
       elrhu=0.0_rp

       call displacementGather
       call ProcHook%Gathers

       !Cartesian derivatives and Jacobian at center of gravity
       call e%elmdcg

       !Element length at center of gravity
       call e%elmlen

       !Physical Properties
       call a%GetPhysicalParameters(1.0_rp,sz,densi,c_elas,ielem)

       dvolt0 = 0.0_rp
       !Gauss Point Loop
       gauss_forces: do igaus=1,e%pgaus

           strain = 0.0_rp
           stress = 0.0_rp

           e%igaus = igaus
           !Time integration
           call ProcHook%InGauss

           dvol = 0.0_rp
           dvol = e%weigp(e%igaus)*e%detjm
           dvolt0 = dvol + dvolt0

           !Interpolate         
           call InterpolateGpDisplacements
           call ProcHook%Interpolates

           !Compute Elext, Temporal Derivatives
           elext=0.0_rp

           !Compute vector of external forces
           call ProcPointer%ExternalForces   
           !Compute contributions to RHS :
           call ProcPointer%sld_elmrhu(e,dvol,elext,elrhu)

           !K_mat and K_geo, and internal forces
           !call ProcHook%InGaussElmats

           !Mass matrix
           call ProcHook%DynamicMass

       enddo gauss_forces

       ! Assembly of jacobian A(d): mass and stiffness mats 
       elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) = elmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) + massmat(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)

       ! Assembly elrhu to elrhs
       elrhs(1:e%ndime,1:e%pnode) = elrhu(1:e%ndime,1:e%pnode) + elrhs(1:e%ndime,1:e%pnode) 

       !Dirichlet Boundary Conditions, set residual for dirichlet nodes
       call sld_elmdir(a,e,elmat,elrhs)

       !Assembly: a = M^-1 f
       call a%LinearSystem%Assembly(e,elmat,elrhs)
       
   enddo elements

   !Hook
   call ProcHook%Finalizations

   call DeallocateAssemblyArrays
   call DeallocateSolidBase
   call DeallocateDynamicArray
   call DeallocateBaseElmopeArrays

   !DeallocateElement
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','sld_preElmope')

   !We solve for initial accelerations
   call a%LinearSystem%Solve(a%unkno)

   a%accel(:,:,1)=a%unkno

   !call a%Mesh%ArrayCommunicator%GhostCommunicate(e%ndime,a%accel(:,:,1))    

end subroutine sld_preElmope
