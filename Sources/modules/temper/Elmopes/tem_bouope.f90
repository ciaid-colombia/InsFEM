module Mod_tem_bouope
   use typre
   use Mod_Element
   use Mod_Temperature
   use Mod_tem_BaseElmope
   use Mod_tem_NonLinearDerivatives
   use Mod_tem_HangingNodes
   use Mod_tem_DustTransport
   use Mod_tem_Advection
   implicit none
   
   real(rp), external :: funcre
   real(rp) :: updbcn
   
contains

   subroutine SetPointers
   
      integer(ip) :: kfl_nonlinear,nelty,kfl_HangingNodes
      
      !External Procedures
      procedure() :: NULLSUB
      
      call ResetProcedureComposition
      
      call SetPointersAndHooksToNULLSUB
      
      call SetPointersNonLinearDerivatives(0)
      call SetPointersHangingNodes(0)
      call SetPointersAdvectionVelocity(0)
      call SetPointersDustTransportBoundary(0)
      
      !-------------------------------------------------------------------
      !Non-Linear Elements
      call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
      if (kfl_nonlinear == 1) then
         call SetPointersNonLinearDerivatives(1)
      endif
      
      !HangingNodes
      call SetPointersHangingNodes(1)
      
      !DustTransport 
      call SetPointersDustTransportBoundary(1)
      
      
      !If more than one element type, then pointers should be reset when ielty changes!
      call a%Mesh%GetNelty(nelty)
      if (nelty > 1) ProcHook%OnIeltyChange => OnIeltyChange
      
      call SetPointersNonLinearDerivatives(100)
      call SetPointersHangingNodes(100)
      call SetPointersAdvectionVelocity(100)
      call SetPointersDustTransportBoundary(100)
   end subroutine
   
   subroutine OnIeltyChange
      implicit none
      
      if (e%ielty /= ielty0) then
         call SetPointers
         ielty0 = e%ielty
      endif
   end subroutine
   
end module


subroutine tem_bouope(TempeProblem)
   use typre
   use Mod_Temperature
   use Mod_tem_bouope
   use Mod_php_Elmdir
   implicit none
   class(TemperatureProblem), target :: TempeProblem
   integer(ip) :: jnode
   
   
   a=>TempeProblem  
   

   !Set the Pointers for execution
   ielem = 1
   call SetPointers
   
   !Initializations
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNboun(nboun)
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','tem_bouope')
   
   call a%Memor%alloc(1,e%mnode,1,e%mnode,elmat,'elmat','tem_bouope')
   call a%Memor%alloc(1,e%mnode,elrhs,'elrhs','tem_bouope')
   
    !Advection
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%alloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%alloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%alloc(e%mnode,AGradV,'AGradV','tem_elmope')
   gpvno = 0.0_rp
   
   !Physical Parameters
   call a%GetPhysicalParameters(1,acden,acsph,actco,acrea,acsou)
   acvis = actco/acsph
   
   !Hook
   call ProcHook%Initializations
   
   ! Loop over boundaries
   boundaries: do iboun=1,nboun
      if (a%NumberofMaterials > 1) then
         call a%Mesh%GetBoundaryIelem(iboun,ielem)
         !Physical Parameters
         call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
         acvis = actco/acsph
         acrcp = acrea/acsph
      endif
      
      if (a%kfl_fixbo(iboun) /= 0) then
         !Load Element
         call a%Mesh%BoundaryLoad(iboun,e)

         !Hook
         call ProcHook%OnIeltyChange
      
         !Hook
         call ProcHook%PreGauss
         
         !Initialize
         elmat=0.0_rp
         elrhs=0.0_rp
         
         !Gathers
         !Hook
         call ProcHook%Gathers
         
         call e%elmdcg
         
         !Element length at center of gravity
         call e%elmlen
         
         !Gauss Point Loop
         gauss_points : do igaub = 1,e%pgaub
            e%igaub = igaub
            
            tract = 0.0_rp
            
            !Hook
            call ProcHook%InGauss

            !Calculate exterior Normal
            !Derivatives at the boundary
            call e%elmderb
            call e%bounor
         
            dsurf=e%weigb(e%igaub)*e%eucta
        
            !Hook
            call ProcHook%Interpolates
            
            !Pointer 
            call ProcPointer%VelocityAndAgradV
      
            !PhysicalProperties
            acvis = actco/acsph
            !Hook
            call ProcHook%PhysicalProp
            
            !Newmann boundary condition
            if (a%kfl_fixbo(iboun) == 2) then
               if (a%kfl_conbc /= 1) then
                  if(a%kfl_funbo(iboun)/=0) then
                     updbcn=funcre(a%funpa(a%kfl_funbo(iboun))%a,&
                     a%kfl_funty(a%kfl_funbo(iboun),2),&
                     a%kfl_funty(a%kfl_funbo(iboun),1),a%bctime)
                  else
                     updbcn = 1.0_rp
                  endif   
               else
                  updbcn=1.0_rp
               end if
               tract = tract - a%bvnat(iboun)%a(1)*updbcn/acsph
               
            !Do nothing boundary condition   
            elseif (a%kfl_fixbo(iboun) == 6) then   
               do inodb = 1,e%pnodb
                  inode = e%lboel(inodb)
                  do jnode = 1,e%pnode
                     elmat(1,inode,1,jnode) = elmat(1,inode,1,jnode) - acvis*e%shapb(inodb,e%igaub)*dot_product(e%cartb(1:e%ndime,jnode),e%baloc(1:e%ndime,e%ndime))
                  enddo
               enddo
            else
               !call runend('tem_bouope: on boundaries boundary condition not implemented yet')
            endif
            
            call ProcHook%InGaussElmats
            
            !Traction to elrhs
            do inodb=1,e%pnodb
               inode = e%lboel(inodb)
               elrhs(1,inode) = elrhs(1,inode) +dsurf*tract*e%shapb(inodb,e%igaub)
            enddo
            
         enddo gauss_points
         
         !PreDirichlet Operations
         call ProcHook%PreDirichlet
      
         !Dirichlet Boundary conditions
         call php_elmdir(a,e,1_ip,a%ndofbc,a%ndofbcstart,1_ip,elmat,elrhs)
         
         !Assembly
         call a%LinearSystem%Assembly(e,elmat,elrhs)
         !call a%Mesh%AssemblyToSystem(e,1_ip,elmat,elrhs,a%LinearSystem)
         
      endif
      
   enddo boundaries
   
   !Finalizations
   !Hook
   call ProcHook%Finalizations
   
   call a%Memor%dealloc(1,e%mnode,1,e%mnode,elmat,'elmat','tem_bouope')
   call a%Memor%dealloc(1,e%mnode,elrhs,'elrhs','tem_bouope')
   
   !Advection
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpvel,'gpvel','tem_elmope')
   call a%Memor%dealloc(e%ndime,gpadv,'gpadv','tem_elmope')
   call a%Memor%dealloc(e%mnode,AGradV,'AGradV','tem_elmope')
   
   call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','tem_bouope')
   
end subroutine
