subroutine plcd_output(a,itask) 
   !-----------------------------------------------------------------------
   !****f* SOLIDS/plcd_output
   ! NAME 
   !    plcd_output
   ! DESCRIPTION
   !    End of a SOLIDS time step 
   !    itask = 0  When timemarching is true. There is output or post-process
   !               of results if required.
   !    itask = 1  When timemarching is false. Output and/or post-process of
   !               results is forced if they have not been written previously.
   ! USES
   !    output
   !    a%FilePostpr%postpr
   ! USED BY
   !    plcd_endste (itask=0)
   !    plcd_turnof (itask=1)
   !***
   !-----------------------------------------------------------------------
   use typre
   use Mod_PLCD
   use Mod_Postpr
   use Mod_Stream
   use Mod_plcd_StrainGenerator
   use Mod_r1pElementAllocation
   implicit none
   class(PLCDProblem) :: a
   integer(ip)             :: itask
   integer(ip)             :: itime,istat,ifiel,ndime,aux4
   integer(ip), save       :: dopost(25)
   real(rp)                :: dummr
   
   type(r2p), allocatable :: PStress(:)
   
   procedure(), pointer :: GetStrain => NULL()
   procedure(), pointer :: AddSphericalComponent => NULL()
   integer(ip) :: vsize
   
   integer(ip) :: igaus,ielem,nelem


   interface
      subroutine plcd_PostprocessMaterialData(b,iiter_char)
         import
         implicit none
         class(PLCDProblem), target :: b
         character(5) :: iiter_char
      end subroutine
   end interface
   
   
   !Decide which postprocesses need to be done
   call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)
   
   !Do the actual postprocess
   !Displacements
   if (dopost(1) == 1) then
      call a%FilePostpr%postpr(a%Displacement(:,:,1),'Displacement',a%istep,a%ctime,a%Mesh)
   endif
   !Stresses
   if (dopost(2) == 1) then
      call a%FilePostpr%postgp(a%Stress,'Stress',a%istep,a%ctime,a%Mesh,'voigtT')
      if (a%UseUPFormulation) call a%FilePostpr%postpr(a%Pressure(:,1),'Pressure',a%istep,a%ctime,a%Mesh)
      
      
      
      
!       !PStresses
!       call a%Mesh%GetNdime(ndime)
!       call AllocateGetStrain(ndime,GetStrain,vsize,AddSphericalComponent)
!       call AllocateR2pElement(a%Mesh,vsize,PStress,a%Memor,'stress')
!       call a%Mesh%GetNelem(nelem)
!       do ielem = 1,nelem
!          do igaus = 1,size(a%Stress(ielem)%a,2)
!             call a%TopologicalDerivative%PT%ApplyPTensor(1.0_rp,a%Stress(ielem)%a(:,igaus),PStress(ielem)%a(:,igaus))
!          enddo
!       enddo
!       
!       call a%FilePostpr%postgp(PStress,'PStress',a%istep,a%ctime,a%Mesh,'voigtT')
      
      
      
      
   endif
   !Strains
   if (dopost(3) == 1) then
      call a%FilePostpr%postgp(a%Strain,'Strain',a%istep,a%ctime,a%Mesh,'voigtT')
   endif   
   !Material Data
   if (dopost(4) == 1) then
      call plcd_PostprocessMaterialData(a,'     ')
   endif
   !Subscales
   if (dopost(5) == 1 .and. a%UseUPFormulation) then
      call a%Mesh%GetNdime(ndime)
      call a%FilePostpr%postgp(a%UPSubscales,'Subscales',a%istep,a%ctime,a%Mesh)
   endif
   
   !TopologicalDerivative
   if (a%kfl_TopologyOptimization == 1) write(*,*) 'hey SIMP output should be fixed'
   if (dopost(6) == 1 .and. a%kfl_TopologyOptimization == 2) then
      call a%TopologicalDerivative%PostprocessChi(a%FilePostpr,a%itera,a%istep,a%ctime)
   endif
   
   !Velocities
   if (dopost(7) == 1 .and. a%kfl_TransientProblem /= 0) then
      call a%FilePostpr%postpr(a%Velocity(:,:,1),'Velocity',a%istep,a%ctime,a%Mesh)
   endif
   
   !Accelerations
   if (dopost(8) == 1 .and. a%kfl_TransientProblem /= 0) then
      call a%FilePostpr%postpr(a%Acceleration(:,:,1),'Acceleration',a%istep,a%ctime,a%Mesh)
   endif
     
   ! Tracking of points.
   if(a%nptra > 0) then
      call a%plcd_outtpo()
   end if   
   
   !Fluid Traction
   if (a%kfl_FSI == 1) then
       call a%FilePostpr%postpr(a%fluidtraction(:,:),'FluidTraction',a%istep,a%ctime,a%Mesh)
       call a%FilePostpr%postpr(a%btraction(:,:),'SolidTraction',a%istep,a%ctime,a%Mesh)
   endif
   

end subroutine plcd_output
