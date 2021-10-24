subroutine sup_output(a,itask)
  !-----------------------------------------------------------------------
  !****f* Nstinc/sup_output
  ! NAME 
  !    sup_output
  ! DESCRIPTION
  !    End of a NSTINC time step 
  !    itask = 0  When timemarching is true. There is output or post-process
  !               of results if required.
  !    itask = 1  When timemarching is false. Output and/or post-process of
  !               results is forced if they have not been written previously.
  ! USES
  !    output
  !    a%FilePostpr%postpr
  ! USED BY
  !    sup_endste (itask=1)
  !    sup_turnof (itask=2)
  !***
  !-----------------------------------------------------------------------
  use typre
  use Mod_Postpr
  use Mod_ThreeField
  use Mod_Stream
  implicit none
  class(ThreeFieldNSProblem) :: a
  integer(ip)                :: itask
  integer(ip)             :: itime,istat,ifiel,ndime,ntens
  real(rp)                :: dummr
  !Todo Multy materials
  integer(ip)             :: imat=1
  integer(ip), save       :: dopost(49)
  
  call a%Mesh%GetNdime(ndime)
  ntens=(ndime-1)*(ndime-1)+2   
  select case(itask)

  !End of a time step.
  
   case(0)
      a%pos_alrea=0
     !------------------------------------------------------------------------- 
     ! Postprocess at a given time
   
      if(a%kfl_exacs/=0 .and. a%kfl_timei/=0) then
         call a%FilePostpr%postpr(a%exsigma(:,:,1),'exSigma',a%istep,a%ctime,a%Mesh)
         call a%FilePostpr%postpr(a%exveloc(:,:,1),'exVeloc',a%istep,a%ctime,a%Mesh)         
         call a%SpecificExaerr    
      end if       
      
  case(1)

      if(a%kfl_exacs/=0 .and. a%kfl_timei==0) then      
         call a%SpecificExaerr    
      end if  

  end select
  
  !Decide which postprocesses need to be done
  call SetupDoPost(itask,a%istep,size(a%npp_stepi),a%npp_inits,a%npp_stepi,a%pos_alrea,dopost)
   
  !Do the actual postprocess

  !non Newtonian case or temperature coupling case
  !if (dopost(5) == 1) then
  !    if(a%MatProp(imat)%lawvi/= 0)then
  !       call a%FilePostpr%postgp(a%viscarray(:),'nnvis',a%istep,a%ctime,a%Mesh)
  !    endif 
  !endif
  
   !Stress
   if (dopost(6) == 1) then
      if (a%LogFormulation==0) then
         call a%FilePostpr%postpr(a%sigma(:,:,1),'Sigma',a%istep,a%ctime,a%Mesh,'voigtT')
      elseif (a%LogFormulation==1) then
         call a%FilePostpr%postpr(a%sigmaold(:,:),'Sigma',a%istep,a%ctime,a%Mesh,'voigtT')
       end if

   endif
   
   if (dopost(29)==1 .and. a%LogFormulation==1) then
      call a%FilePostpr%postpr(a%sigma(:,:,1),'Psi',a%istep,a%ctime,a%Mesh,'voigtT')
    end if    
    
   if (dopost(30)==1 .and. a%LogFormulation==1) then
      call a%FilePostpr%postpr(a%sigma(:,:,1),'PsiReal',a%istep,a%ctime,a%Mesh,'voigtT')
    end if  
   
   !Residual Projections
   if (dopost(19) == 1) then
      call a%FilePostpr%postpr(a%repro(1:ntens,:),'Proj_residual_S',a%istep,a%ctime,a%Mesh)
   end if
     
   !Subscales
   if (dopost(14) == 1 .and. a%kfl_tacsg==1) then    
      if(a%kfl_repro>=2)then
        call a%FilePostpr%postgp(a%vesgs3,'Velocity_SGS3',a%istep,a%ctime,a%Mesh)
      end if
      
      if(a%MatProp(imat)%lawvi< 0)then
         call a%FilePostpr%postgp(a%sisgs,'Stress_SGS1',a%istep,a%ctime,a%Mesh,'voigtT')
         if(a%kfl_repro == 4)then
            call a%FilePostpr%postgp(a%sisgs2,'Stress_SGS2',a%istep,a%ctime,a%Mesh,'voigtT')
            call a%FilePostpr%postgp(a%sisgs3,'Stress_SGS3',a%istep,a%ctime,a%Mesh,'voigtT')
         end if
      end if
   end if
       
   !RelaxTime
   if (dopost(23) == 1) then
      if (a%kfl_cotem_WLF ==1 .or. a%kfl_cotem_Arrhenius==1) call a%FilePostpr%postgp(a%lambdarray(:),'RelaxTime',a%istep,a%ctime,a%Mesh)
   end if    
  
   if (dopost(28) == 1) then
      call a%FilePostpr%postgp(a%tau_mom(:),'Tau_mom',a%istep,a%ctime,a%Mesh)
      call a%FilePostpr%postgp(a%tau_sig(:),'Tau_sig',a%istep,a%ctime,a%Mesh)
   end if  
  
  
end subroutine sup_output

