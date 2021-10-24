subroutine nsc_output(a,itask)
  !-----------------------------------------------------------------------
  !****f* Nscomp/nsc_output
  ! NAME 
  !    nsc_output
  ! DESCRIPTION
  !    End of a NSCOMP time step 
  !    itask = 0  When timemarching is true. There is output or post-process
  !               of results if required.
  !    itask = 1  When timemarching is false. Output and/or post-process of
  !               results is forced if they have not been written previously.
  ! USES
  !    output
  !    a%FilePostpr%postpr
  ! USED BY
  !    nsc_endste (itask=0)
  !    nsc_turnof (itask=1)
  !***
  !-----------------------------------------------------------------------
  use typre
  use Mod_NSCompressible
  use Mod_Postpr
  use Mod_Stream
  implicit none
  class(NSCompressibleProblem) :: a
  integer(ip)                :: itask

  
  real(rp), allocatable   :: error(:)
  real(rp)                :: TotalEstimatedError
  integer(ip)             :: itime,ndime,nelem
  integer(ip), save       :: pos_alrea(25)

   interface
      subroutine nsc_outtpo(a)
         import NSCompressibleProblem
         implicit none
         class(NSCompressibleProblem) :: a
      end subroutine
   end interface


  select case(itask)

  !End of a time step.
  case(0)

     pos_alrea=0

     ! Tracking of points.
     if(a%nptra > 0) then
        call nsc_outtpo(a)
     end if

     ! Postprocess at a given time step
     if(a%istep>=a%npp_inits.and.a%npp_stepi(1)>0) then     
        if(mod(a%istep,a%npp_stepi(1))==0) then
           pos_alrea(1)=1
        end if
     end if
     
     if(a%istep>=a%npp_inits.and.a%npp_stepi(2)>0) then     
        if(mod(a%istep,a%npp_stepi(2))==0) then
           pos_alrea(2)=1
        end if
     end if

     if(a%istep>=a%npp_inits.and.a%npp_stepi(3)>0) then     
        if(mod(a%istep,a%npp_stepi(3))==0) then
           pos_alrea(3)=1
        end if
     end if
     
     if(a%istep>=a%npp_inits.and.a%npp_stepi(4)>0) then     
        if(mod(a%istep,a%npp_stepi(4))==0) then
           pos_alrea(4)=1
        end if
     end if
     if(a%istep>=a%npp_inits.and.a%npp_stepi(5)>0) then     
        if(mod(a%istep,a%npp_stepi(5))==0) then
           pos_alrea(5)=1
        end if
     end if
     if(a%istep>=a%npp_inits.and.a%npp_stepi(6)>0) then     
        if(mod(a%istep,a%npp_stepi(6))==0) then
           pos_alrea(6)=1
        end if
     end if
     if(a%istep>=a%npp_inits.and.a%npp_stepi(7)>0) then     
         if(mod(a%istep,a%npp_stepi(7))==0) then     
            pos_alrea(7)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(9)>0) then     
         if(mod(a%istep,a%npp_stepi(9))==0) then
            pos_alrea(9)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(10)>0) then     
         if(mod(a%istep,a%npp_stepi(10))==0) then
            pos_alrea(10)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(11)>0) then     
         if(mod(a%istep,a%npp_stepi(11))==0) then
            pos_alrea(11)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(12)>0) then     
         if(mod(a%istep,a%npp_stepi(12))==0) then
            pos_alrea(12)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(13)>0) then     
         if(mod(a%istep,a%npp_stepi(13))==0) then
            pos_alrea(13)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(14)>0) then     
         if(mod(a%istep,a%npp_stepi(14))==0) then
            pos_alrea(14)=1
         end if
      end if     
      if(a%istep>=a%npp_inits.and.a%npp_stepi(15)>0) then     
         if(mod(a%istep,a%npp_stepi(15))==0) then
            pos_alrea(15)=1
         end if
      end if     
 
    ! write(*,*) 'nsc_outpu: heat flux not ready for a%FilePostpr%postprocess'

     ! Postprocess at a given time
     if(a%ctime>=a%pos_tinit.and.pos_alrea(1)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,1)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,1)>0.0_rp) then
               pos_alrea(1)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(2)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,2)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,2)>0.0_rp) then
               pos_alrea(2)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(3)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,3)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,3)>0.0_rp) then
               pos_alrea(3)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(4)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,4)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,4)>0.0_rp) then
               pos_alrea(4)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(5)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,5)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,5)>0.0_rp) then
               pos_alrea(5)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(6)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,6)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,6)>0.0_rp) then
               pos_alrea(6)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(7)==0) then  
         do itime=1,10
            if(   abs(a%pos_times(itime,7)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,7)>0.0_rp) then
               pos_alrea(7)=1
            end if
         end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(9)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,9)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,9)>0.0_rp) then
               pos_alrea(9)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(10)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,10)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,10)>0.0_rp) then
               pos_alrea(10)=1
           end if
        end do
     end if
     if(a%ctime>=a%pos_tinit.and.pos_alrea(11)==0) then  
        do itime=1,10
           if(   abs(a%pos_times(itime,10)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,10)>0.0_rp) then
               pos_alrea(11)=1
            end if 
         end do
      end if     
      if(a%ctime>=a%pos_tinit.and.pos_alrea(12)==0) then  
         do itime=1,10
            if(   abs(a%pos_times(itime,12)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,12)>0.0_rp) then
               pos_alrea(12)=1
            end if
         end do
      end if
      if(a%ctime>=a%pos_tinit.and.pos_alrea(13)==0) then  
         do itime=1,10
            if(   abs(a%pos_times(itime,13)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,13)>0.0_rp) then
               pos_alrea(13)=1
            end if
         end do
      end if
      if(a%ctime>=a%pos_tinit.and.pos_alrea(14)==0) then  
         do itime=1,10
            if(   abs(a%pos_times(itime,14)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,14)>0.0_rp) then
               pos_alrea(14)=1
            end if
         end do
      end if
      if(a%ctime>=a%pos_tinit.and.pos_alrea(15)==0) then  
         do itime=1,10
            if(   abs(a%pos_times(itime,15)-a%ctime)<(0.5_rp*a%dtime).and.a%pos_times(itime,15)>0.0_rp) then
               pos_alrea(15)=1
            end if
         end do
      end if

      if(a%kfl_exacs/=0) then      
         if (a%RefinerErrorEstimator == 'SUBSC') then
            call a%Mesh%GetNelem(nelem)
            call a%Memor%alloc(nelem,error,'error','nsc_output')
            call a%SpecificExactSolRefCriteria(error)
            call a%FilePostpr%postgp(error,'ExactError',a%istep,a%ctime,a%Mesh)
            call a%Memor%dealloc(nelem,error,'error','nsc_output')
         end if
      end if       
  !End of the run.
  case(1)

      if(a%npp_stepi(1)/=0)then
         if(pos_alrea(1)==1) then
            pos_alrea(1)=0
         else 
            pos_alrea(1)=1
         end if
      end if
      if(a%npp_stepi(2)/=0)then
         if(pos_alrea(2)==1) then
            pos_alrea(2)=0
         else 
            pos_alrea(2)=1
         end if
      end if      
      if(a%npp_stepi(3)/=0)then
         if(pos_alrea(3)==1) then
            pos_alrea(3)=0
         else 
            pos_alrea(3)=1
         end if
      end if      
      if(a%npp_stepi(4)/=0)then
         if(pos_alrea(4)==1) then
            pos_alrea(4)=0
         else 
            pos_alrea(4)=1
         end if
      end if      
      if(a%npp_stepi(5)/=0)then
         if(pos_alrea(5)==1) then
            pos_alrea(5)=0
         else 
            pos_alrea(5)=1
         end if
      end if
      if(a%npp_stepi(6)/=0)then
         if(pos_alrea(6)==1) then
            pos_alrea(6)=0
         else 
            pos_alrea(6)=1
         end if
      end if      
      if(a%npp_stepi(7)/=0)then
         if(pos_alrea(7)==1) then
            pos_alrea(7)=0
         else 
            pos_alrea(7)=1
         end if
      end if      
      if(a%npp_stepi(9)/=0)then
         if(pos_alrea(9)==1) then
            pos_alrea(9)=0
         else 
            pos_alrea(9)=1
         end if
      end if
      if(a%npp_stepi(10)/=0)then
         if(pos_alrea(10)==1) then
            pos_alrea(10)=0
         else 
            pos_alrea(10)=1
         end if
      end if
      if(a%npp_stepi(11)/=0)then
         if(pos_alrea(11)==1) then
            pos_alrea(11)=0
         else 
            pos_alrea(11)=1
         end if
      end if
      if(a%npp_stepi(12)/=0)then
         if(pos_alrea(12)==1) then
            pos_alrea(12)=0
         else 
            pos_alrea(12)=1
         end if
      end if
      if(a%npp_stepi(13)/=0)then
         if(pos_alrea(13)==1) then
            pos_alrea(13)=0
         else 
            pos_alrea(13)=1
         end if
      end if
      if(a%npp_stepi(14)/=0)then
         if(pos_alrea(14)==1) then
            pos_alrea(14)=0
         else 
            pos_alrea(14)=1
         end if
      end if
      if(a%npp_stepi(15)/=0) then      
         call a%FinalizeStats
         call a%FilePostpr%postpr(a%rmsd,'RMSDensf',a%istep,a%ctime,a%Mesh)
         call a%FilePostpr%postpr(a%rmsv,'RMSVeloc',a%istep,a%ctime,a%Mesh)
         call a%FilePostpr%postpr(a%rmsp,'RMSPress',a%istep,a%ctime,a%Mesh)
         call a%FilePostpr%postpr(a%rmst,'RMSTempe',a%istep,a%ctime,a%Mesh)
         call a%FilePostpr%postpr(a%turbi,'TurbulenceIntensity',a%istep,a%ctime,a%Mesh)
         if(pos_alrea(15)==1) pos_alrea(15)=0
      endif
      
      if(a%kfl_exacs/=0) then      
         call a%SpecificExaerr    
      end if       

      if (a%RefinerErrorEstimator == 'SUBSC') then
         if(pos_alrea(12)==1) then
            pos_alrea(12)=0
         else 
            pos_alrea(12)=1
         end if
         call a%Mesh%GetNelem(nelem)
         call a%Memor%alloc(nelem,error,'error','nsc_output')
         call a%SpecificSubscalesRefCriteria(error,TotalEstimatedError)
         call a%FilePostpr%postgp(error,'Error',a%istep,a%ctime,a%Mesh)
         call a%Memor%dealloc(nelem,error,'error','nsc_output')
      end if

  end select

  call a%Mesh%GetNdime(ndime)
  !Velocity
  if(pos_alrea(1)==1) then
     call a%FilePostpr%postpr(a%veloc(:,:,1),'Velocity',a%istep,a%ctime,a%Mesh)
  end if
  !Pressure
  if(pos_alrea(2)==1) then
     call a%FilePostpr%postpr(a%press(:,1),'Pressure',a%istep,a%ctime,a%Mesh)
  end if      
  !Temperature
  if(pos_alrea(3)==1) then
     call a%FilePostpr%postpr(a%tempe(:,1),'Temperature',a%istep,a%ctime,a%Mesh)
  end if      
  !Density
  if(pos_alrea(4)==1) then
     call a%FilePostpr%postpr(a%densf(:,1),'Density',a%istep,a%ctime,a%Mesh)
  end if      
  !Momentum
  if(pos_alrea(5)==1) then
     call a%FilePostpr%postpr(a%momen(:,:,1),'Momentum',a%istep,a%ctime,a%Mesh)
  end if
  !Energy
  if(pos_alrea(6)==1) then
     call a%FilePostpr%postpr(a%energ(:,1),'Energy',a%istep,a%ctime,a%Mesh)
  end if      
  !streamlines
  if(pos_alrea(7)==1) then
     call stream(a%veloc(:,:,1),ndime,a%istep,a%ctime,a%Mesh,a%Memor,a%FilePostpr,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPISize)
  end if      
  if(pos_alrea(9)==1) then      
     if(a%kfl_shock > 0)then
        call a%FilePostpr%postgp(a%scdiffarray(1,:),'arvis',a%istep,a%ctime,a%Mesh)
        call a%FilePostpr%postgp(a%scdiffarray(2,:),'artco',a%istep,a%ctime,a%Mesh)
     end if 
  end if
  if(pos_alrea(10)==1) then      
     call a%FilePostpr%postgp(a%sgdiffarray(1,:),'sgvis',a%istep,a%ctime,a%Mesh)
     call a%FilePostpr%postgp(a%sgdiffarray(2,:),'sgtco',a%istep,a%ctime,a%Mesh)
  end if
  if(pos_alrea(11)==1) then      
     if(a%kfl_outfm>0)then
!        call a%FilePostpr%postpr(a%forfl(:,:),'FORFL',a%istep,a%ctime,a%Mesh)
!        call a%FilePostpr%postpr(a%momfl(:,:),'MOMFL',a%istep,a%ctime,a%Mesh)
!        call a%FilePostpr%postpr(a%hfxfl(:),'HFXFL',a%istep,a%ctime,a%Mesh)
     endif
  end if
  if(pos_alrea(12)==1) then      
     if(a%kfl_trasg>0)then
        call a%FilePostpr%postgp(a%cosgs,'COSGS',a%istep,a%ctime,a%Mesh,'scalar')
        call a%FilePostpr%postgp(a%mosgs,'MOSGS',a%istep,a%ctime,a%Mesh)
        call a%FilePostpr%postgp(a%ensgs,'ENSGS',a%istep,a%ctime,a%Mesh,'scalar')
     end if 
  end if
  if(pos_alrea(13)==1) then      
        call a%FilePostpr%postgp(a%coResGP,'CORES',a%istep,a%ctime,a%Mesh)
        call a%FilePostpr%postgp(a%moResGP,'MORES',a%istep,a%ctime,a%Mesh)
        call a%FilePostpr%postgp(a%enResGP,'ENRES',a%istep,a%ctime,a%Mesh)
  end if
  if(pos_alrea(14)==1) then      
        if(a%kfl_repro>0)then
           call a%FilePostpr%postpr(-a%repro(1,:),'DERPR',a%istep,a%ctime,a%Mesh)
           call a%FilePostpr%postpr(-a%repro(2:ndime+1,:),'MORPR',a%istep,a%ctime,a%Mesh)
           call a%FilePostpr%postpr(-a%repro(ndime+2,:),'ENRPR',a%istep,a%ctime,a%Mesh)
        end if 
  end if
  if(pos_alrea(15)==1) then    
     call a%FilePostpr%postpr(a%timad(:,1),'TimeAveraDensf',a%istep,a%ctime,a%Mesh)
     call a%FilePostpr%postpr(a%timav(:,:,1),'TimeAveraVeloc',a%istep,a%ctime,a%Mesh)
     call a%FilePostpr%postpr((a%timap(:,1)),'TimeAveraPress',a%istep,a%ctime,a%Mesh)
     call a%FilePostpr%postpr((a%timat(:,1)),'TimeAveraTempe',a%istep,a%ctime,a%Mesh)
  endif

end subroutine nsc_output
