subroutine lmn_cvgunk(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    low Mach number equations
   !-----------------------------------------------------------------------
   use MPI
   use typre
   use def_parame
   use Mod_LowMach
   implicit none
   class(LowMachProblem)   :: a
   integer(ip), intent(in) :: itask
   integer(ip), save       :: ipass=0
   real(rp)                :: rplmn,rvlmn,rtlmn,ripre,rivel,ritem,rilmn
   integer(ip)             :: npoinLocal,ndime,ierr
   
   select case(itask)
      
   !Check convergence of the inner iterations, always in the L2 norm:
   !|| u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
   case(1)
    
      call a%InnerResiduals(rvlmn,rplmn,rtlmn)
      rilmn = max(rvlmn,rtlmn,rplmn)
      a%rilmn = rilmn
      call MPI_ALLREDUCE(MPI_IN_PLACE,a%rilmn,1,MPI_REAL8,MPI_MAX,a%MPIcomm,ierr)

      if(((a%rilmn<a%cotol) .or. (a%itera>=a%maxit))) a%kfl_goite = 0       
      if (a%MPIrank == a%MPIroot) then
         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
            write(a%lun_outph,200)
            write(a%lun_outph,201)
            write(a%lun_nolin,300)
         end if

         write(a%lun_conve,101) a%istep,a%itera,a%ctime,rvlmn,rtlmn,rplmn

         write(a%lun_outph,202) a%istep,a%itera,a%ctime,& 
            a%tamin,a%tamax,a%tamea,&
            a%remin,a%remax,a%remea
         ! Write in the log file.
         if (a%kfl_flush == 1) then
            call flush(a%lun_outph)
            call flush(a%lun_conve)
         endif
      end if
  
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
   case(2)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPI(a%kfl_normc,ndime*npoinLocal,a%veloc(1,1,1),a%veloc(1,1,2),a%resid,zelmn,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(a%kfl_normc,npoinLocal,      a%press(1,1),  a%press(1,2),  a%resip,zelmn,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(a%kfl_normc,npoinLocal,      a%tempe(1,1),  a%tempe(1,2),  a%resit,zelmn,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      a%resid = max(a%resiv,a%resit,a%resip) 
      
   
   !Check residual of the time evolution, always in the L2 norm:
   !|| u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
   case(3)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      call vecresMPI(two,ndime*npoinLocal,a%veloc(1,1,1),a%veloc(1,1,3),rivel,zelmn,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npoinLocal,      a%press(1,1),  a%press(1,3),  ripre,zelmn,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npoinLocal,      a%tempe(1,1),  a%tempe(1,3),  ritem,zelmn,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      rilmn = max(rivel,ritem,ripre)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,rivel,ritem,ripre 
         if(rilmn.le.a%sstol) then
            a%kfl_stead = 1
            write(a%lun_conve,102) a%istep
         end if
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   end select       

   !Formats
100 format('$ ','  Time','      Inner',&
       &      '       Current','      Velocity','    Temperature','     Pressure',/,&
       & '$ ','  step','  iteration',&
       &      '          time','      residual','      residual','      residual')
103 format('$ ','  Time','     Global','      Inner',' Non-Linear',&
       &      '       Current','      Velocity','    Temperature','     Pressure',/,&
       & '$ ','  step','  iteration','  iteration',&
       &      '          time','      residual','      residual','      residual')

101 format(i9,1x,i9,5(2x,e12.6),1(i9))
104 format(2x,4(2x,i9),3(2x,e12.6),4x,a)
102 format('$ >>>  VELOCITY IS STATIONARY AT TIME STEP ',i5)

200 format(///,&
       5x, '>>>  EVOLUTION OF NUMERICAL PARAMETERS:',///,&
       10x,a)
201 format('$ ','       Time','      Inner',&
       &      '       Current','    Tau_min       Tau_max',&
       &      '       Tau_mean       Re_min        Re_max',&
       &      '        Re_mean       yp_min        yp_max',&
       &      '        yp_mean',/,&
       & '$ ','       step','  iteration',&
       &      '          time')
202 format(4x,i9,2x,i9,7(2x,e12.6))
300 format('$ ','   Step','   Iteration','     SGS iterations',&
       &        '    Non Linear Subscales Maximum Residual',/,&
       &        45x,  ' Velocity','     Temperature')


end subroutine lmn_cvgunk

subroutine lmn_InnerResiduals(a,rvlmn,rplmn,rtlmn)
   use typre
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   real(rp) :: rvlmn,rplmn,rtlmn
   
   integer(ip) :: ndime, npoinLocal
   
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   call vecresMPIHeterogeneous(2_ip,ndime+2,ndime,npoinLocal,1,1,ndime,a%unkno,a%veloc(1,1,1),rvlmn,zelmn,a%MPIcomm,a%MPIroot)
   call vecresMPIHeterogeneous(2_ip,ndime+2,1,npoinLocal,ndime+2,1,1,a%unkno,a%press(1,1),rplmn,zelmn,a%MPIcomm,a%MPIroot)     
   call vecresMPIHeterogeneous(2_ip,ndime+2,1,npoinLocal,ndime+1,1,1,a%unkno,a%tempe(1,1),rtlmn,zelmn,a%MPIcomm,a%MPIroot)     
   
end subroutine
   
   


