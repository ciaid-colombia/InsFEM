subroutine supf_cvgunkS(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    incoma%pressible NS equations
   !-----------------------------------------------------------------------
   use typre
   use Mod_SUPFractionalStep
   use MPI
   use def_parame
   implicit none
   class(SUPFractionalStepProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip), save       :: ipass=0
   integer(ip)             :: npoinLocal,ndime,idime,ipoin
   real(rp)                :: va,vo,numer,denom,vnume,vdeno,pnume,pdeno,prnsi,rinsi,ripre, &
                                srnsi,risig
   real(rp)                :: rwa(4),rwa2(4)
   integer                 :: ierr,auxtens
   real(rp)                :: zensi
   
   interface

      subroutine supf_InnerResidualsS(a,rinsi,rprnsi,rsnsi)
         use typre      
         import SUPFractionalStepProblem
         implicit none
         class(SUPFractionalStepProblem) :: a
         real(rp) :: rinsi,rprnsi,rsnsi           
      end subroutine    
   
   end interface  
   
   
   zensi=0.0_rp
   
   select case(itask)
      
   !Check convergence of the inner iterations, always in the L2 norm:
   !|| u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
   case(1)
      
     call supf_InnerResidualsS(a,rinsi,prnsi,srnsi)
         
      if (a%MPIrank == a%MPIroot) then
         if((srnsi<a%cotol).or.(a%iteraS>=a%maxit)) a%kfl_goiteS = 0        
         
         if(ipass==0) then
            ipass=1
            write(a%lun_conve,100)
            write(a%lun_outph,200)
            write(a%lun_outph,201)
            !write(a%lun_nolin,300)
         end if
         write(a%lun_conve,101) a%istep,a%iteraS,a%ctime,rinsi,prnsi,srnsi
         write(a%lun_outph,202) a%istep,a%iteraS,a%ctime,& 
               a%tamin,a%tamax,a%tamea,&
               a%remin,a%remax,a%remea,&
               a%ypmin,a%ypmax,a%ypmean
         ! Write in the log file.
         if (a%kfl_flush == 1) then
            call flush(a%lun_outph)
            call flush(a%lun_conve)
         endif
      endif
      CALL MPI_BCAST(a%kfl_goiteS, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
   case(2)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      auxtens=(ndime-1)*(ndime-1)+2      
      call vecresMPI(a%kfl_normc,ndime*npoinLocal,a%veloc(1,1,1),a%veloc(1,1,2),a%resid,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(a%kfl_normc,npoinLocal,      a%press(1,1),  a%press(1,2),  a%resip,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(a%kfl_normc,auxtens*npoinLocal,a%sigma(1,1,1),a%sigma(1,1,2),a%resis,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)      
      
   
   !Check residual of the time evolution, always in the L2 norm:
   !|| u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
   case(3)
      call a%Mesh%GetNpoinLocal(npoinLocal)
      call a%Mesh%GetNdime(ndime)
      auxtens=(ndime-1)*(ndime-1)+2      
      call vecresMPI(two,ndime*npoinLocal,a%veloc(1,1,1),a%veloc(1,1,3),rinsi,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npoinLocal,      a%press(1,1),  a%press(1,3),  ripre,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,auxtens*npoinLocal,a%sigma(1,1,1),a%sigma(1,1,3),risig,zensi,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,rinsi,ripre,risig  
         if((rinsi.le.a%sstol)) then         
!         if((rinsi.le.a%sstol).and.(ripre.le.a%sstol).and.(risig.le.a%sstol)) then
            a%kfl_stead = 1
            write(a%lun_conve,102) a%istep
         end if
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
      
   end select       

   !Formats
100 format('$ ','       Time','      Inner',&
       &      '       Current','      Velocity','      Pressure',/,&
       & '$ ','       step','  iteration',&
       &      '          time','      residual','      residual')
103 format('$ ','       Time','     Global','      Inner',' Non-Linear',&
       &      '       Current','      Velocity','      Pressure',/,&
       & '$ ','       step','  iteration','  iteration','  iteration',&
       &      '          time','      residual','      residual')

101 format(4x,i9,2x,i9,20(2x,e12.6))
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
202 format(4x,i9,2x,i9,20(2x,e12.6))

300 format('$ ','       Time','     Global','      Inner',&
       &      '       Current','      Max Vel1','      Arg Max1',/,&
       & '$ ','       step','  iteration','  iteration',&
       &      '          time','      residual','      residual')
301 format(4x,i9,2x,i9,2x,i9,2x,e12.6,20(2x,e12.6,2x,i6))

end subroutine supf_cvgunkS

subroutine supf_InnerResidualsS(a,rinsi,rprnsi,rsnsi)
   use typre
   use Mod_SUPFractionalStep
   implicit none
   class(SUPFractionalStepProblem) :: a
   real(rp) :: rinsi,rprnsi,rsnsi,zensi
   
   integer(ip) :: ndime,npoinLocal,auxtens
  
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%Mesh%GetNdime(ndime)
   zensi=0.0_rp
   auxtens=(ndime-1)*(ndime-1)+2  
   
   call vecresMPIHeterogeneous(2_ip,auxtens,auxtens,npoinLocal,1,1,auxtens,a%unknoS,a%sigma(1,1,1),rsnsi,zensi,a%MPIcomm,a%MPIroot) 
   rinsi  = 0.0_rp
   rprnsi = 0.0_rp
   
end subroutine