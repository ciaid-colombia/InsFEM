subroutine sldup_cvgunk(a,itask)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   !    This routine performs several convergence checks for the 
   !    Cauchy equations
   !-----------------------------------------------------------------------
   use MPI
   use typre
   use def_parame
   use Mod_UPSolids
   implicit none
   class(UPSolidsProblem) :: a
   integer(ip), intent(in)    :: itask
   
   integer(ip), save       :: ipass=0, ipasscp=0
   integer(ip)             :: npl,nd,idime,ipoin
   real(rp)                :: ursld,prsld,srsld
   integer                 :: ierr
   logical                 :: isALE

   call a%Mesh%Getndime(nd)
   call a%Mesh%GetnpoinLocal(npl)

   select case(itask)

   !Check convergence of the inner iterations, always in the L2 norm:
   !|| u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
   case(1)
      
       call a%InnerResiduals(srsld,ursld,prsld)

       if (a%MPIrank == a%MPIroot) then
         if(((ursld<a%cotol).and.(prsld<a%cotol)).or.(a%itera>=a%maxit)) a%kfl_goite = 0      

         if(ipass==0) then
             ipass=1
             write(a%lun_conve,100)
             write(a%lun_outph,200)
         end if

         write(a%lun_conve,101) a%istep,a%itera,a%ctime,ursld,prsld
         
         ! Write in the log file.
         if (a%kfl_flush == 1) then
            call flush(a%lun_outph)
            call flush(a%lun_conve)
         endif
      endif
      CALL MPI_BCAST(a%kfl_goite, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)
   
   !Check convergence of the outer iterations in the norm selected by the user:
   !|| u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
   case(2)

      call vecresMPI(a%kfl_normc,nd*npl,a%disp(:,:,1),a%disp(:,:,2),ursld,zesup,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npl        ,a%press(:,1)  ,a%press(:,2)  ,prsld,zesup,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

   !Check residual of the time evolution, always in the L2 norm:
   !|| u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
   case(3)

      call vecresMPI(two,nd*npl  ,a%disp(:,:,1) ,a%disp(:,:,3) ,ursld,zesup,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
      call vecresMPI(two,npl        ,a%press(:,1)  ,a%press(:,3)  ,prsld,zesup,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

      if (a%MPIrank == a%MPIroot) then
         write(a%lun_conve,101) a%istep,0,a%ctime,ursld,prsld
         if((ursld.le.a%sstol).and.(prsld.le.a%sstol)) then
            a%kfl_stead = 1
            write(a%lun_conve,102) a%istep
         end if
         call a%Mesh%GetALE(isALE)
         if (isALE) a%kfl_stead=0
      endif
      CALL MPI_BCAST(a%kfl_stead, 1, MPI_INTEGER4, a%MPIroot, a%MPIcomm, ierr)

      !Check residual of case coupling, for example FSI
   case(4)

       call vecresMPI(two,nd*npl,a%disp(:,:,1) ,a%disp_cp(:,:) ,ursld,zesup,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)
       call vecresMPI(two,npl   ,a%press(:,1)  ,a%press_cp(:)  ,prsld,zesup,a%kfl_MPIComType,a%MPIcomm,a%MPIroot,a%MPIrank,a%MPIsize)

       a%cpres = max(ursld,prsld)

      if (a%MPIrank == a%MPIroot) then

          if(ipasscp==0) then
              ipasscp=1
              write(a%lun_cpconve,100)
              write(a%lun_outph,200)
          end if

          write(a%lun_cpconve,101) a%istep,a%cpiter,a%ctime,ursld,prsld

          if (a%kfl_flush == 1) then
              call flush(a%lun_cpconve)
          endif     

          if(a%cpres.le.a%cptol) then
              a%kfl_coupconv= .true.
          end if
      endif

      CALL MPI_BCAST(a%kfl_coupconv, 1, MPI_LOGICAL, a%MPIroot, a%MPIcomm, ierr)
      
      
   end select       

   !Formats
100 format('$ ','       Time','      Inner','       Current','         Displ','      Pressure',/,&
         & '$ ','       step','  iteration','          time','      residual','      residual')

101 format(4x,i9,2x,i9,20(2x,e12.6),20(2x,e12.6),20(2x,e12.6))
102 format('$ >>>  DISPLACEMENTS STATIONARY AT TIME STEP ',i5)

200 format(///,&
       5x, '>>>  EVOLUTION OF NUMERICAL PARAMETERS:',///,&
       10x,a)

end subroutine sldup_cvgunk

subroutine sldup_InnerResiduals(a,srsld,ursld,prsld)
   use typre
   use Mod_SUPSolids
   implicit none
   class(SUPSolidsProblem) :: a
   real(rp)    :: ursld,srsld,prsld
   integer(ip) :: nd,npl,np,idime,ipoin
   integer(ip) :: u1,uf,s1,sf,p1,bc
   real(rp), allocatable   :: dis_aux(:,:)
   real(rp), allocatable   :: pre_aux(:)
   
   call a%Mesh%GetnpoinLocal(npl)
   call a%Mesh%Getnpoin(np)
   call a%Mesh%Getndime(nd)
   srsld = 0.0_rp

   call a%GetMatrixOrganization(u1,uf,s1,sf,p1,bc)


   call a%Memor%alloc(nd,npl,dis_aux,'dis','sldsup_aux')
   call a%Memor%alloc(   npl,pre_aux,'pre','sldsup_aux')

   dis_aux = a%unkno(u1:uf,:)!+a%disp(:,:,1)
   pre_aux = a%unkno(p1,:)   !+a%press(:,1)

   call vecresMPIHeterogeneous(2_ip,nd,nd,npl,1,1,nd,dis_aux,a%disp(:,:,1) ,ursld,zesup,a%MPIcomm,a%MPIroot)
   call vecresMPIHeterogeneous(2_ip,1 ,1 ,npl,1,1,1 ,pre_aux,a%press(:,1)  ,prsld,zesup,a%MPIcomm,a%MPIroot)

   call a%Memor%dealloc(nd,npl,dis_aux,'dis','sldsup_aux')
   call a%Memor%dealloc(   npl,pre_aux,'pre','sldsup_aux')
end subroutine
