subroutine nsc_exaerr(a)    
!-----------------------------------------------------------------------
!****f* Nscomp/nsc_exaerr
! NAME 
!    nsc_norms
! DESCRIPTION
!    This routine computes the FEM errors.
!    refered to the exact analytical solution nsc_exacso
!    The norms are normalized by the appropriate norm of the exact sol. 
! USES
! USED BY
!    nsc_output
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_NSCompressible
   use Mod_Mesh
   use Mod_Element
   use Mod_memor
   use MPI
   use Mod_NscExacso
   
   implicit none
   class(NSCompressibleProblem) :: a   

   type(MemoryMan) :: Memor
   class(FiniteElement), pointer :: e => NULL() 
   type(NscExacso) :: exacso 

   !Elemental Values
   real(rp), allocatable :: elden(:), elmom(:,:), elene(:)
   real(rp), allocatable :: elpre(:), elvel(:,:), eltem(:)

   !Numerical Gauss Values
   real(rp), allocatable :: gppre(:), gpvel(:), gptem(:)
   real(rp), allocatable :: gpcod(:),gradp(:),gradv(:,:),gradt(:) 

   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:),exprg(:),exteg(:)
   real(rp)                :: expre,extem
   real(rp)    :: dvolu
   
   !Error Values
   real(rp)    :: diffp,abpre,difeu,abvel,difft,abtem,zero_rp,weightfactor

   !Errors 
   real(rp)    :: erp01(2),erp02(2),erp0i(2),erp11(2),erp12(2),erp1i(2)
   real(rp)    :: eru01(2),eru02(2),eru0i(2),eru11(2),eru12(2),eru1i(2) 
   real(rp)    :: ert01(2),ert02(2),ert0i(2),ert11(2),ert12(2),ert1i(2)   
   
   !Reduced Errors 
   real(rp)    :: erp01r(2),erp02r(2),erp0ir(2),erp11r(2),erp12r(2),erp1ir(2)
   real(rp)    :: eru01r(2),eru02r(2),eru0ir(2),eru11r(2),eru12r(2),eru1ir(2)
   real(rp)    :: ert01r(2),ert02r(2),ert0ir(2),ert11r(2),ert12r(2),ert1ir(2)   

   integer :: ierr   
   integer(ip) :: nelem,npoin,nelty,ndime,idime,jdime,ielem,igaus
   integer(ip)              :: icount,inode,npoinLocal
   integer(ip), pointer     :: nnode(:) => NULL(),ngaus(:) => NULL(),ltopo(:) => NULL() 
   real(rp)    :: acvis,actco,accph,accvh
   

! Initializations

   call a%SpecificNSCompExaerr

   erp01=0.0_rp
   erp02=0.0_rp
   erp0i=0.0_rp
   erp11=0.0_rp
   erp12=0.0_rp
   erp1i=0.0_rp
   eru01=0.0_rp
   eru02=0.0_rp
   eru0i=0.0_rp
   eru11=0.0_rp
   eru12=0.0_rp
   eru1i=0.0_rp
   ert01=0.0_rp
   ert02=0.0_rp
   ert0i=0.0_rp
   ert11=0.0_rp
   ert12=0.0_rp
   ert1i=0.0_rp

   erp01r=0.0_rp
   erp02r=0.0_rp
   erp0ir=0.0_rp
   erp11r=0.0_rp
   erp12r=0.0_rp
   erp1ir=0.0_rp
   eru01r=0.0_rp
   eru02r=0.0_rp
   eru0ir=0.0_rp
   eru11r=0.0_rp
   eru12r=0.0_rp
   eru1ir=0.0_rp
   ert01r=0.0_rp
   ert02r=0.0_rp
   ert0ir=0.0_rp
   ert11r=0.0_rp
   ert12r=0.0_rp
   ert1ir=0.0_rp
   zero_rp=0.0_rp

   call a%Mesh%GetNdime(ndime)   
   call a%Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNpoinLocal(npoinLocal)
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)

   ! Element object Initialization
   call a%Mesh%ElementAlloc(e,Memor,'DefaultRule','nsc_exaerr')   
   
   ! Allocate elemental values
   call a%Memor%alloc(e%mnode,elden,'elden','nsc_exaerr')     
   call a%Memor%alloc(ndime,e%mnode,elmom,'elmom','nsc_exaerr')     
   call a%Memor%alloc(e%mnode,elene,'elene','nsc_exaerr')     
   call a%Memor%alloc(e%mnode,elpre,'elpre','nsc_exaerr')     
   call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','nsc_exaerr')     
   call a%Memor%alloc(e%mnode,eltem,'eltem','nsc_exaerr')     

   ! Allocate gauss point values
   call a%Memor%alloc(ndime,gpcod,'gpcod','nsc_exaerr')
   call a%Memor%alloc(1,gppre,'gppre','nsc_exaerr')     
   call a%Memor%alloc(ndime,gpvel,'gpvel','nsc_exaerr')     
   call a%Memor%alloc(1,gptem,'gptem','nsc_exaerr')     
   call a%Memor%alloc(ndime,gradp,'gradp','nsc_exaerr')     
   call a%Memor%alloc(ndime,ndime,gradv,'gradv','nsc_exaerr')     
   call a%Memor%alloc(ndime,gradt,'gradt','nsc_exaerr')     

   ! Allocate Exact Values
   call a%Memor%alloc(ndime,exvel,'exvel','nsc_exaerr')     
   call a%Memor%alloc(ndime,exprg,'exprg','nsc_exaerr')     
   call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsc_exaerr')     
   call a%Memor%alloc(ndime,exteg,'exteg','nsc_exaerr')     

   elements: do ielem=1,nelem  
      !Load Element   
      call a%Mesh%ElementLoad(ielem,e)

      icount = 0
      do inode = 1,e%pnode
         if (e%lnods(inode) <= npoinLocal) then
            icount = icount +1
         endif
      enddo

      weightfactor = real(icount)/real(e%pnode)
   
      call e%elmdcg
      
      dvolu = 0.0_rp

      ! Gather operations    
      call e%gather(1,elpre,a%press(:,1))
      call e%gather(ndime,elvel,a%veloc(:,:,1))
      call e%gather(1,eltem,a%tempe(:,1))

      gauss_points: do igaus=1,e%pgaus
         e%igaus = igaus      

         dvolu = e%weigp(e%igaus)*e%detjm
         
         !Pressure and Velocity gradient
         call e%gradient(1,elpre,gradp)          
         call e%gradient(ndime,elvel,gradv) 
         call e%gradient(1,eltem,gradt)          

         !Gauss point values
         call e%interpg(ndime,e%elcod,gpcod)
         call e%interpg(1,elpre,gppre) 
         call e%interpg(ndime,elvel,gpvel)       
         call e%interpg(1,eltem,gptem) 
                                            
         ! Exact solution
         call exacso%nsc_ComputeSolution(ndime,gpcod,a)
         call exacso%nsc_GetPressure(ndime,expre,exprg)         
         call exacso%nsc_GetVelocity(ndime,exvel,exveg) 
         call exacso%nsc_GetTemperature(ndime,extem,exteg)         

         !Errors
         diffp = abs(gppre(1)-expre)
!         a%coResGP(ielem)%a(igaus)=diffp
         abpre = abs(expre)
         erp01(1) = erp01(1) + weightfactor*diffp*dvolu
         erp02(1) = erp02(1) + weightfactor*diffp*diffp*dvolu
         erp0i(1) = max(erp0i(1),diffp)
         erp01(2) = erp01(2) + weightfactor*abpre*dvolu
         erp02(2) = erp02(2) + weightfactor*abpre*abpre*dvolu       
         erp0i(2) = max(erp0i(2),expre)
         difft = abs(gptem(1)-extem)
!         a%enResGP(ielem)%a(igaus)=difft
         abtem = abs(extem)
         ert01(1) = ert01(1) + weightfactor*difft*dvolu
         ert02(1) = ert02(1) + weightfactor*difft*difft*dvolu
         ert0i(1) = max(ert0i(1),difft)
         ert01(2) = ert01(2) + weightfactor*abtem*dvolu
         ert02(2) = ert02(2) + weightfactor*abtem*abtem*dvolu       
         ert0i(2) = max(ert0i(2),extem)
         do idime=1,ndime
               difeu = abs(gpvel(idime)-exvel(idime))
!               a%moResGP(ielem)%a(idime,igaus)=difeu
               abvel = abs(exvel(idime))
               eru01(1) = eru01(1) + weightfactor*difeu*dvolu
               eru02(1) = eru02(1) + weightfactor*difeu*difeu*dvolu
               eru0i(1) = max(eru0i(1),difeu)
               eru01(2) = eru01(2) + weightfactor*abvel*dvolu
               eru02(2) = eru02(2) + weightfactor*abvel*abvel*dvolu
               eru0i(2) = max(eru0i(2),abvel)
               diffp = abs(gradp(idime)-exprg(idime))
               abpre = abs(exprg(idime)) 
               erp11(1) = erp11(1) + weightfactor*diffp*dvolu
               erp12(1) = erp12(1) + weightfactor*diffp*diffp*dvolu
               erp1i(1) = max(erp1i(1),diffp)
               erp11(2) = erp11(2) + weightfactor*abpre*dvolu
               erp12(2) = erp12(2) + weightfactor*abpre*abpre*dvolu
               erp1i(2) = max(erp1i(2),abpre)
               difft = abs(gradt(idime)-exteg(idime))
               abtem = abs(exteg(idime)) 
               ert11(1) = ert11(1) + weightfactor*difft*dvolu
               ert12(1) = ert12(1) + weightfactor*difft*difft*dvolu
               ert1i(1) = max(ert1i(1),difft)
               ert11(2) = ert11(2) + weightfactor*abtem*dvolu
               ert12(2) = ert12(2) + weightfactor*abtem*abtem*dvolu
               ert1i(2) = max(ert1i(2),abtem)
               do jdime=1,ndime
                  difeu = abs(gradv(idime,jdime)-exveg(idime,jdime))
                  abvel = abs(exveg(idime,jdime))
                  eru11(1) = eru11(1) + weightfactor*difeu*dvolu
                  eru12(1) = eru12(1) + weightfactor*difeu*difeu*dvolu
                  eru1i(1) = max(eru1i(1),difeu)
                  eru11(2) = eru11(2) + weightfactor*abvel*dvolu
                  eru12(2) = eru12(2) + weightfactor*abvel*abvel*dvolu
                  eru1i(2) = max(eru1i(2),abvel)
               end do
         end do
      end do gauss_points
   
   end do elements

   call  MPI_REDUCE( erp01,  erp01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( erp02,  erp02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( erp0i,  erp0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
   
   call  MPI_REDUCE( ert01,  ert01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( ert02,  ert02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( ert0i,  ert0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

   call  MPI_REDUCE( eru01,  eru01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( eru02,  eru02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( eru0i,  eru0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

   call  MPI_REDUCE( erp11,  erp11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( erp12,  erp12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( erp1i,  erp1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

   call  MPI_REDUCE( ert11,  ert11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( ert12,  ert12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( ert1i,  ert1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

   call  MPI_REDUCE( eru11,  eru11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( eru12,  eru12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   call  MPI_REDUCE( eru1i,  eru1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )       

   erp02r(1) = sqrt(erp02r(1))
   erp12r(1) = sqrt(erp12r(1))
   eru02r(1) = sqrt(eru02r(1))
   eru12r(1) = sqrt(eru12r(1))
   ert02r(1) = sqrt(ert02r(1))
   ert12r(1) = sqrt(ert12r(1))
   
   if(erp01r(2).gt.zero_rp) erp01r(1) = erp01r(1)/erp01r(2) 
   if(erp02r(2).gt.zero_rp) erp02r(1) = erp02r(1)/sqrt(erp02r(2))
   if(erp0ir(2).gt.zero_rp) erp0ir(1) = erp0ir(1)/erp0ir(2)

   if(eru01r(2).gt.zero_rp) eru01r(1) = eru01r(1)/eru01r(2) 
   if(eru02r(2).gt.zero_rp) eru02r(1) = eru02r(1)/sqrt(eru02r(2))
   if(eru0ir(2).gt.zero_rp) eru0ir(1) = eru0ir(1)/eru0ir(2) 

   if(ert01r(2).gt.zero_rp) ert01r(1) = ert01r(1)/ert01r(2) 
   if(ert02r(2).gt.zero_rp) ert02r(1) = ert02r(1)/sqrt(ert02r(2))
   if(ert0ir(2).gt.zero_rp) ert0ir(1) = ert0ir(1)/ert0ir(2) 

   if(erp11r(2).gt.zero_rp) erp11r(1) = erp11r(1)/erp11r(2) 
   if(erp12r(2).gt.zero_rp) erp12r(1) = erp12r(1)/sqrt(erp12r(2))
   if(erp1ir(2).gt.zero_rp) erp1ir(1) = erp1ir(1)/erp1ir(2)

   if(eru11r(2).gt.zero_rp) eru11r(1) = eru11r(1)/eru11r(2) 
   if(eru12r(2).gt.zero_rp) eru12r(1) = eru12r(1)/sqrt(eru12r(2))
   if(eru1ir(2).gt.zero_rp) eru1ir(1) = eru1ir(1)/eru1ir(2) 

   if(ert11r(2).gt.zero_rp) ert11r(1) = ert11r(1)/ert11r(2) 
   if(ert12r(2).gt.zero_rp) ert12r(1) = ert12r(1)/sqrt(ert12r(2))
   if(ert1ir(2).gt.zero_rp) ert1ir(1) = ert1ir(1)/ert1ir(2) 

   if (a%MPIrank == a%MPIroot) then 

        write(a%lun_outph,100)  &
        erp01r(1),erp02r(1),erp0ir(1),erp11r(1),erp12r(1),erp1ir(1), &
        eru01r(1),eru02r(1),eru0ir(1),eru11r(1),eru12r(1),eru1ir(1), &
        ert01r(1),ert02r(1),ert0ir(1),ert11r(1),ert12r(1),ert1ir(1)
      if (a%kfl_flush == 1) call flush(a%lun_outph)
   endif
      
   ! Deallocate gauss point values
   call a%Memor%dealloc(ndime,gpcod,'gpcod','nsc_exaerr')
   call a%Memor%dealloc(1,gppre,'gppre','nsc_exaerr')     
   call a%Memor%dealloc(ndime,gpvel,'gpvel','nsc_exaerr')     
   call a%Memor%dealloc(1,gptem,'gptem','nsc_exaerr')     
   call a%Memor%dealloc(ndime,gradp,'gradp','nsc_exaerr')     
   call a%Memor%dealloc(ndime,ndime,gradv,'gradv','nsc_exaerr')     
   call a%Memor%dealloc(ndime,gradt,'gradt','nsc_exaerr')     

   ! Deallocate Exact Values
   call a%Memor%dealloc(ndime,exvel,'exvel','nsc_exaerr')     
   call a%Memor%dealloc(ndime,exprg,'exprg','nsc_exaerr')     
   call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsc_exaerr')     
   call a%Memor%dealloc(ndime,exteg,'exteg','nsc_exaerr')     

   ! Deallocate elemental values
   call a%Memor%dealloc(e%mnode,elden,'elden','nsc_exaerr')     
   call a%Memor%dealloc(ndime,e%mnode,elmom,'elmom','nsc_exaerr')     
   call a%Memor%dealloc(e%mnode,elene,'elene','nsc_exaerr')     
   call a%Memor%dealloc(e%mnode,elpre,'elpre','nsc_exaerr')     
   call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','nsc_exaerr')     
   call a%Memor%dealloc(e%mnode,eltem,'eltem','nsc_exaerr')     

   call a%Mesh%ElementDealloc(e,Memor,'DefaultRule','nsc_exaerr')   

  100 format(///,10X,'FINITE ELEMENT ERRORS',                              &
&              /,10X,'=====================',//,                           &
&  '          W(0,1)P ',E12.5,/, &
&  '          W(0,2)P ',E12.5,/, &
&  '          W(0,i)P ',E12.5,/, &
&  '          W(1,1)P ',E12.5,/, &
&  '          W(1,2)P ',E12.5,/, &
&  '          W(1,i)P ',E12.5,/, &
&  '          W(0,1)V ',E12.5,/, &
&  '          W(0,2)V ',E12.5,/, &
&  '          W(0,i)V ',E12.5,/, &
&  '          W(1,1)V ',E12.5,/, &
&  '          W(1,2)V ',E12.5,/, &
&  '          W(1,i)V ',E12.5,/, &
&  '          W(0,1)T ',E12.5,/, &
&  '          W(0,2)T ',E12.5,/, &
&  '          W(0,i)T ',E12.5,/, &
&  '          W(1,1)T ',E12.5,/, &
&  '          W(1,2)T ',E12.5,/, &
&  '          W(1,i)T ',E12.5)

end subroutine nsc_exaerr
