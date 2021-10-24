subroutine plcd_exaerr(a)    
!-----------------------------------------------------------------------
!****f* Nstinc/plcd_exaerr
! NAME 
!    plcd_exaerr
! DESCRIPTION
!    This routine computes the FEM errors (referred to an analytical
!    solution defined in exacso.f90). The errors are normalized by the 
!    appropriate norm of the exact solution, except when this norm is zero.      
! USES
! USED BY
!    plcd_output
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_NavierStokes   
   use Mod_Mesh
   use Mod_Element
   use Mod_memor
   use MPI
   use Mod_NsiExacso 
   
   implicit none
   class(NavierStokesProblem) :: a   
   type(MemoryMan) :: Memor
   type(FemMesh)    :: Mesh  
   class(FiniteElement), pointer :: e => NULL() 
   type(NsiExacso)  :: exacso      

   !Numerical Gaus Values
   real(rp), allocatable   :: vegau(:),grave(:,:),gradp(:),gpcod(:)
   real(rp)    :: prgau(1),dvolu
   
   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:),exprg(:)
   real(rp)                :: expre
   
   real(rp)    :: difeu,abvel,diffp,abpre,dummr,zero_rp,weightfactor
   !Errors 
   real(rp)    :: erp01(2),erp02(2),erp0i(2),erp11(2),erp12(2),erp1i(2)
   real(rp)    :: eru01(2),eru02(2),eru0i(2),eru11(2),eru12(2),eru1i(2) 
   !Reduced Errors 
   real(rp)    :: erp01r(2),erp02r(2),erp0ir(2),erp11r(2),erp12r(2),erp1ir(2)
   real(rp)    :: eru01r(2),eru02r(2),eru0ir(2),eru11r(2),eru12r(2),eru1ir(2)    

   integer(ip) :: ievab     ! Indices and dimensions
   integer(ip) :: pevab,pevat
   integer(ip) :: pelty,nelem,npoin,nelty,ndime,dummi,idime,jdime,ielem,igaus
   integer(ip), pointer     :: nnode(:),ngaus(:),ltopo(:) 

   
   !Elemental Values
   real(rp), allocatable :: elveloc(:,:),elpress(:) 
   integer :: ierr   
   integer(ip)              :: i,root,icount,inode,npoinLocal
   
!
! Initializations

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
   zero_rp=0.0_rp

   call a%Mesh%GetNdime(ndime)   
   call a%Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNpoinLocal(npoinLocal)
 ! Element object Initialization
   call a%Mesh%ElementAlloc(e,Memor,'DefaultRule','plcd_exaerr')   
   

   ! Allocate exact components
   call a%Memor%alloc(ndime,gpcod,'gpcod','plcd_exaerr') 
   call a%Memor%alloc(ndime,exvel,'exvel','plcd_exaerr')
   call a%Memor%alloc(ndime,vegau,'vegau','plcd_exaerr')
   call a%Memor%alloc(ndime,exprg,'exprg','plcd_exaerr')  
   call a%Memor%alloc(ndime,gradp,'gradp','plcd_exaerr')
   call a%Memor%alloc(ndime,ndime,exveg,'exveg','plcd_exaerr')
   call a%Memor%alloc(ndime,ndime,grave,'grave','plcd_exaerr')
   
   !Allocate exact Values
   call a%Memor%alloc(ndime,e%mnode,elveloc,'elveloc','plcd_exaerr')     
   call a%Memor%alloc(e%mnode,elpress,'elpress','plcd_exaerr')     
    
   
   
!   ielty=1
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
   
      
      ! Gather operations    
      call e%gather(ndime,elveloc,a%veloc)
      call e%gather(1_ip,elpress,a%press)

      call e%elmdcg

      gauss_points: do igaus=1,e%pgaus
         e%igaus = igaus      

         dvolu = e%weigp(e%igaus)*e%detjm
         

         ! Velocity, pressure, velocity gradient and pressure gradient
         
         !Pressure and Velocity gradient
         call e%gradient(1,elpress,gradp)          
         call e%gradient(ndime,elveloc,grave)           

         !Gauss point values
         call e%interpg(ndime,elveloc,vegau)       
         call e%interpg(ndime,e%elcod,gpcod)
         call e%interpg(1,elpress,prgau) 
         
         ! Exact solution
         call exacso%plcd_ComputeSolution(ndime,gpcod,a)
         call exacso%plcd_GetPressure(ndime,expre,exprg)         
         call exacso%plcd_GetVelocity(ndime,exvel,exveg)      

                                            
        ! Errors
         diffp = abs(prgau(1)-expre)
         abpre = abs(expre)
         erp01(1) = erp01(1) + weightfactor*diffp*dvolu
         erp02(1) = erp02(1) + weightfactor*diffp*diffp*dvolu
         erp0i(1) = max(erp0i(1),diffp)
         erp01(2) = erp01(2) + weightfactor*abpre*dvolu
         erp02(2) = erp02(2) + weightfactor*expre*expre*dvolu       
         erp0i(2) = max(erp0i(2),expre)
         do idime=1,ndime
            difeu = abs(vegau(idime)-exvel(idime))
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
            do jdime=1,ndime
               difeu = abs(grave(idime,jdime)-exveg(idime,jdime))
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

   ! deallocate elemental values
   call a%Memor%dealloc(ndime,e%mnode,elveloc,'elveloc','plcd_exaerr')     
   call a%Memor%dealloc(e%mnode,elpress,'elpress','plcd_exaerr')     
   
   ! deallocate exact components
   call a%Memor%dealloc(ndime,gpcod,'gpcod','plcd_exaerr') 
   call a%Memor%dealloc(ndime,exvel,'exvel','plcd_exaerr')
   call a%Memor%dealloc(ndime,vegau,'vegau','plcd_exaerr')
   call a%Memor%dealloc(ndime,exprg,'exprg','plcd_exaerr')  
   call a%Memor%dealloc(ndime,gradp,'gradp','plcd_exaerr')
   call a%Memor%dealloc(ndime,ndime,exveg,'exveg','plcd_exaerr')
   call a%Memor%dealloc(ndime,ndime,grave,'grave','plcd_exaerr')     
   
   
   call a%Mesh%ElementDealloc(e,Memor,'DefaultRule','plcd_exaerr')   
     
      call  MPI_REDUCE( erp01,  erp01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp02,  erp02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp0i,  erp0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
       
      call  MPI_REDUCE( eru01,  eru01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru02,  eru02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru0i,  eru0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
  
      call  MPI_REDUCE( erp11,  erp11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp12,  erp12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp1i,  erp1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
 
      call  MPI_REDUCE( eru11,  eru11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru12,  eru12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru1i,  eru1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )       
      
      erp02r(1) = sqrt(erp02r(1))
      erp12r(1) = sqrt(erp12r(1))
      eru02r(1) = sqrt(eru02r(1))
      eru12r(1) = sqrt(eru12r(1))
      
      if(erp01r(2).gt.zero_rp) erp01r(1) = erp01r(1)/erp01r(2) 
      if(erp02r(2).gt.zero_rp) erp02r(1) = erp02r(1)/sqrt(erp02r(2))
      if(erp0ir(2).gt.zero_rp) erp0ir(1) = erp0ir(1)/erp0ir(2)
      if(eru01r(2).gt.zero_rp) eru01r(1) = eru01r(1)/eru01r(2) 
      if(eru02r(2).gt.zero_rp) eru02r(1) = eru02r(1)/sqrt(eru02r(2))
      if(eru0ir(2).gt.zero_rp) eru0ir(1) = eru0ir(1)/eru0ir(2) 
      if(erp11r(2).gt.zero_rp) erp11r(1) = erp11r(1)/erp11r(2) 
      if(erp12r(2).gt.zero_rp) erp12r(1) = erp12r(1)/sqrt(erp12r(2))
      if(erp1ir(2).gt.zero_rp) erp1ir(1) = erp1ir(1)/erp1ir(2)
      if(eru11r(2).gt.zero_rp) eru11r(1) = eru11r(1)/eru11r(2) 
      if(eru12r(2).gt.zero_rp) eru12r(1) = eru12r(1)/sqrt(eru12r(2))
      if(eru1ir(2).gt.zero_rp) eru1ir(1) = eru1ir(1)/eru1ir(2) 



   if (a%MPIrank == a%MPIroot) then 
   
!        write(lun_outpu_nsi,100)  &
!        eru01r(1),erp01r(1),eru02r(1),erp02r(1),eru0ir(1),erp0ir(1),  &
!        eru11r(1),erp11r(1),eru12r(1),erp12r(1),eru1ir(1),erp1ir(1)  

      write(*,*) 'diferencia u   diferencia p'
      write(*,*)  eru01r(1), erp01r(1)
      write(*,*) 'L2 u     L2 p'   
      write(*,*) eru02r(1), erp02r(1)
      write(*,*) 'Max(abs_exvel)  Max(abs_exp)' 
      write(*,*) eru0ir(1), erp0ir(1)
      write(*,*) 'dif_gradu    dif_gradu'
      write(*,*) eru11r(1), erp11r(1)
      write(*,*) 'L2 gradu     L2 gradp'
      write(*,*) eru12r(1), erp12r(1)
      write(*,*) 'Max(abs_exgradu)  Max(abs_exgradp)'
      write(*,*) eru1ir(1), erp1ir(1) 
   endif
      
  100 format(///,10X,'FINITE ELEMENT ERRORS',                              &
&              /,10X,'=====================',//,                           &
&  '          TIME STEP NO.',I5,',  ITERATION NO. ',I5,/,5X,40('-'),/,    &
&  '          NORM       VELOCITY             PRESSURE ',/,5X,40('-'),/,  &
&  '          W(0,1) ',E12.5,9X,E12.5 ,/, &
&  '          W(0,2) ',E12.5,9X,E12.5 ,/, &
&  '          W(0,i) ',E12.5,9X,E12.5 ,/, &
&  '          W(1,1) ',E12.5,9X,E12.5 ,/, &
&  '          W(1,2) ',E12.5,9X,E12.5 ,/, &
&  '          W(1,i) ',E12.5,9X,E12.5 ,/,5X,40('-'))


end subroutine plcd_exaerr
