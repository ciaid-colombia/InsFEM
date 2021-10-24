subroutine nsi_exaerr(b)    
!-----------------------------------------------------------------------
!****f* Nstinc/nsi_exaerr
! NAME 
!    nsi_exaerr
! DESCRIPTION
!    This routine computes the FEM errors (referred to an analytical
!    solution defined in exacso.f90). The errors are normalized by the 
!    appropriate norm of the exact solution, except when this norm is zero.      
! USES
! USED BY
!    nsi_output
!***
!-----------------------------------------------------------------------
   use MPI
   use typre
   use Mod_Element
   use Mod_NavierStokes   
   use Mod_NsiExacso 
   use Mod_nsm_BaseElmope
   
   implicit none
   class(NavierStokesProblem), target  :: b
   integer(ip), save :: ipass=0
   integer           :: ierr   

   !Numerical Gaus Values
   real(rp), allocatable :: vegau(:),grave(:,:),gradp(:),gpcod(:)
   real(rp)              :: prgau(1),dvolu,divve
   
   !Exact Values
   real(rp), allocatable :: exvel(:),exveg(:,:),exprg(:)
   real(rp)              :: expre,exdiv
   
   !Errors 
   real(rp), allocatable :: diff_adv(:),diff_dif(:,:)
   real(rp)              :: diff_div,err_adv,err_div,err_dif,terror,err_adv0,err_div0,err_dif0
   real(rp), allocatable :: error(:)

   integer(ip)           :: npoin,nelty
   integer(ip), pointer  :: nnode(:)=>NULL(),ngaus(:)=>NULL(),ltopo(:)=>NULL() 

   !Elemental Values
   real(rp), allocatable :: elveloc(:,:),elpress(:) 
   integer(ip)           :: icount,inode,npoinLocal
   real(rp)              :: weightfactor
   
   a => b

   call a%Mesh%GetNdime(ndime)   
   call a%Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNpoinLocal(npoinLocal)
   !Element object Initialization
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsi_exaerr')   
   
   !Allocate exact components
   call a%Memor%alloc(ndime,gpcod,'gpcod','nsi_exaerr') 
   call a%Memor%alloc(ndime,exvel,'exvel','nsi_exaerr')
   call a%Memor%alloc(ndime,vegau,'vegau','nsi_exaerr')
   call a%Memor%alloc(ndime,exprg,'exprg','nsi_exaerr')  
   call a%Memor%alloc(ndime,gradp,'gradp','nsi_exaerr')
   call a%Memor%alloc(ndime,ndime,exveg,'exveg','nsi_exaerr')
   call a%Memor%alloc(ndime,ndime,grave,'grave','nsi_exaerr')
   
   !Allocate elemental values
   call a%Memor%alloc(ndime,e%mnode,elveloc,'elveloc','nsi_exaerr')     
   call a%Memor%alloc(e%mnode,elpress,'elpress','nsi_exaerr')     
    
   !Allocate errors
   call a%Memor%alloc(ndime,diff_adv,'diff_adv','nsi_exaerr') 
   call a%Memor%alloc(ndime,ndime,diff_dif,'diff_dif','nsi_exaerr') 
   call a%Memor%alloc(nelem,error,'error','nsi_exaerr')

! Initializations
   diff_div = 0.0_rp
   diff_dif = 0.0_rp
   diff_adv = 0.0_rp

   err_div = 0.0_rp
   err_dif = 0.0_rp
   err_adv = 0.0_rp
   terror   = 0.0_rp

   err_div0 = 0.0_rp
   err_dif0 = 0.0_rp
   err_adv0 = 0.0_rp

   acden=a%MatProp(1)%densi
   
   if(a%MatProp(1)%lawvi==0)then
      acvis = a%MatProp(1)%visco
   elseif(a%MatProp(1)%lawvi==1)then
      acvis = a%MatProp(1)%LawViParam(1)
   elseif(a%MatProp(1)%lawvi==2)then
      call runend('nsi_ComputeSolution: exact solution not implemented for Carreau model')    
   end if
     
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
   
      error(ielem) = 0.0_rp
   
      !Gather operations    
      call e%gather(ndime,elveloc,a%veloc)
      call e%gather(1_ip,elpress,a%press)

      call e%elmdcg
      call e%elmlen
      call elmchl(e,a%kfl_advec,elveloc,chale,a%kfl_hdifumin)

      gauss_points: do igaus=1,e%pgaus
         e%igaus = igaus      

         dvolu = e%weigp(e%igaus)*e%detjm
      
         !Velocity, pressure, velocity gradient and pressure gradient
      
         !Pressure and Velocity gradient
         call e%gradient(1,elpress,gradp)          
         call e%gradient(ndime,elveloc,grave)           
         call e%divergence(elveloc,divve)           

         !Gauss point values
         call e%interpg(ndime,elveloc,vegau)       
         call e%interpg(ndime,e%elcod,gpcod)
         call e%interpg(1,elpress,prgau) 
      
         !Exact solution
         call exacso%nsi_ComputeSolution(ndime,gpcod,a%ctime,a)
         call exacso%nsi_GetPressure(ndime,expre,exprg)         
         call exacso%nsi_GetVelocity(ndime,exvel,exveg)
         call trace(ndime,exveg,exdiv)

         if (a%kfl_advco == 1) then 
            call vecnor(a%advco,e%ndime,gpvno,2)
            call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)
         endif
         call ComputeTauDiv(e,a%staco(4),chale,timom,tidiv)

         !Errors
         diff_adv = abs(matmul(vegau-exvel,grave-exveg)+gradp-exprg)
         err_adv0 = err_adv + weightfactor*acden*tidiv*dot_product(diff_adv,diff_adv)*dvolu
         diff_div = abs(divve - exdiv)
         err_div0 = err_div + weightfactor*acden*timom*diff_div*diff_div*dvolu
         diff_dif = abs(grave-exveg)
         err_dif0 = err_dif + weightfactor*acvis*sum(diff_dif*diff_dif)*dvolu
         error(ielem) = error(ielem) + (acvis*sum(diff_dif*diff_dif) + acden*timom*diff_div*diff_div + acden*tidiv*dot_product(diff_adv,diff_adv))*dvolu
      end do gauss_points
   
   end do elements

   error = sqrt(error)
   call a%FilePostpr%postgp(error,'Error',a%istep,a%ctime,a%Mesh)
   
   call  MPI_ALLREDUCE(err_div0,err_div,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
   call  MPI_ALLREDUCE(err_dif0,err_dif,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
   call  MPI_ALLREDUCE(err_adv0,err_adv,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
   terror = sqrt(err_adv+err_div+err_dif)
   
   if (a%MPIrank == a%MPIroot) then 
      if(ipass==0) then
         ipass=1
         write(a%lun_error,200)
      end if
      write(a%lun_error,201) a%ctime,terror
      if (a%kfl_flush == 1) call flush(a%lun_error)
   endif
   
   !Deallocate errors
   call a%Memor%dealloc(ndime,diff_adv,'diff_adv','nsi_exaerr') 
   call a%Memor%dealloc(ndime,ndime,diff_dif,'diff_dif','nsi_exaerr') 
   call a%Memor%dealloc(nelem,error,'error','nsi_exaerr')

   !Deallocate elemental values
   call a%Memor%dealloc(ndime,e%mnode,elveloc,'elveloc','nsi_exaerr')     
   call a%Memor%dealloc(e%mnode,elpress,'elpress','nsi_exaerr')     
   
   !Deallocate exact components
   call a%Memor%dealloc(ndime,gpcod,'gpcod','nsi_exaerr') 
   call a%Memor%dealloc(ndime,exvel,'exvel','nsi_exaerr')
   call a%Memor%dealloc(ndime,vegau,'vegau','nsi_exaerr')
   call a%Memor%dealloc(ndime,exprg,'exprg','nsi_exaerr')  
   call a%Memor%dealloc(ndime,gradp,'gradp','nsi_exaerr')
   call a%Memor%dealloc(ndime,ndime,exveg,'exveg','nsi_exaerr')
   call a%Memor%dealloc(ndime,ndime,grave,'grave','nsi_exaerr')     
   
   call a%Mesh%ElementDealloc(e,a%Memor,'DefaultRule','nsi_exaerr')      
      
   200   format('$ ',5x,' Time',5x,'Error 3')
   201   format(2(2x,E12.5))

end subroutine nsi_exaerr
