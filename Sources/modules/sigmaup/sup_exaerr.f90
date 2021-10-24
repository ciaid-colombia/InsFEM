subroutine sup_exaerr(a)    
!-----------------------------------------------------------------------
!****f* ThreeField/sup_exaerr
! NAME 
!    sup_exaerr
! DESCRIPTION
!    This routine computes the FEM errors (referred to an analytical
!    solution defined in exacso.f90). The errors are normalized by the 
!    appropriate norm of the exact solution, except when this norm is zero.      
! USES
! USED BY
!    sup_output
!***
!-----------------------------------------------------------------------
   use typre
   use Mod_ThreeField
   use Mod_Mesh
   use Mod_Element
   use Mod_memor
   use MPI
   use Mod_SupExacso    
   
   implicit none
   class(ThreeFieldNSProblem) :: a   
   type(MemoryMan) :: Memor
   type(FemMesh)    :: Mesh  
   class(FiniteElement), pointer :: e => NULL()
   type(SupExacso)  :: exacso     
   
   
   !Numerical Gaus Values
   real(rp), allocatable    :: vegau(:),grave(:,:),grsig(:,:)
   real(rp), allocatable    :: gradp(:),gpcod(:),sigau(:) 
   real(rp)    :: prgau(1)
   
   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:),exprg(:),exsig(:),exsigr(:,:)
   real(rp)                :: expre,dvolu
   
   real(rp)    :: difeu,abvel,diffp,abpre,dummr,zero_rp,weightfactor,difes,absig,Press_average
   !Errors 
   real(rp)    :: erp01(2),erp02(2),erp0i(2),erp11(2),erp12(2),erp1i(2)
   real(rp)    :: eru01(2),eru02(2),eru0i(2),eru11(2),eru12(2),eru1i(2) 
   real(rp)    :: ers01(2),ers02(2),ers0i(2),ers11(2),ers12(2),ers1i(2)   
   
   !Reduced Errors 
   real(rp)    :: erp01r(2),erp02r(2),erp0ir(2),erp11r(2),erp12r(2),erp1ir(2)
   real(rp)    :: eru01r(2),eru02r(2),eru0ir(2),eru11r(2),eru12r(2),eru1ir(2)
   real(rp)    :: ers01r(2),ers02r(2),ers0ir(2),ers11r(2),ers12r(2),ers1ir(2)   

   integer(ip) :: ievab     ! Indices and dimensions
   integer(ip) :: pevab,pevat
   integer(ip) :: pelty,nelem,npoin,nelty,ndime,dummi,idime,jdime,ielem,igaus,ntens
   integer(ip), pointer     :: nnode(:),ngaus(:),ltopo(:) 

   
   !Elemental Values
   real(rp), allocatable :: elveloc(:,:),elpress(:),elsig(:,:) 
   integer :: ierr   
   integer(ip)              :: i,root,icount,inode,npoinLocal,auxstep
   real(rp)  :: auxdt1,auxdt2,auxdt3
   
   real(rp), pointer :: coord(:)
   real(rp), allocatable :: exactPressureError(:)
   integer(ip) :: ipoin
   
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
   ers01=0.0_rp
   ers02=0.0_rp
   ers0i=0.0_rp
   ers11=0.0_rp
   ers12=0.0_rp
   ers1i=0.0_rp

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
   ers01r=0.0_rp
   ers02r=0.0_rp
   ers0ir=0.0_rp
   ers11r=0.0_rp
   ers12r=0.0_rp
   ers1ir=0.0_rp

   zero_rp=0.0_rp
   
   call a%Mesh%GetNdime(ndime)
   ntens=(ndime-1)*(ndime-1)+2

   call a%Mesh%GetElementTypeSizes(nelty,nnode,ngaus,ltopo)
   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNpoin(npoin) 
   call a%Mesh%GetNpoinLocal(npoinLocal)
   !Element object Initialization
   call a%Mesh%ElementAlloc(e,Memor,'DefaultRule','sup_exaerr')   
   

   ! Allocate exact components
   call a%Memor%alloc(ndime,gpcod,'gpcod','sup_exaerr')
   call a%Memor%alloc(ndime,exvel,'exvel','sup_exaerr') 
   call a%Memor%alloc(ndime,vegau,'vegau','sup_exaerr')
   call a%Memor%alloc(ndime,exprg,'exprg','sup_exaerr')  
   call a%Memor%alloc(ndime,gradp,'gradp','sup_exaerr')
   call a%Memor%alloc(ndime,ndime,exveg,'exveg','sup_exaerr')
   call a%Memor%alloc(ndime,ndime,grave,'grave','sup_exaerr')
   call a%Memor%alloc(ntens,exsig,'exsig','sup_exaerr') 
   call a%Memor%alloc(ntens,sigau,'sigau','sup_exaerr')
   call a%Memor%alloc(ntens,ndime,exsigr,'exsigr','sup_exaerr') 
   call a%Memor%alloc(ntens,ndime,grsig,'grsig','sup_exaerr')
   
   !Allocate exact Values
   call a%Memor%alloc(ndime,e%mnode,elveloc,'elveloc','sup_exaerr')     
   call a%Memor%alloc(e%mnode,elpress,'elpress','sup_exaerr') 
   call a%Memor%alloc(ntens,e%mnode,elsig,'elsig','sup_exaerr')       
    
   if(a%kfl_timei==0)then
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
         call e%gather(ntens,elsig,a%sigma)

         call e%elmdcg

         gauss_points: do igaus=1,e%pgaus
            e%igaus = igaus      

            call e%elmder
            dvolu = e%weigp(e%igaus)*e%detjm
            

            ! Velocity, pressure, velocity gradient and pressure gradient
            
            !Pressure and Velocity gradient
            call e%gradient(1,elpress,gradp)          
            call e%gradient(ndime,elveloc,grave) 
            call e%gradient(ntens,elsig,grsig)          

            !Gauss point values
            call e%interpg(ndime,elveloc,vegau)  
            call e%interpg(ndime,e%elcod,gpcod)
            call e%interpg(1,elpress,prgau)  
            call e%interpg(ntens,elsig,sigau)       

      
            ! Exact solution
            call exacso%sup_ComputeSolution(ndime,gpcod,a%ctime,a%LogFormulation,a)
            call exacso%sup_GetPressure(ndime,expre,exprg)         
            call exacso%sup_GetVelocity(ndime,exvel,exveg) 
            if (a%LogFormulation/=0) then
               call exacso%sup_GetPsi(ndime,exsig,exsigr)
            else   
               call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
            end if  

               
         ! Errors
            diffp = abs(prgau(1)-expre)
            abpre = abs(expre)
            erp01(1) = erp01(1) + weightfactor*diffp*dvolu
            erp02(1) = erp02(1) + weightfactor*diffp*diffp*dvolu
            erp0i(1) = max(erp0i(1),diffp)
            erp01(2) = erp01(2) + weightfactor*abpre*dvolu
            erp02(2) = erp02(2) + weightfactor*abpre*abpre*dvolu       
            erp0i(2) = max(erp0i(2),expre)

            do idime=1,ntens
               difes = abs(sigau(idime)-exsig(idime))
               absig = abs(exsig(idime)) 
               ers01(1) = ers01(1) + weightfactor*difes*dvolu
               ers02(1) = ers02(1) + weightfactor*difes*difes*dvolu
               ers0i(1) = max(ers0i(1),difes)
               ers01(2) = ers01(2) + weightfactor*absig*dvolu
               ers02(2) = ers02(2) + weightfactor*absig*absig*dvolu
               ers0i(2) = max(ers0i(2),absig)
               do jdime=1,ndime
                  difes = abs(grsig(idime,jdime)-exsigr(idime,jdime))
                  absig = abs(exsigr(idime,jdime))
                  ers11(1) = ers11(1) + weightfactor*difes*dvolu
                  ers12(1) = ers12(1) + weightfactor*difes*difes*dvolu
                  ers1i(1) = max(ers1i(1),difes)
                  ers11(2) = ers11(2) + weightfactor*absig*dvolu
                  ers12(2) = ers12(2) + weightfactor*absig*absig*dvolu
                  ers1i(2) = max(ers1i(2),absig)
               end do
            end do


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
               
      call  MPI_REDUCE( erp01,  erp01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp02,  erp02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp0i,  erp0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )
      
      call  MPI_REDUCE( eru01,  eru01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru02,  eru02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru0i,  eru0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

      call  MPI_REDUCE( ers01,  ers01r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( ers02,  ers02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( ers0i,  ers0ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

      call  MPI_REDUCE( erp11,  erp11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp12,  erp12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( erp1i,  erp1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

      call  MPI_REDUCE( ers11,  ers11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( ers12,  ers12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( ers1i,  ers1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )

      call  MPI_REDUCE( eru11,  eru11r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru12,  eru12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
      call  MPI_REDUCE( eru1i,  eru1ir, 2, MPI_REAL8, MPI_MAX, a%MPIroot,a%MPIcomm, ierr )       
      
      erp02r(1) = sqrt(erp02r(1))
      erp12r(1) = sqrt(erp12r(1))
      eru02r(1) = sqrt(eru02r(1))
      eru12r(1) = sqrt(eru12r(1))
      ers02r(1) = sqrt(ers02r(1))
      ers12r(1) = sqrt(ers12r(1))
      
      if(erp01r(2).gt.zero_rp) erp01r(1) = erp01r(1)/erp01r(2) 
      if(erp02r(2).gt.zero_rp) erp02r(1) = erp02r(1)/sqrt(erp02r(2))
      if(erp0ir(2).gt.zero_rp) erp0ir(1) = erp0ir(1)/erp0ir(2)

      if(eru01r(2).gt.zero_rp) eru01r(1) = eru01r(1)/eru01r(2) 
      if(eru02r(2).gt.zero_rp) eru02r(1) = eru02r(1)/sqrt(eru02r(2))
      if(eru0ir(2).gt.zero_rp) eru0ir(1) = eru0ir(1)/eru0ir(2) 

      if(ers01r(2).gt.zero_rp) ers01r(1) = ers01r(1)/ers01r(2) 
      if(ers02r(2).gt.zero_rp) ers02r(1) = ers02r(1)/sqrt(ers02r(2))
      if(ers0ir(2).gt.zero_rp) ers0ir(1) = ers0ir(1)/ers0ir(2) 

      if(erp11r(2).gt.zero_rp) erp11r(1) = erp11r(1)/erp11r(2) 
      if(erp12r(2).gt.zero_rp) erp12r(1) = erp12r(1)/sqrt(erp12r(2))
      if(erp1ir(2).gt.zero_rp) erp1ir(1) = erp1ir(1)/erp1ir(2)

      if(eru11r(2).gt.zero_rp) eru11r(1) = eru11r(1)/eru11r(2) 
      if(eru12r(2).gt.zero_rp) eru12r(1) = eru12r(1)/sqrt(eru12r(2))
      if(eru1ir(2).gt.zero_rp) eru1ir(1) = eru1ir(1)/eru1ir(2) 

      if(ers11r(2).gt.zero_rp) ers11r(1) = ers11r(1)/ers11r(2) 
      if(ers12r(2).gt.zero_rp) ers12r(1) = ers12r(1)/sqrt(ers12r(2))
      if(ers1ir(2).gt.zero_rp) ers1ir(1) = ers1ir(1)/ers1ir(2) 
      
      
      
   else
      
      auxstep=1
      if(a%istep>=auxstep)then

         
         !temporal initialization 
         if(a%istep==auxstep)then
            a%erplit=0.0_rp
            a%erulit=0.0_rp
            a%erslit=0.0_rp
         end if          
         
         do ielem=1,nelem  
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
            call e%gather(ntens,elsig,a%sigma)

            call e%elmdcg

            do igaus=1,e%pgaus
               e%igaus = igaus    
               
               call e%elmder

               dvolu = e%weigp(e%igaus)*e%detjm         
             
               !Pressure and Velocity gradient
               call e%gradient(1,elpress,gradp)          
               call e%gradient(ndime,elveloc,grave) 
               call e%gradient(ntens,elsig,grsig)          

               !Gauss point values
               call e%interpg(ndime,elveloc,vegau)  
               call e%interpg(ndime,e%elcod,gpcod)
               call e%interpg(1,elpress,prgau)  
               call e%interpg(ntens,elsig,sigau)       

         
               ! Exact solution
               call exacso%sup_ComputeSolution(ndime,gpcod,a%ctime,a%LogFormulation,a)
               call exacso%sup_GetPressure(ndime,expre,exprg)         
               call exacso%sup_GetVelocity(ndime,exvel,exveg) 
               
               if (a%LogFormulation/=0) then
                  call exacso%sup_GetPsi(ndime,exsig,exsigr)
               else   
                  call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
               end if    

                  
               !Errors
               diffp = abs(prgau(1)-expre)
               abpre = abs(expre)
               erp02(1) = erp02(1) + weightfactor*diffp*diffp*dvolu
               erp02(2) = erp02(2) + weightfactor*abpre*abpre*dvolu       
               
               do idime=1,ntens
                  difes = abs(sigau(idime)-exsig(idime))
                  absig = abs(exsig(idime)) 
                  ers02(1) = ers02(1) + weightfactor*difes*difes*dvolu
                  ers02(2) = ers02(2) + weightfactor*absig*absig*dvolu
                  do jdime=1,ndime
                     difes = abs(grsig(idime,jdime)-exsigr(idime,jdime))
                     absig = abs(exsigr(idime,jdime))
                     ers12(1) = ers12(1) + (1.0_rp/a%dtinv)*weightfactor*difes*difes*dvolu
                     ers12(2) = ers12(2) + (1.0_rp/a%dtinv)*weightfactor*absig*absig*dvolu
                  end do
               end do 

               do idime=1,ndime
                  difeu = abs(vegau(idime)-exvel(idime))
                  abvel = abs(exvel(idime))
                  eru02(1) = eru02(1) + weightfactor*difeu*difeu*dvolu
                  eru02(2) = eru02(2) + weightfactor*abvel*abvel*dvolu
                  
                  diffp = abs(gradp(idime)-exprg(idime))
                  abpre = abs(exprg(idime)) 
                  erp12(1) = erp12(1) + (1.0_rp/a%dtinv)*weightfactor*diffp*diffp*dvolu
                  erp12(2) = erp12(2) + (1.0_rp/a%dtinv)*weightfactor*abpre*abpre*dvolu
                  
                  do jdime=1,ndime
                     difeu = abs(grave(idime,jdime)-exveg(idime,jdime))
                     abvel = abs(exveg(idime,jdime))
                     eru12(1) = eru12(1) + (1.0_rp/a%dtinv)*weightfactor*difeu*difeu*dvolu
                     eru12(2) = eru12(2) + (1.0_rp/a%dtinv)*weightfactor*abvel*abvel*dvolu
                  end do
               end do
               
            end do 
            
         end do 
            
            !L2
            call  MPI_REDUCE( eru02,  eru02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )      
            call  MPI_REDUCE( erp02,  erp02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )       
            call  MPI_REDUCE( ers02,  ers02r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
            !H1
            call  MPI_REDUCE( eru12,  eru12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )           
            call  MPI_REDUCE( erp12,  erp12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
            call  MPI_REDUCE( ers12,  ers12r, 2, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr ) 
            
            !numerador
            !L2
            eru02r(1) = sqrt(eru02r(1))      
            erp02r(1) = sqrt(erp02r(1))
            ers02r(1) = sqrt(ers02r(1))      
            !Grad
            eru12r(1) = eru12r(1)      
            erp12r(1) = erp12r(1)      
            ers12r(1) = ers12r(1)
            !denominador
            !L2
            eru02r(2) = sqrt(eru02r(2))      
            erp02r(2) = sqrt(erp02r(2))
            ers02r(2) = sqrt(ers02r(2))
            !H1
            eru12r(2) = eru12r(2)      
            erp12r(2) = erp12r(2)
            ers12r(2) = ers12r(2)
            
            !temporal errors Linf(L2)
            
            auxdt1=a%erplit(1)
            auxdt2=a%erulit(1)
            auxdt3=a%erslit(1)
            
            a%erulit(1)=max(a%erulit(1),eru02r(1))      
            a%erplit(1)=max(a%erplit(1),erp02r(1)) 
            a%erslit(1)=max(a%erslit(1),ers02r(1))
            
            a%erulit(2)=max(a%erulit(2),eru02r(2))
            a%erplit(2)=max(a%erplit(2),erp02r(2)) 
            a%erslit(2)=max(a%erslit(2),ers02r(2))
            
            !You need to take the square root of this value at the final of the evolution problem
            a%eruh1t(1)=a%eruh1t(1) + eru12r(1)
            a%erph1t(1)=a%erph1t(1) + erp12r(1)       
            a%ersh1t(1)=a%ersh1t(1) + ers12r(1)
            
            a%eruh1t(2)=a%eruh1t(2) + eru12r(2)
            a%erph1t(2)=a%erph1t(2) + erp12r(2)
            a%ersh1t(2)=a%ersh1t(2) + ers12r(2)  
            
            if(a%erulit(1).gt.auxdt2) write(*,*) a%istep   
         
      end if
      
   end if
      

   ! deallocate elemental values
   call a%Memor%dealloc(ndime,e%mnode,elveloc,'elveloc','sup_exaerr')
   call a%Memor%dealloc(ntens,e%mnode,elsig,'elsig','sup_exaerr')     
   call a%Memor%dealloc(e%mnode,elpress,'elpress','sup_exaerr')     
   
   ! deallocate exact components
   call a%Memor%dealloc(ndime,gpcod,'gpcod','sup_exaerr')
   call a%Memor%dealloc(ndime,exvel,'exvel','sup_exaerr') 
   call a%Memor%dealloc(ndime,vegau,'vegau','sup_exaerr')
   call a%Memor%dealloc(ndime,exprg,'exprg','sup_exaerr')  
   call a%Memor%dealloc(ndime,gradp,'gradp','sup_exaerr')
   call a%Memor%dealloc(ndime,ndime,exveg,'exveg','sup_exaerr')
   call a%Memor%dealloc(ndime,ndime,grave,'grave','sup_exaerr')
   call a%Memor%dealloc(ntens,exsig,'exsig','sup_exaerr') 
   call a%Memor%dealloc(ntens,sigau,'sigau','sup_exaerr')
   call a%Memor%dealloc(ntens,ndime,exsigr,'exsigr','sup_exaerr') 
   call a%Memor%dealloc(ntens,ndime,grsig,'grsig','sup_exaerr')
   
   
   call a%Mesh%ElementDealloc(e,Memor,'DefaultRule','sup_exaerr') 
   
   if (a%MPIrank == a%MPIroot) then 
   
      if(a%kfl_timei==0)then

         write(*,*) 'diferencia u   diferencia p   diferencia s'  
         write(*,*)  eru01r(1),  erp01r(1),  ers01r(1)
         write(*,*) 'L2 u     L2 p     L2 s'   
         write(*,*) eru02r(1),   erp02r(1),  ers02r(1)
         write(*,*) 'Max(abs_exvel) Max(abs_exp)   Max(abs_exsig)' 
         write(*,*) eru0ir(1),   erp0ir(1),  ers0ir(1)
         write(*,*) 'dif_gradu    dif_gradu  dif_grads'
         write(*,*) eru11r(1),   erp11r(1),  ers11r(1)
         write(*,*) 'L2 gradu     L2 gradp   L2 grads'
         write(*,*) eru12r(1),   erp12r(1),  ers12r(1)
         write(*,*) 'Max(abs_exgradu)  Max(abs_exgradp)  Max(abs_exgradp)'
         write(*,*) eru1ir(1),   erp1ir(1),  ers1ir(1) 
         
      elseif(a%kfl_timei/=0)then
         write(*,*)    'Num_Linf(L2)_u         Num_Linf(L2)_p         Num_Linf(L2)_s'  
         write(*,*)  a%erulit(1),  a%erplit(1),  a%erslit(1)
         write(*,*)    'Den_Linf(L2)_u         Den_Linf(L2)_p         Den_Linf(L2)_s'  
         write(*,*)  a%erulit(2),  a%erplit(2),  a%erslit(2)      
         write(*,*)    'Num_L2(H1)_u         Num_L2(H1)_p         Num_L2(H1)_s'  
         write(*,*)  a%eruh1t(1),  a%erph1t(1),  a%ersh1t(1)
         write(*,*)    'Den_L2(H1)_u         Den_L2(H1)_p         Den_L2(H1)_s'  
         write(*,*)  a%eruh1t(2),  a%erph1t(2),  a%ersh1t(2)        
         
      end if
      
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


end subroutine sup_exaerr
