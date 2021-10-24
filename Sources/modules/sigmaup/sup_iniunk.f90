subroutine sup_iniunk(a)
!NAME 
!   sup_iniunk
!DESCRIPTION
!   This routine sets up the initial conditions for the a%velocity and a%sigma.
!   If this is a restart, initial conditions are loaded somewhere else
!   but Dirichlet boundary conditions are still loaded here.
!-----------------------------------------------------------------------
   use typre
   use Mod_ThreeField
   use Mod_Mesh  
   use Mod_memor
   use Mod_SupExacso 
   use Mod_LogOperations
   implicit none
   class(ThreeFieldNSProblem) :: a
   type(FemMesh)    :: Mesh  
   type(MemoryMan)  :: Memor
   type(SupExacso)  :: exacso 
   
   real(rp)    :: vmodu,vfreq,dummr
   integer(ip) :: dummi
   integer(ip) :: icomp,ipoin,idime,ndime,npoin,ntens
   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:),exprg(:),exsig(:),exsigr(:,:),psi(:), sigma_aux(:)
   real(rp)                :: auxL,beta, lamb,lamb0, acvis
   real(rp)                :: expre,wdtime,wtimebdf3,wtimebdf2,aux(3)
   !Nodal coordinates
   real(rp), pointer       :: coord(:)
   
   !Todo multy materials
   integer(ip) :: imat=1
   
   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNpoin(npoin)
   
   ntens=(ndime-1)*(ndime-1)+2    
   
   
      if(a%MatProp(imat)%lawvi>=0)then   
         do ipoin=1,npoin
            do idime=1,ndime
               if(a%kfl_fixno(idime,ipoin)==1) then               
                  a%veloc(idime,ipoin,a%ncomp) = a%bvess(idime,ipoin,1)
                  if(a%ncomp>=4) a%veloc(idime,ipoin,4) = a%bvess(idime,ipoin,1)  ! bdf3
                  if(a%ncomp>=3) a%veloc(idime,ipoin,3) = a%bvess(idime,ipoin,1)  ! bdf2 & bdf3 
               end if
            end do
         end do
      
     
      elseif(a%MatProp(imat)%lawvi<0)then
      
         if (a%LogFormulation/=0) then 
            call a%Memor%alloc(ntens,psi,'psi','sup_iniunk')  
            call a%Memor%alloc(ntens,sigma_aux,'sigma_aux','sup_iniunk') 
            acvis = a%MatProp(1)%LawViParam(1) 
            beta  = a%MatProp(1)%LawViParam(2)     
            
            if (a%kfl_exacs/=0) then
               lamb0=a%MatProp(1)%LawViParam(5)
            else   

               lamb0=a%MatProp(imat)%LawViParam(3)*a%MatProp(1)%LawviParam(5) 
            end if 

            if (lamb0==0.0_rp) lamb0=0.01_rp
            auxL=((1.0_rp-beta)*acvis)/lamb0
         end if   
         
         
         do ipoin=1,npoin
            do idime=1,ndime
               if(a%kfl_fixno(ntens+idime,ipoin)==1) then
                  a%veloc(idime,ipoin,a%ncomp) = a%bvess(ntens+idime,ipoin,1)
                  if(a%ncomp>=4) a%veloc(idime,ipoin,4) = a%bvess(ntens+idime,ipoin,1)  ! bdf3
                  if(a%ncomp>=3) a%veloc(idime,ipoin,3) = a%bvess(ntens+idime,ipoin,1)  ! bdf2 & bdf3
               end if
            end do
            
            if(a%kfl_fixno(1,ipoin)==1) then
                
               if (a%LogFormulation/=0) then 
                  psi=0.0_rp
                  sigma_aux(:)= a%bvess(1:ntens,ipoin,1)
                  
                  call sup_LogarithmicTransformation(ndime,ntens,auxL,sigma_aux,psi)
                  
                  do idime=1,ntens  
                     a%bvess(idime,ipoin,1) = psi(idime)
                     a%sigma(idime,ipoin,a%ncomp) = psi(idime)
                     if(a%ncomp>=4) a%sigma(idime,ipoin,4) = psi(idime) ! bdf3
                     if(a%ncomp>=3) a%sigma(idime,ipoin,3) = psi(idime) ! bdf2 & bdf3 
                  end do   

               else
                  do idime=1,ntens
                     a%sigma(idime,ipoin,a%ncomp) = a%bvess(idime,ipoin,1)
                     if(a%ncomp>=4) a%sigma(idime,ipoin,4) = a%bvess(idime,ipoin,1)  ! bdf3
                     if(a%ncomp>=3) a%sigma(idime,ipoin,3) = a%bvess(idime,ipoin,1)  ! bdf2 & bdf3   
                  end do
               end if 
            end if
            
         end do
         
         if (a%LogFormulation/=0) then
            call a%Memor%dealloc(ntens,psi,'psi','sup_iniunk')  
            call a%Memor%dealloc(ntens,sigma_aux,'sigma_aux','sup_iniunk')          
         end if
         
     end if
      
      if(a%kfl_incnd==1) then
         call runend('sup_iniunk: Initial conditions not implemented')
      end if

      if(a%kfl_exacs/=0) then

         !Allocate exact components
         call a%Memor%alloc(ndime,exvel,'exvel','sup_iniunk')  
         call a%Memor%alloc(ndime,ndime,exveg,'exveg','sup_iniunk')     
         call a%Memor%alloc(ntens,exsig,'exsig','sup_iniunk')  
         call a%Memor%alloc(ntens,ndime,exsigr,'exsigr','sup_iniunk')
         call a%Memor%alloc(ndime,exprg,'exprg','sup_iniunk')    
      
         if(a%kfl_timei==0)then   
            do ipoin = 1,npoin 
               call a%Mesh%GetPointCoord(ipoin,coord)

               call exacso%sup_ComputeSolution(ndime,coord,a%ctime,a%LogFormulation,a)
               call exacso%sup_GetVelocity(ndime,exvel,exveg)
                                                            
               do idime=1,ndime
                  a%veloc(idime,ipoin,a%ncomp) = exvel(idime)
                  
               end do
            end do
         elseif(a%kfl_timei/=0)then
            do ipoin = 1,npoin 
               call a%Mesh%GetPointCoord(ipoin,coord)

               ! Exact solution
               call exacso%sup_ComputeSolution(ndime,coord,a%ctime,a%LogFormulation,a)
               call exacso%sup_GetPressure(ndime,expre,exprg)
               call exacso%sup_GetVelocity(ndime,exvel,exveg) 
               if (a%LogFormulation/=0) then
                  call exacso%sup_GetPsi(ndime,exsig,exsigr)
               else   
                  call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
               end if                
               
               !Pressure analytical value
               a%press(ipoin,1)   = expre 
               a%press(ipoin,2)   = expre
               do idime=1,ndime
                  a%veloc(idime,ipoin,1) = exvel(idime)
                  a%veloc(idime,ipoin,2) = exvel(idime)                     
               end do
               do idime=1,ntens
                  a%sigma(idime,ipoin,1) = exsig(idime)
                  a%sigma(idime,ipoin,2) = exsig(idime)                
               end do
               
               if(a%ncomp>=3)then
                  wdtime=1.0_rp/a%dtinv
                  wtimebdf2=a%ctime-wdtime                  
                  ! Exact solution  bdf2 => ctime-dtime
                  call exacso%sup_ComputeSolution(ndime,coord,wtimebdf2,a%LogFormulation,a)
                  call exacso%sup_GetVelocity(ndime,exvel,exveg) 
                  if (a%LogFormulation/=0) then
                     call exacso%sup_GetPsi(ndime,exsig,exsigr)
                  else   
                     call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
                  end if 
                  do idime=1,ndime
                     a%veloc(idime,ipoin,3) = exvel(idime)  ! bdf2 & bdf3
                  end do
                  do idime=1,ntens
                     a%sigma(idime,ipoin,3) = exsig(idime)  ! bdf2 & bdf3
                  end do
                  a%press(ipoin,3) = expre                  
               elseif(a%ncomp>=4)then               
                  wtimebdf3=a%ctime-2.0_rp*wdtime                  
                  ! Exact solution  bdf2 => ctime-dtime
                  call exacso%sup_ComputeSolution(ndime,coord,wtimebdf3,a%LogFormulation,a)
                  call exacso%sup_GetVelocity(ndime,exvel,exveg) 
                  if (a%LogFormulation/=0) then
                     call exacso%sup_GetPsi(ndime,exsig,exsigr)
                  else   
                     call exacso%sup_GetStress(ndime,a%LogFormulation,exsig,exsigr)
                  end if 
                  
                  do idime=1,ndime
                     a%veloc(idime,ipoin,4) = exvel(idime)  ! bdf3
                  end do
                  
                  do idime=1,ntens
                     a%sigma(idime,ipoin,4) = exsig(idime)  ! bdf3
                  end do
                  a%press(ipoin,4) = expre 
                  
               endif              
            end do
         end if       

         ! Deallocate exact components
         call a%Memor%dealloc(ndime,exvel,'exvel','sup_iniunk') 
         call a%Memor%dealloc(ndime,ndime,exveg,'exveg','sup_iniunk')
         call a%Memor%dealloc(ntens,exsig,'exsig','sup_iniunk')  
         call a%Memor%dealloc(ntens,ndime,exsigr,'exsigr','sup_iniunk')
         call a%Memor%dealloc(ndime,exprg,'exprg','sup_iniunk')    

      end if
   
   !Assign u(n,i,*) <-- u(n-1,*,*), initial guess after initialization (or reading restart)
   do icomp = 1,a%ncomp-1
      a%veloc(1:ndime,1:npoin,icomp) = a%veloc(1:ndime,1:npoin,a%ncomp)
      a%press(1:npoin,icomp) = a%press(1:npoin,a%ncomp) 
      a%sigma(1:ntens,1:npoin,icomp) = a%sigma(1:ntens,1:npoin,a%ncomp)      
   enddo  
   
   
   call a%Ifconf

end subroutine sup_iniunk

