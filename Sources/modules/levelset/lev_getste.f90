subroutine lev_getste(a,dtinv)
   use typre
   use Mod_LevelSet
   use Mod_Element
   use Mod_GatherScatterDtcri
   use MPI
   implicit none
   class(LevelSetProblem) :: a
   real(rp) :: dtinv
   
   real(rp), allocatable :: elvel(:,:)
   real(rp), allocatable :: gpvel(:)
   class(FiniteElement), pointer :: e => NULL()
   
   real(rp) :: dtcri,dtmin,gpvno,hclen,rdinv
   integer(ip) :: ielem,nelem,irank,ndime
   
   interface
      subroutine lev_coupling_outerr(a)
         use Mod_LevelSet
         implicit none
         class(LevelSetProblem) :: a
      end subroutine
   
   end interface
   
   !Check if the coupling between modules is done
   call lev_coupling_outerr(a)
   
   
   call a%Timer%Total%Tic
   call a%Timer%Getste%Tic

   !Critical time step
   if(a%kfl_timei/=0.and.a%kfl_stead/=1) then
   
      !Dimensions and general variables
      call a%Mesh%GetNdime(ndime) 
      call a%Mesh%GetNelem(nelem)
      rdinv=1.0_rp/real(ndime)
      gpvno = 0.0_rp
      dtmin = 1e6
      
      !Element Initialization
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lev_getste')

      call a%Memor%alloc(ndime,e%mnode,elvel,'elvel','lev_getste')
      call a%Memor%alloc(ndime,gpvel,'gpvel','lev_getste') 

      
      do ielem = 1,nelem
         
         !Load Element
         call a%Mesh%ElementLoad(ielem,e)
         
         !Derivative (and detjm) at the center of gravity
         call e%elmdcg
         
         !Element length
         hclen=(e%weicg*e%detjm)**rdinv
         
         !Critical Time step computation
         dtcri=0.0_rp

         if (a%kfl_advec == 1) then
            call e%gather(ndime,elvel,a%veloc(:,:))          
            call e%interpc(ndime,elvel,gpvel)
            call vecnor(gpvel,e%ndime,gpvno,2)
            !2u/h (Typicall value for the advection stabilized constant) 
         else
            gpvno = 0.0_rp
         endif
         dtcri=dtcri+2.0_rp*gpvno/hclen
         if(dtcri <= zelev) then 
             dtcri = 1e9
         else 
             dtcri=1.0_rp/dtcri
         end if
         dtmin=min(dtmin,dtcri)
      end do
      
      a%dtcri = dtmin
      
      !Gather Dtcri from all processes, find minimum and scatter
      call php_GatherScatterDtcri(a)
      
      a%dtinv = 1.0_rp/(a%dtcri*a%safet)
      dtinv = a%dtinv
      
     !Deallocate variables
     call a%Memor%dealloc(ndime,e%mnode,elvel,'elvel','lev_getste')
     call a%Memor%dealloc(ndime,gpvel,'gpvel','lev_getste')
     
     call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','lev_getste')
!      call e%dealloc(a%Memor,'lev_getste')     
   end if
   
   call a%Timer%Getste%Toc
   call a%Timer%Total%Toc
  
   
end subroutine
