subroutine nsc_ExactSolRefCriteria(b,error)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_NSCompressibleImplicit
   use Mod_nsc_BaseElmope
   use Mod_NscExacso
   use MPI
   implicit none
   class(NSCompressibleImplicitProblem), target :: b
   real(rp) :: error(:)
   
   integer(ip) :: ndime,ierr
   integer(ip) :: idime,npoinLocal
   real(rp) :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp) :: subscalesnorm
   real(rp) :: auxtm, auxte,invcp,invcv,denominator,entropy,denominator1,entropy1
   real(rp) :: gppress(1)
   real(rp) :: gpdens(1), gpenerg(1)
   real(rp), allocatable :: gpcod(:) 
   real(rp), allocatable :: eleden(:), elemom(:,:), eleene(:),gpmome(:)
   real(rp), allocatable :: elepre(:)
   
   !Exact Values
   real(rp), allocatable   :: exmom(:),exmog(:,:),exdeg(:),exeng(:)
   real(rp)                :: exden,exene
   real(rp), allocatable   :: exprg(:)
   real(rp)                :: expre

   !Error Values
   real(rp)    :: diffp,diffd,momnorm,diffe
   real(rp), allocatable   ::diffmom(:)

   a=> b


   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%mnode,eleden,'eleden','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elemom,'elemom','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%mnode,eleene,'eleene','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%mnode,elepre,'elepre','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,gpmome,'gpmome','nsc_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,gpcod,'gpcod','nsc_ExactSolRefCriteria')     
   
   ! Allocate Exact Values
   call a%Memor%alloc(e%ndime,exmom,'exmom','nsc_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,exdeg,'exdeg','nsc_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,e%ndime,exmog,'exmog','nsc_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,exeng,'exeng','nsc_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,diffmom,'diffmom','nsc_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,exprg,'exprg','nsc_ExactSolRefCriteria')     

   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)

   if(accph>zensi) then
      invcp = 1.0_rp/accph
   else
      invcp = 0.0_rp
   end if
   if(accvh>zensi) then
      invcv = 1.0_rp/accvh
   else
      invcv = 0.0_rp
   end if

   do ielem = 1,nelem
   
      call a%Mesh%ElementLoad(ielem,e)  
      
      call e%elmdcg
      call e%elmlen
      
      call e%gather(1,eleden,a%densf(:,1))
      call e%gather(e%ndime,elemom,a%momen(:,:,1))
      call e%gather(1,eleene,a%energ(:,1))
      call e%gather(1,elepre,a%press(:,1))

      call nsc_ComputeElementVelocity(e,eleden,elemom,elvel)
      ! Compute the characteristic length chale
      call elmchl(e,1_ip,elvel,chale)

      error(ielem) = 0.0_rp

      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         call e%elmlen
         dvol = e%weigp(e%igaus)*e%detjm
   
         call e%interpg(e%ndime,e%elcod,gpcod)
         call e%interpg(1_ip,eleden,gpdens(1))
         call e%interpg(e%ndime,elemom,gpmome(:))
         call e%interpg(1_ip,eleene,gpenerg(1))
         call e%interpg(1_ip,elepre,gppress(1))

         ! Exact solution
         call exacso%nsc_ComputeSolution(e%ndime,gpcod,a)
         call exacso%nsc_GetDensity(e%ndime,exden,exdeg)         
         call exacso%nsc_GetMomentum(e%ndime,exmom,exmog) 
         call exacso%nsc_GetEnergy(e%ndime,exene,exeng)         
         call exacso%nsc_GetPressure(e%ndime,expre,exprg)         

         !Scaled L_2 norm
         if (a%ErrorEstimatorTypeOfSubscales == 1) then
      
            !Advection velocity norm 
            call vecnor(gpmome,e%ndime,gpmno,2)
            gptem = (gpenerg(1) - (gpmno*gpmno/(gpdens(1)*2.0_rp)))*invcv/gpdens(1)
            !Sound speed
            call nsc_ComputeSoundSpeed(accph,accvh,gptem,gpspd)
      
            call nsc_ComputeTau(e,gpspd,acvis/gpdens(1),actco*invcp/gpdens(1),gpmno/gpdens(1),a%staco,chale,timom)
            
            diffd = abs(gpdens(1)-exden)
            diffe = abs(gpenerg(1)-exene)
            do idime=1,e%ndime
               diffmom(idime) = abs(gpmome(idime)-exmom(idime))
            end do 
      
            call vecnor(diffmom(1:e%ndime),e%ndime,momnorm,2)
      
            a%referencevelocity = gpspd + gpmno/gpdens(1)

            error(ielem) = error(ielem) + (diffd**2)*(a%referencevelocity**2)*dvol/timom(1)
            error(ielem) = error(ielem) + (momnorm**2)*dvol/timom(2)
            error(ielem) = error(ielem) + (diffe**2)*dvol/((a%referencevelocity**2)*timom(3))
   
         !Entropy function
         elseif (a%ErrorEstimatorTypeOfSubscales == 2) then

            diffd = abs(gpdens(1)-exden)
            diffp = abs(gppress(1)-expre)

            acgamma = accph/accvh

            denominator = exden**acgamma 

            entropy = accvh*log(expre/denominator)
               
            denominator1 = gpdens(1)**acgamma 

            entropy1 = accvh*log((gppress(1)+a%relpre)/denominator1) 

            error(ielem) = error(ielem) + (entropy-entropy1)*(entropy-entropy1)*dvol/(entropy1*entropy1)
         endif

      enddo

   enddo
   
   
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','nsc_ExactSolRefCriteria')
   call a%Memor%dealloc(e%mnode,eleden,'eleden','nsc_ExactSolRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elemom,'elemom','nsc_ExactSolRefCriteria')
   call a%Memor%dealloc(e%mnode,eleene,'eleene','nsc_ExactSolRefCriteria')
   call a%Memor%dealloc(e%mnode,elepre,'elepre','nsc_ExactSolRefCriteria')
   call a%Memor%dealloc(e%ndime,gpmome,'gpmome','nsc_ExactSolRefCriteria')
   call a%Memor%dealloc(e%ndime,gpcod,'gpcod','nsc_ExactSolRefCriteria')     
   
   ! Allocate Exact Values
   call a%Memor%dealloc(e%ndime,exmom,'exmom','nsc_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,exdeg,'exdeg','nsc_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,e%ndime,exmog,'exmog','nsc_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,exeng,'exeng','nsc_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,diffmom,'diffmom','nsc_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,exprg,'exprg','nsc_ExactSolRefCriteria')     


   !Compute the total error
   TotalEstimatedError0 = 0.0_rp
   call a%Mesh%GetNpoinLocal(npoinLocal)
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      weightfactor = 1.0-real(count(e%lnods(1:e%pnode)>npoinLocal))/real(e%pnode)
   
      TotalEstimatedError0 = TotalEstimatedError0 + error(ielem)*weightfactor

   enddo
   call  MPI_REDUCE(TotalEstimatedError0, TotalEstimatedError, 1, MPI_REAL8, MPI_SUM, a%MPIroot,a%MPIcomm, ierr )
   
   TotalEstimatedError = sqrt(TotalEstimatedError)
   if (a%MPIrank == a%MPIroot) write (*,"(F8.4)",advance='no')TotalEstimatedError
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','nsc_ExactSolRefCriteria')
   
end subroutine   
