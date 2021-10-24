subroutine nsc_pr_ExactSolRefCriteria(b,error)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_NSCompressiblePrimitive
   use Mod_nsc_pr_BaseElmope
   use Mod_NscExacso
   use MPI
   implicit none
   class(NSCompressiblePrimitiveProblem), target :: b
   real(rp) :: error(:)
   
   integer(ip) :: ndime,ierr
   integer(ip) :: idime,npoinLocal
   real(rp) :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp) :: subscalesnorm, denominator,entropy, acgamma,denominator1,entropy1
   real(rp) :: auxtm, auxte
   real(rp) :: gppress(1), gptempe(1)
   real(rp) :: gpdens(1)
   real(rp), allocatable :: gpcod(:) 
   real(rp), allocatable :: elepre(:), elevel(:,:), eletem(:),gpvele(:)
   real(rp), allocatable :: eleden(:)
   
   !Exact Values
   real(rp), allocatable   :: exvel(:),exveg(:,:),exprg(:),exteg(:)
   real(rp)                :: expre,extem
   real(rp), allocatable   :: exdeg(:)
   real(rp)                :: exden

   !Error Values
   real(rp)    :: diffp,diffd,velnorm,difft
   real(rp), allocatable   ::diffvel(:)


   a => b

   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   

   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_pr_ExactSolRefCriteria')
   call a%Memor%alloc(e%mnode,eleden,'eleden','nsc_pr_ExactSolRefCriteria')
   call a%Memor%alloc(e%mnode,elepre,'elepre','nsc_pr_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elevel,'elevel','nsc_pr_ExactSolRefCriteria')
   call a%Memor%alloc(e%mnode,eletem,'eletem','nsc_pr_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,gpvele,'gpvele','nsc_pr_ExactSolRefCriteria')
   call a%Memor%alloc(e%ndime,gpcod,'gpcod','nsc_pr_ExactSolRefCriteria')     
   
   ! Allocate Exact Values
   call a%Memor%alloc(e%ndime,exvel,'exvel','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,exprg,'exprg','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,e%ndime,exveg,'exveg','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,exteg,'exteg','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,diffvel,'diffvel','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%alloc(e%ndime,exdeg,'exdeg','nsc_pr_ExactSolRefCriteria')     

   call a%Mesh%GetNelem(nelem)
   call a%Mesh%GetNdime(ndime)
   
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)
   do ielem = 1,nelem
   
      call a%Mesh%ElementLoad(ielem,e)  
      
      call e%elmdcg
      call e%elmlen
      
      call e%gather(1,eleden,a%densf(:,1))
      call e%gather(1,elepre,a%press(:,1))
      call e%gather(e%ndime,elevel,a%veloc(:,:,1))
      call e%gather(1,eletem,a%tempe(:,1))

      ! Compute the characteristic length chale
      call elmchl(e,1_ip,elevel,chale)

      error(ielem) = 0.0_rp

      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         call e%elmlen
         dvol = e%weigp(e%igaus)*e%detjm
   
         call e%interpg(e%ndime,e%elcod,gpcod)
         call e%interpg(1_ip,eleden,gpdens(1))
         call e%interpg(1_ip,elepre,gppress(1))
         call e%interpg(e%ndime,elevel,gpvele(:))
         call e%interpg(1_ip,eletem,gptempe(1))

         ! Exact solution
         call exacso%nsc_ComputeSolution(e%ndime,gpcod,a)
         call exacso%nsc_GetDensity(e%ndime,exden,exdeg)         
         call exacso%nsc_GetPressure(e%ndime,expre,exprg)         
         call exacso%nsc_GetVelocity(e%ndime,exvel,exveg) 
         call exacso%nsc_GetTemperature(e%ndime,extem,exteg)         

         !Scaled L_2 norm
         if (a%ErrorEstimatorTypeOfSubscales == 1) then

            ! Ideal Gas State Law 
            gpden = (gppress(1)+a%relpre)/((accph - accvh) * (gptempe(1)+a%reltem))
            !Advection velocity norm 
            call vecnor(gpvele,e%ndime,gpvno,2)
            !Sound speed
            call nsc_ComputeSoundSpeed(accph,accvh,gptempe(1)+a%reltem,gpspd)

            call a%ComputeNscompNstincVelocity(gpvno,gpspd,gpvst)
            call ComputeTau(e,gpden,acvis,gpvst,a%staco,chale,auxtm)
            call ComputeTauCDR(e,gpden*accvh,actco,gpden*a%grnor,gpvst,a%staco,chale,auxte) 

            timom(2) = auxtm
            timom(3) = auxte
            timom(1) = a%staco(4)*chale(2)*chale(2)/(gpden*auxtm*e%npol*e%npol)

            
            diffp = gppress(1)-expre
            difft = gptempe(1)-extem
            do idime=1,e%ndime
               diffvel(idime) = gpvele(idime)-exvel(idime)
            end do 

            call vecnor(diffvel(1:e%ndime),e%ndime,velnorm,2)

            error(ielem) = error(ielem) + (diffp**2)*dvol/(gpden*timom(1))
            error(ielem) = error(ielem) + (velnorm**2)*dvol/timom(2)
            error(ielem) = error(ielem) + (difft**2)*dvol/((gptempe(1)+a%reltem)*timom(3))
         !Entropy function

         elseif (a%ErrorEstimatorTypeOfSubscales == 2) then
   
            diffd = abs(gpdens(1)-exden)
            diffp = abs(gppress(1)-expre)

            acgamma = accph/accvh

            denominator = exden**acgamma 

            entropy = accvh*log((expre+a%relpre)/denominator) 

            denominator1 = gpdens(1)**acgamma 

            entropy1 = accvh*log((gppress(1)+a%relpre)/denominator1) 

            error(ielem) = error(ielem) + (entropy-entropy1)*(entropy-entropy1)*dvol/(entropy1*entropy1)

         endif

      enddo

   enddo
   
   
   call a%Memor%dealloc(e%mnode,eleden,'eleden','nsc_pr_ExactSolRefCriteria')
   call a%Memor%dealloc(e%mnode,elepre,'elepre','nsc_pr_ExactSolRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elevel,'elevel','nsc_pr_ExactSolRefCriteria')
   call a%Memor%dealloc(e%mnode,eletem,'eletem','nsc_pr_ExactSolRefCriteria')
   call a%Memor%dealloc(e%ndime,gpvele,'gpvele','nsc_pr_ExactSolRefCriteria')
   call a%Memor%dealloc(e%ndime,gpcod,'gpcod','nsc_pr_ExactSolRefCriteria')     
   

   ! Deallocate Exact Values
   call a%Memor%dealloc(e%ndime,exvel,'exvel','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,exprg,'exprg','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,e%ndime,exveg,'exveg','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,exteg,'exteg','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,diffvel,'diffvel','nsc_pr_ExactSolRefCriteria')     
   call a%Memor%dealloc(e%ndime,exdeg,'exdeg','nsc_pr_ExactSolRefCriteria')     

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
  
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','nsc_pr_ExactSolRefCriteria')
   
end subroutine   
