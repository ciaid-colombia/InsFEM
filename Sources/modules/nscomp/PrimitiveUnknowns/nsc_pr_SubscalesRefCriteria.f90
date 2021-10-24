subroutine nsc_pr_SubscalesRefCriteria(b,error,TotalEstimatedError)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_NSCompressiblePrimitive
   use Mod_nsc_pr_BaseElmope
   use MPI
   implicit none
   class(NSCompressiblePrimitiveProblem), target  :: b
   real(rp) :: error(:)
   
   integer(ip) :: ndime,ierr
   integer(ip) :: idime,npoinLocal
   real(rp) :: ortStressNorm, ortHFluxNorm,taue(2)
   real(rp) :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp) :: subscalesnorm,denominator,entropy,acgamma,denominator1,entropy1
   real(rp) :: auxtm, auxte,invvis,invtco
   real(rp) :: gppress(1), gptempe(1)
   real(rp), allocatable :: elepre(:), elevel(:,:), eletem(:),gpvele(:)
   real(rp), allocatable :: elgpv(:,:,:), elgpt(:,:)
   real(rp), allocatable :: gprjv(:,:), gprjt(:)
   real(rp), allocatable :: ortGradVel(:,:), ortGradTempe(:)
   real(rp), allocatable :: ortStress(:,:), ortHFlux(:)


   a => b

   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,elepre,'elepre','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elevel,'elevel','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,eletem,'eletem','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gpvele,'gpvele','nsc_pr_SubscalesRefCriteria')
   
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,grtem,'grtem','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elgpt,'elgpt','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,gprjv,'gprjv','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gprjt,'gprjt','nsc_pr_SubscalesRefCriteria')

   call a%Memor%alloc(e%ndime,e%ndime,ortGradVel,'ortGradVel','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortGradTempe,'ortGradTempe','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,ortStress,'ortStress','nsc_pr_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortHFlux,'ortHFlux','nsc_pr_SubscalesRefCriteria')

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)
   do ielem = 1,nelem
     
      call a%Mesh%ElementLoad(ielem,e)  
      
      call e%elmdcg
      call e%elmlen
      
      call e%gather(1,elepre,a%press(:,1))
      call e%gather(e%ndime,elevel,a%veloc(:,:,1))
      call e%gather(1,eletem,a%tempe(:,1))

      ! Compute the characteristic length chale
      call elmchl(e,1_ip,elevel,chale)

      error(ielem) = 0.0_rp

      !Scaled L_2 norm
      if (a%ErrorEstimatorTypeOfSubscales == 1) then
        
         if(acvis>zensi) then
            invvis = 1.0_rp/acvis
         else
            invvis = 0.0_rp
         end if
         if(actco>zensi) then
            invtco = 1.0_rp/actco
         else
            invtco = 0.0_rp
         end if
         
         !Gradient projection for stress calculation 
         do idime=1,e%ndime
            elgpv(idime,:,1:e%pnode) = a%grprj(idime,:,e%lnods(1:e%pnode))
         enddo
         elgpt(:,1:e%pnode) = a%grprj(e%ndime+1,:,e%lnods(1:e%pnode))

         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%elmder
            call e%elmlen
            dvol = e%weigp(e%igaus)*e%detjm
   
            call e%interpg(1_ip,elepre,gppress(1))
            call e%interpg(e%ndime,elevel,gpvele(:))
            call e%interpg(1_ip,eletem,gptempe(1))

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

            call vecnor(a%mosgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,subscalesnorm,2)

            !Interior subscales
            error(ielem) = error(ielem) + (a%cosgs(ielem)%a(1,e%igaus)**2)*dvol/(gpden*timom(1))
            error(ielem) = error(ielem) + (subscalesnorm**2)*dvol/timom(2)
            error(ielem) = error(ielem) + (a%ensgs(ielem)%a(1,e%igaus)**2)*dvol/((gptempe(1)+a%reltem)*timom(3))   
            !Subscales on the element boundaries

            ! Compute element variables gradient 
            call e%gradient(e%ndime,elevel,grvel)
            call e%gradient(1_ip,eletem,grtem)
            
            divvel = 0.0_rp

            do idime=1,e%ndime
               gprjv(idime,:) = matmul(elgpv(idime,:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
               ortGradVel(idime,:)   = grvel(idime,:)-gprjv(idime,:)
               divvel = divvel + ortGradVel(idime,idime) 
            enddo
            gprjt(:) = matmul(elgpt(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
            ortGradTempe = grtem-gprjt
            
            !compute stresses
            ortStress = acvis*(ortGradVel + transpose(ortGradVel))
            do idime=1,e%ndime
               ortStress(idime,idime) = ortStress(idime,idime) + 2.0_rp*acvis*divvel/3.0_rp
            enddo
            ortHFlux = actco*ortGradTempe 
            call matNormF(e%ndime,ortStress,ortStressNorm)
            call vecnor(ortHFlux,e%ndime,ortHFluxNorm,2)

            !compute taue
            taue(1) = e%hleng(1)*invvis/4
            taue(2) = e%hleng(1)*invtco/4

            error(ielem) = error(ielem) + 0.1*(ortStressNorm**2)*taue(1)*dvol/e%hleng(1)
            error(ielem) = error(ielem) + 0.1*(ortHFluxNorm**2)*taue(2)*dvol/((gptempe(1)+a%reltem)*e%hleng(1))
         enddo

      !Entropy function
      elseif (a%ErrorEstimatorTypeOfSubscales == 2) then

         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%elmder
            call e%elmlen
            dvol = e%weigp(e%igaus)*e%detjm

            call e%interpg(1_ip,elepre,gppress(1))
            call e%interpg(1_ip,eletem,gptempe(1))

            acgamma = accph/accvh

            denominator = ((a%relpre+gppress(1)+a%cosgs(ielem)%a(1,e%igaus))/((accph-accvh)*(a%reltem+gptempe(1)+a%ensgs(ielem)%a(1,e%igaus))))**acgamma 
            entropy = accvh*log((a%relpre+gppress(1)+a%cosgs(ielem)%a(1,e%igaus))/denominator)

            denominator1 = ((a%relpre+gppress(1))/((accph-accvh)*(a%reltem+gptempe(1))))**acgamma 
            entropy1 = accvh*log((a%relpre+gppress(1))/denominator1)

            error(ielem) = error(ielem) + (entropy-entropy1)*(entropy-entropy1)*dvol/(entropy1*entropy1)

         enddo

      endif

   enddo
   
   
   call a%Memor%dealloc(e%mnode,elepre,'elepre','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elevel,'elevel','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%mnode,eletem,'eletem','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gpvele,'gpvele','nsc_pr_SubscalesRefCriteria')
   
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,grtem,'grtem','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elgpt,'elgpt','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,gprjv,'gprjv','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gprjt,'gprjt','nsc_pr_SubscalesRefCriteria')

   call a%Memor%dealloc(e%ndime,e%ndime,ortGradVel,'ortGradVel','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortGradTempe,'ortGradTempe','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,ortStress,'ortStress','nsc_pr_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortHFlux,'ortHFlux','nsc_pr_SubscalesRefCriteria')

   !Compute the total error
   TotalEstimatedError0 = 0.0_rp
   call a%Mesh%GetNpoinLocal(npoinLocal)
   do ielem = 1,nelem
      call a%Mesh%ElementLoad(ielem,e)
      weightfactor = 1.0-real(count(e%lnods(1:e%pnode)>npoinLocal))/real(e%pnode)
   
      TotalEstimatedError0 = TotalEstimatedError0 + error(ielem)*weightfactor

   enddo
   call MPI_ALLREDUCE(TotalEstimatedError0,TotalEstimatedError,1,MPI_REAL8,MPI_SUM,a%MPIcomm,ierr)
   
   TotalEstimatedError = sqrt(TotalEstimatedError)
   if (a%MPIrank == a%MPIroot .and. a%kfl_adap) write(a%lun_adapt,*) ' Estimated error: ',TotalEstimatedError
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','nsc_pr_SubscalesRefCriteria')
   
end subroutine   
