subroutine nsc_SubscalesRefCriteria(b,error,TotalEstimatedError)
   use typre
   use Mod_Memor
   use Mod_NSCompressibleImplicit
   use Mod_Element
   use Mod_nsc_BaseElmope
   use MPI
   implicit none
   class(NSCompressibleImplicitProblem), target :: b
   real(rp) :: error(:)
   
   integer(ip) :: ndime,ierr
   integer(ip) :: idime,npoinLocal
   real(rp) :: ortStressNorm, ortHFluxNorm,taue(2)
   real(rp) :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp) :: subscalesnorm,numerator,denominator,entropy,momnorm,numerator1,denominator1,entropy1
   real(rp) :: auxtm, auxte,invvis,invtco,invcp,invcv
   real(rp) :: gpdens(1), gpenerg(1)
   real(rp), allocatable :: eleden(:), elemom(:,:), eleene(:),gpmome(:)
   real(rp), allocatable :: elgpd(:,:), elgpm(:,:,:), elgpe(:,:)
   real(rp), allocatable :: gprjd(:), gprjm(:,:), gprje(:)
   real(rp), allocatable :: ortGradDen(:), ortGradMom(:,:), ortGradEne(:)
   real(rp), allocatable :: ortStress(:,:), ortHFlux(:)
   
   a=> b


   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elvel,'elvel','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,eleden,'eleden','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elemom,'elemom','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,eleene,'eleene','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gpmome,'gpmome','nsc_SubscalesRefCriteria')
   
   call a%Memor%alloc(e%ndime,grden,'grden','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,grmom,'grmom','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,grene,'grene','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elgpd,'elgpd','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgpm,'elgpm','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elgpe,'elgpe','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gprjd,'gprjd','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,gprjm,'gprjm','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gprje,'gprje','nsc_SubscalesRefCriteria')

   call a%Memor%alloc(e%ndime,ortGradDen,'ortGradDen','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,ortGradMom,'ortGradMom','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortGradEne,'ortGradEne','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,ortStress,'ortStress','nsc_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortHFlux,'ortHFlux','nsc_SubscalesRefCriteria')

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,actco,accph,accvh)
   do ielem = 1,nelem
     
      call a%Mesh%ElementLoad(ielem,e)  
      
      call e%elmdcg
      call e%elmlen
      
      call e%gather(1,eleden,a%densf(:,1))
      call e%gather(e%ndime,elemom,a%momen(:,:,1))
      call e%gather(1,eleene,a%energ(:,1))

      call nsc_ComputeElementVelocity(e,eleden,elemom,elvel)
      ! Compute the characteristic length chale
      call elmchl(e,1_ip,elvel,chale)

      !Scaled L_2 norm
      if (a%ErrorEstimatorTypeOfSubscales == 1) then
  
        
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
         elgpd(:,1:e%pnode) = a%dgrprj(:,e%lnods(1:e%pnode))
         do idime=1,e%ndime
            elgpm(idime,:,1:e%pnode) = a%grprj(idime,:,e%lnods(1:e%pnode))
         enddo
         elgpe(:,1:e%pnode) = a%grprj(e%ndime+1,:,e%lnods(1:e%pnode))

         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%elmder
            call e%elmlen
            dvol = e%weigp(e%igaus)*e%detjm
   
            call e%interpg(1_ip,eleden,gpdens(1))
            call e%interpg(e%ndime,elemom,gpmome(:))
            call e%interpg(1_ip,eleene,gpenerg(1))

            !Advection velocity norm 
            call vecnor(gpmome,e%ndime,gpmno,2)
            gptem = (gpenerg(1) - (gpmno*gpmno/(gpdens(1)*2.0_rp)))*invcv/gpdens(1)
            !Sound speed
            call nsc_ComputeSoundSpeed(accph,accvh,gptem,gpspd)

            call nsc_ComputeTau(e,gpspd,acvis/gpdens(1),actco*invcp/gpdens(1),gpmno/gpdens(1),a%staco,chale,timom)

            a%referencevelocity = gpspd + gpmno/gpdens(1)

            call vecnor(a%mosgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,subscalesnorm,2)

            !Interior subscales
            error(ielem) = error(ielem) + (a%cosgs(ielem)%a(1,e%igaus)**2)*(a%referencevelocity**2)*dvol/timom(1)
            error(ielem) = error(ielem) + (subscalesnorm**2)*dvol/timom(2)
            error(ielem) = error(ielem) + (a%ensgs(ielem)%a(1,e%igaus)**2)*dvol/((a%referencevelocity**2)*timom(3))
            !Subscales on the element boundaries

            ! Compute element variables gradient 
            call e%gradient(1_ip,eleden,grden)
            call e%gradient(e%ndime,elemom,grmom)
            call e%gradient(1_ip,eleene,grene)
            
            divmom = 0.0_rp

            gprjd(:) = matmul(elgpd(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
            ortGradDen = grden-gprjd
            do idime=1,e%ndime
               gprjm(idime,:) = matmul(elgpm(idime,:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
               ortGradMom(idime,:)   = grmom(idime,:)-gprjm(idime,:)
               divmom = divmom + ortGradMom(idime,idime) 
            enddo
            gprje(:) = matmul(elgpe(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
            ortGradEne = grene-gprje
            
            !compute stresses
            ortStress = acvis*(ortGradMom + transpose(ortGradMom))/gpdens(1)
            do idime=1,e%ndime
               ortStress(idime,1:e%ndime) = ortStress(idime,1:e%ndime) - acvis*(gpmome(idime)*ortGradDen(1:e%ndime) + gpmome(1:e%ndime)*ortGradDen(idime))/(gpdens(1)**2)
               ortStress(idime,idime) = ortStress(idime,idime) + 2.0_rp*acvis*divmom/(3.0_rp*gpdens(1))
               ortStress(idime,idime) = ortStress(idime,idime) - 2.0_rp*acvis*dot_product(gpmome,ortGradDen)/(3.0_rp*gpdens(1))
            enddo

            ortHFlux = -matmul(ortGradMom,gpmome)*actco*invcv/(gpdens(1)**2)
            ortHFlux = ortHflux + actco*ortGradEne*invcv/gpdens(1)  
            ortHFlux = ortHFlux - (gpmno**2)*actco*ortGradDen*invcv/(gpdens(1)**3)
            ortHFlux = ortHFlux - gpenerg(1)*actco*ortGradDen*invcv/(gpdens(1)**2)
            call matNormF(e%ndime,ortStress,ortStressNorm)
            call vecnor(ortHFlux,e%ndime,ortHFluxNorm,2)

            !compute taue
            taue(1) = e%hleng(1)*invvis*gpdens(1)/4
            taue(2) = e%hleng(1)*invtco*gpdens(1)*accph/4

            error(ielem) = error(ielem) + (ortStressNorm**2)*taue(1)*dvol/e%hleng(1)
            error(ielem) = error(ielem) + (ortHFluxNorm**2)*taue(2)*dvol/((a%referencevelocity**2)*e%hleng(1))

         enddo
      !Entropy function
      elseif (a%ErrorEstimatorTypeOfSubscales == 2) then

         do igaus = 1,e%pgaus
            e%igaus = igaus
            call e%elmder
            call e%elmlen
            dvol = e%weigp(e%igaus)*e%detjm

            call e%interpg(1_ip,eleden,gpdens(1))
            call e%interpg(e%ndime,elemom,gpmome(:))
            call e%interpg(1_ip,eleene,gpenerg(1))
   
            acgamma = accph/accvh

            call vecnor((gpmome(1:e%ndime)+a%mosgs(ielem)%a(1:e%ndime,1,e%igaus)),e%ndime,subscalesnorm,2)
            call vecnor((gpmome(1:e%ndime)),e%ndime,momnorm,2)

            numerator = (acgamma - 1.0_rp) * ((gpenerg(1)+a%ensgs(ielem)%a(1,e%igaus)) - subscalesnorm*subscalesnorm*0.5_rp/(gpdens(1)+a%cosgs(ielem)%a(1,e%igaus)))
            denominator = (gpdens(1)+a%cosgs(ielem)%a(1,e%igaus))**acgamma 

            entropy = accvh*log(numerator/denominator)

            numerator1 = (acgamma - 1.0_rp) * (gpenerg(1) - momnorm*momnorm*0.5_rp/gpdens(1))
            denominator1 = gpdens(1)**acgamma 
            entropy1 = accvh*log(numerator1/denominator1)

            error(ielem) = error(ielem) + (entropy-entropy1)*(entropy-entropy1)*dvol/(entropy1*entropy1)

         enddo

      endif

   enddo
   
   
   call a%Memor%dealloc(e%ndime,e%mnode,elvel,'elvel','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%mnode,eleden,'eleden','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elemom,'elemom','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%mnode,eleene,'eleene','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gpmome,'gpmome','nsc_SubscalesRefCriteria')

   call a%Memor%dealloc(e%ndime,grden,'grden','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,grmom,'grmom','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,grene,'grene','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elgpd,'elgpd','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgpm,'elgpm','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elgpe,'elgpe','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gprjd,'gprjd','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,gprjm,'gprjm','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gprje,'gprje','nsc_SubscalesRefCriteria')
                
   call a%Memor%dealloc(e%ndime,ortGradDen,'ortGradDen','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,ortGradMom,'ortGradMom','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortGradEne,'ortGradEne','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,ortStress,'ortStress','nsc_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortHFlux,'ortHFlux','nsc_SubscalesRefCriteria')

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
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','nsc_SubscalesRefCriteria')
   
end subroutine   
