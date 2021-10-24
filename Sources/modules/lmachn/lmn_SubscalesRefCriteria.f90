subroutine lmn_SubscalesRefCriteria(b,error,TotalEstimatedError)
   use MPI
   use typre
   use Mod_Element
   use Mod_LowMach
   use Mod_lmn_BaseElmope
   implicit none
   class(LowMachProblem), target  :: b
   real(rp) :: error(:)
   
   integer(ip) :: ndime,ierr
   integer(ip) :: idime,npoinLocal
   real(rp) :: ortStressNorm, ortHFluxNorm,taue(2),subscalesnorm
   real(rp) :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp) :: auxtm, auxte,invvis,invtco,acgas
   real(rp) :: gppress(1), gptempe(1)
   real(rp), allocatable :: elepre(:), elevel(:,:), eletem(:),gpvele(:)
   real(rp), allocatable :: elgpv(:,:,:), elgpt(:,:)
   real(rp), allocatable :: gprjv(:,:), gprjt(:)
   real(rp), allocatable :: ortGradVel(:,:), ortGradTempe(:)
   real(rp), allocatable :: ortStress(:,:), ortHFlux(:)
   integer(ip), save     :: ipass1=0,ipass2=0

   a => b

   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,elepre,'elepre','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elevel,'elevel','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,eletem,'eletem','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gpvele,'gpvele','lmn_SubscalesRefCriteria')
   
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(1_ip,e%ndime,grtem,'grtem','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elgpt,'elgpt','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,gprjv,'gprjv','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gprjt,'gprjt','lmn_SubscalesRefCriteria')

   call a%Memor%alloc(e%ndime,e%ndime,ortGradVel,'ortGradVel','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortGradTempe,'ortGradTempe','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,ortStress,'ortStress','lmn_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortHFlux,'ortHFlux','lmn_SubscalesRefCriteria')

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   
   !Physical Parameters
   call a%GetPhysicalParameters(acvis,acsph,actco,acgas)
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
     
      if(acvis>zelmn) then
         invvis = 1.0_rp/acvis
      else
         invvis = 0.0_rp
      end if
      if(actco>zelmn) then
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
         acden = a%pther(1)/acgas/gptempe(1)
         !Advection velocity norm 
         call vecnor(gpvele,e%ndime,gpvno,2)

         call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timom)
         call ComputeTau(e,acden,actco,gpvno,a%staco,chale,tiene)

         ticon = chale(2)*chale(2)/(a%staco(1)*acden*auxtm*e%npol*e%npol)

         call vecnor(a%vesgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,subscalesnorm,2)

         !Interior subscales
         error(ielem) = error(ielem) + (a%prsgs(ielem)%a(e%igaus)**2)*dvol/(acden*ticon)
         error(ielem) = error(ielem) + (subscalesnorm**2)*dvol/timom
         error(ielem) = error(ielem) + (a%tesgs(ielem)%a(1,e%igaus)**2)*dvol/(gptempe(1)*tiene)
         !Subscales on the element boundaries

         ! Compute element variables gradient 
         call e%gradient(e%ndime,elevel,grvel)
         call e%gradient(1_ip,eletem,grtem(1,:))
         
         divvel = 0.0_rp

         do idime=1,e%ndime
            gprjv(idime,:) = matmul(elgpv(idime,:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
            ortGradVel(idime,:)   = grvel(idime,:)-gprjv(idime,:)
            divvel = divvel + ortGradVel(idime,idime) 
         enddo
         gprjt(:) = matmul(elgpt(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
         ortGradTempe = grtem(1,:)-gprjt
         
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
         error(ielem) = error(ielem) + 0.1*(ortHFluxNorm**2)*taue(2)*dvol/(gptempe(1)*e%hleng(1))
      enddo
   enddo
   
   call a%Memor%dealloc(e%mnode,elepre,'elepre','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elevel,'elevel','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%mnode,eletem,'eletem','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gpvele,'gpvele','lmn_SubscalesRefCriteria')
   
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(1_ip,e%ndime,grtem,'grtem','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elgpt,'elgpt','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,gprjv,'gprjv','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gprjt,'gprjt','lmn_SubscalesRefCriteria')

   call a%Memor%dealloc(e%ndime,e%ndime,ortGradVel,'ortGradVel','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortGradTempe,'ortGradTempe','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,ortStress,'ortStress','lmn_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortHFlux,'ortHFlux','lmn_SubscalesRefCriteria')

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
   if (a%MPIrank == a%MPIroot .and. a%kfl_adap) then
      if(ipass1==0) then
         ipass1=1
         write(a%lun_adapt,400)
      end if
      write(a%lun_adapt,401) a%istep, a%ctime, nelem
      if (a%kfl_flush == 1) call flush(a%lun_adapt)
   endif
   if (a%MPIrank == a%MPIroot .and. a%kfl_error) then
      if(ipass2==0) then
         ipass2=1
         write(a%lun_ersgs,200)
      end if
      write(a%lun_ersgs,201) a%istep, a%ctime, TotalEstimatedError, maxval(error), minval(error)
      if (a%kfl_flush == 1) call flush(a%lun_error)
   endif
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','lmn_SubscalesRefCriteria')
   
   !Formats
200 format('$ ','     Time','     Current','       Error','        Error','        Error',/,&
         & '$ ','     step','        time','         SGS','          max','          min')
201 format(i9,1x,4(2x,e12.6))
400 format('$ ','     Time','     Current','      Number of',/,&
         & '$ ','     step','        time','       elements')
401 format(i9,1x,e12.6,1x,i9)
   
end subroutine   
