subroutine nsi_SubscalesRefCriteria(b,error,TotalEstimatedError)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_NavierStokes
   use Mod_nsm_BaseElmope
   use MPI
   implicit none
   class(NavierStokesProblem), target  :: b
   real(rp) :: error(:)
   
   integer(ip) :: npoinLocal,ierr
   real(rp)    :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp)    :: ortStressNorm,ortHFluxNorm, taue1,taue2, subscalesnorm,auxtm,invvis,invtco
   real(rp), allocatable :: elevel(:,:), gpvele(:), grtem(:,:)
   real(rp), allocatable :: elgpv(:,:,:), elgpt(:,:), gprjv(:,:), gprjt(:)
   real(rp), allocatable :: ortGradVel(:,:), ortGradTempe(:)
   real(rp), allocatable :: ortStress(:,:), ortHFlux(:)
   real(rp)    :: gptempe(1),tiene,actco
   real(rp), allocatable :: eletem(:)
   integer(ip), save     :: ipass1=0,ipass2=0

   a => b

   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elevel,'elevel','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gpvele,'gpvele','nsm_SubscalesRefCriteria')
   
   call a%Memor%alloc(e%ndime,e%ndime,grvel,'grvel','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,gprjv,'gprjv','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(1_ip,e%ndime,grtem,'grtem','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elgpt,'elgpt','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gprjt,'gprjt','nsm_SubscalesRefCriteria')

   call a%Memor%alloc(e%ndime,e%ndime,ortGradVel,'ortGradVel','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortGradTempe,'ortGradTempe','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%ndime,ortStress,'ortStress','nsm_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortHFlux,'ortHFlux','nsm_SubscalesRefCriteria')

   if (a%kfl_cotem == 1) then
      call a%Memor%alloc(e%mnode,eletem,'eletem','nsm_SubscalesRefCriteria')
   end if

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   
   !Physical Parameters
   call a%GetPhysicalParameters(1,acden,acvis)
   do ielem = 1,nelem
     
      call a%Mesh%ElementLoad(ielem,e)  
      
      call e%elmdcg
      call e%elmlen
      
      call e%gather(e%ndime,elevel,a%veloc(:,:,1))
      if (a%kfl_cotem == 1) call e%gather(1,eletem,a%tempe)

      ! Compute the characteristic length chale
      call elmchl(e,a%kfl_advec,elevel,chale,a%kfl_hdifumin)

      error(ielem) = 0.0_rp

      !Scaled L_2 norm
        
      if(acvis>zensi) then
         invvis = 1.0_rp/acvis
      else
         invvis = 0.0_rp
      end if
      
      !Gradient projection for stress calculation 
      do idime=1,e%ndime
         elgpv(idime,:,1:e%pnode) = a%grprj(idime,:,e%lnods(1:e%pnode))
      enddo

      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         call e%elmlen
         dvol = e%weigp(e%igaus)*e%detjm
   
         call e%interpg(e%ndime,elevel,gpvele(:))
         if (a%kfl_cotem == 1) call e%interpg(1_ip,eletem,gptempe(1))

         !Advection velocity norm 
         call vecnor(gpvele,e%ndime,gpvno,2)

         call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,auxtm)
         if (a%kfl_cotem == 1) call ComputeTauCDR(e,acden,0.0035196_rp,0.0_rp,gpvno,a%staco,chale,tiene)

         timom = auxtm
         tidiv = a%staco(4)*chale(2)*chale(2)/(acden*auxtm*e%npol*e%npol)

         call vecnor(a%vesgs(ielem)%a(1:e%ndime,1,e%igaus),e%ndime,subscalesnorm,2)

         !Interior subscales
         error(ielem) = error(ielem) + (a%prsgs(ielem)%a(e%igaus)**2)*dvol/(acden*tidiv)
         error(ielem) = error(ielem) + (subscalesnorm**2)*dvol/timom
         if (a%kfl_cotem == 1) error(ielem) = error(ielem) + (a%tesgs(ielem)%a(1,e%igaus)**2)*dvol/(gptempe(1)*tiene)
         
         !Subscales on the element boundaries

         ! Compute element variables gradient 
         call e%gradient(e%ndime,elevel,grvel)
         
         divvel = 0.0_rp

         do idime=1,e%ndime
            gprjv(idime,:) = matmul(elgpv(idime,:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
            ortGradVel(idime,:)   = grvel(idime,:)-gprjv(idime,:)
            divvel = divvel + ortGradVel(idime,idime) 
         enddo

         !compute stresses
         ortStress = acvis*(ortGradVel + transpose(ortGradVel))
         do idime=1,e%ndime
            ortStress(idime,idime) = ortStress(idime,idime) + 2.0_rp*acvis*divvel/3.0_rp
         enddo
         call matNormF(e%ndime,ortStress,ortStressNorm)

         !compute taue
         taue1 = e%hleng(1)*invvis/4

         error(ielem) = error(ielem) + (ortStressNorm**2)*taue1*dvol/e%hleng(1)

         !compute T
         if (a%kfl_cotem == 1) then
            actco=3.59
            gprjt(:) = matmul(elgpt(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
            ortGradTempe = grtem(1,:)-gprjt
            ortHFlux = actco*ortGradTempe 
            call vecnor(ortHFlux,e%ndime,ortHFluxNorm,2)
            taue2 = e%hleng(1)*invtco/4
            error(ielem) = error(ielem) + (ortHFluxNorm**2)*taue2*dvol/(gptempe(1)*e%hleng(1))
         end if
      enddo

   enddo
   
   
   call a%Memor%dealloc(e%ndime,e%mnode,elevel,'elevel','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gpvele,'gpvele','nsm_SubscalesRefCriteria')
   
   call a%Memor%dealloc(e%ndime,e%ndime,grvel,'grvel','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,elgpv,'elgpv','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,gprjv,'gprjv','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(1_ip,e%ndime,grtem,'grtem','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elgpt,'elgpt','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gprjt,'gprjt','nsm_SubscalesRefCriteria')

   call a%Memor%dealloc(e%ndime,e%ndime,ortGradVel,'ortGradVel','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortGradTempe,'ortGradTempe','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%ndime,ortStress,'ortStress','nsm_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortHFlux,'ortHFlux','nsm_SubscalesRefCriteria')

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
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','nsm_SubscalesRefCriteria')

   !Formats
200 format('$ ','     Time','     Current','       Error','        Error','        Error',/,&
         & '$ ','     step','        time','         SGS','          max','          min')
201 format(i9,1x,4(2x,e12.6))
400 format('$ ','     Time','     Current','      Number of',/,&
         & '$ ','     step','        time','       elements')
401 format(i9,1x,e12.6,1x,i9)
   
end subroutine   
