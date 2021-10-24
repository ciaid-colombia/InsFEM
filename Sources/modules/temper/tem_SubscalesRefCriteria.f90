subroutine tem_SubscalesRefCriteria(b,error,TotalEstimatedError)
   use typre
   use Mod_Memor
   use Mod_Element
   use Mod_Temperature
   use Mod_tem_BaseElmope
   use MPI
   implicit none
   class(TemperatureProblem), target  :: b
   real(rp) :: error(:)
   
   integer(ip) :: ndime,ierr
   integer(ip) :: idime,npoinLocal
   real(rp) :: TotalEstimatedError,TotalEstimatedError0,weightfactor
   real(rp) :: invtco,gptempe(1),ortHFluxNorm,taue
   real(rp), allocatable :: elevel(:,:), eletem(:),gpvele(:)
   real(rp), allocatable :: elgpt(:,:),gprjt(:),grtem(:)
   real(rp), allocatable :: ortGradTempe(:),ortHFlux(:)
   integer(ip), save     :: ipass1=0,ipass2=0

   a => b

   !Subscales as error indicator
   error = 0.0_rp
   
   if (a%istep == 0) return
   
   call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','tem_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elevel,'elevel','tem_SubscalesRefCriteria')
   call a%Memor%alloc(e%mnode,eletem,'eletem','tem_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gpvele,'gpvele','tem_SubscalesRefCriteria')
   
   call a%Memor%alloc(e%ndime,grtem,'grtem','tem_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,e%mnode,elgpt,'elgpt','tem_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,gprjt,'gprjt','tem_SubscalesRefCriteria')

   call a%Memor%alloc(e%ndime,ortGradTempe,'ortGradTempe','tem_SubscalesRefCriteria')
   call a%Memor%alloc(e%ndime,ortHFlux,'ortHFlux','tem_SubscalesRefCriteria')

   call a%Mesh%GetNdime(ndime)
   call a%Mesh%GetNelem(nelem)
   

   do ielem = 1,nelem

      !Physical Parameters
      call a%GetPhysicalParameters(ielem,acden,acsph,actco,acrea,acsou)
      acvis = actco/acsph
      acrcp = acrea/acsph
     
      call a%Mesh%ElementLoad(ielem,e)  
      
      call e%elmdcg
      call e%elmlen
      
      call e%gather(e%ndime,elevel,a%veloc(:,:))
      call e%gather(1,eletem,a%tempe(:,1))

      ! Compute the characteristic length chale
      call elmchl(e,1_ip,elevel,chale)

      error(ielem) = 0.0_rp

      !Scaled L_2 norm
        
      if(actco>zetem) then
         invtco = 1.0_rp/actco
      else
         invtco = 0.0_rp
      end if
      
      !Gradient projection for stress calculation 
      elgpt(:,1:e%pnode) = a%grprj(:,e%lnods(1:e%pnode))

      do igaus = 1,e%pgaus
         e%igaus = igaus
         call e%elmder
         call e%elmlen
         dvol = e%weigp(e%igaus)*e%detjm

         call e%interpg(e%ndime,elevel,gpvele(:))
         call e%interpg(1_ip,eletem,gptempe(1))

         !Advection velocity norm 
         call vecnor(gpvele,e%ndime,gpvno,2)

         call ComputeTauCDR(e,acden,acvis,acrcp,gpvno,a%staco,chale,timom)

         !Interior subscales
         error(ielem) = error(ielem) + (a%tesgs(ielem)%a(1,e%igaus)**2)*dvol/(gptempe(1)*timom)   
         !Subscales on the element boundaries

         ! Compute element variables gradient 
         call e%gradient(1_ip,eletem,grtem)
         
         gprjt(:) = matmul(elgpt(:,1:e%pnode),e%shape(1:e%pnode,e%igaus))
         ortGradTempe = grtem-gprjt
         
         ortHFlux = actco*ortGradTempe 
         call vecnor(ortHFlux,e%ndime,ortHFluxNorm,2)

         !compute taue
         taue = e%hleng(1)*invtco/4

         error(ielem) = error(ielem) + 0.1*(ortHFluxNorm**2)*taue*dvol/(gptempe(1)*e%hleng(1))
      enddo

   enddo
   
   
   call a%Memor%dealloc(e%ndime,e%mnode,elevel,'elevel','tem_SubscalesRefCriteria')
   call a%Memor%dealloc(e%mnode,eletem,'eletem','tem_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gpvele,'gpvele','tem_SubscalesRefCriteria')
   
   call a%Memor%dealloc(e%ndime,grtem,'grtem','tem_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,e%mnode,elgpt,'elgpt','tem_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,gprjt,'gprjt','tem_SubscalesRefCriteria')

   call a%Memor%dealloc(e%ndime,ortGradTempe,'ortGradTempe','tem_SubscalesRefCriteria')
   call a%Memor%dealloc(e%ndime,ortHFlux,'ortHFlux','tem_SubscalesRefCriteria')

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
      if (a%kfl_flush == 1) call flush(a%lun_error)
   endif
   if (a%MPIrank == a%MPIroot .and. a%kfl_error) then
      if(ipass2==0) then
         ipass2=1
         write(a%lun_ersgs,200)
      end if
      write(a%lun_ersgs,201) a%istep, a%ctime, TotalEstimatedError, maxval(error), minval(error)
      if (a%kfl_flush == 1) call flush(a%lun_error)
   endif
   
   call a%Mesh%ElementdeAlloc(e,a%Memor,'DefaultRule','tem_SubscalesRefCriteria')
   
   !Formats
200 format('$ ','     Time','     Current','       Error','        Error','        Error',/,&
         & '$ ','     step','        time','         SGS','          max','          min')
201 format(i9,1x,4(2x,e12.6))
400 format('$ ','     Time','     Current','      Number of',/,&
         & '$ ','     step','        time','       elements')
401 format(i9,1x,e12.6,1x,i9)
   
end subroutine   
