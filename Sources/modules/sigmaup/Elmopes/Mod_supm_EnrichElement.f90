module Mod_supm_EnrichElement
 use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersEnrichElement
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersEnrichElement(itask)
      integer(ip) :: itask
      !procedure() :: NULL()
      
      select case (itask)   
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            if(a%kfl_EnrichElem == 1)then
               call PrependProcedure(ProcHook_PreAllocate,MatricesReallocTF)
               call ConcatenateProcedures(ProcHook_PreGauss,ComputeEnrichNodeTF)
               call PrependProcedure(ProcHook_InGauss,RecoverPnodeTF)               
               call ConcatenateProcedures(ProcHook_InGauss,ExtraNodeDerivativesTF)
               call ConcatenateProcedures(ProcHook_PostInterpolate,ChangePnodeTF)
               call PrependProcedure( ProcPointer%PreAssembly_sup,RecoverPnodeTF)                
               call PrependProcedure( ProcPointer%PreAssembly_sup,CondensationMatrixTF)
            end if        
         end if
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select   
      
   end subroutine  
   
   

   subroutine MatricesReallocTF
      implicit none
      e%mnode = e%mnode+1
      call a%Memor%prealloc(e%ndime,e%mnode,e%cartd,'e%cartd','MatricesRealloc') 
   end subroutine   
   
   subroutine ComputeEnrichNodeTF
      implicit none
      integer(ip) :: inode,ipoin,pStatus
      !Array to modify the shapes functions in the enriched node
      real(rp)    :: kshape(e%pnode)
      integer(ip) :: nodeStatus(e%pnode),ninters
      integer(ip)  :: elemStatus       
      
      call a%CutMesh%GetElementType(ielem,elemStatus)      
      if(elemStatus ==0)then         
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         call a%CutMesh%GetEnrichBuble(ielem,e%pnode,kshape)    
         do inode=1,e%pnode
            ipoin=e%lnods(inode)
            call a%CutMesh%GetPointType(ipoin,pStatus)
            nodeStatus(inode)= pStatus      
         end do            
         call e%SetEnrichedElement(a%Memor,e%pnode+1,kshape,ngauss_minus,ngauss_plus,nodeStatus)        
      end if
   end subroutine      
   
   subroutine ExtraNodeDerivativesTF
      implicit none
      integer(ip) :: inode,ipoin,pStatus
      !Array to modify the shapes functions in the enriched node
      real(rp)    :: kshape(e%pnode)
      integer(ip) :: nodeStatus(e%pnode),ninters,idime
      integer(ip)  :: elemStatus,auxinit
      integer(ip)  :: gauss_points(2,e%pnode)     
      
      call a%CutMesh%GetElementType(ielem,elemStatus)      
      if(elemStatus ==0)then
         call a%CutMesh%GetNgaussSide(ielem,ngauss_plus,ngauss_minus)
         call a%CutMesh%GetEnrichBuble(ielem,e%pnode,kshape) 
         !To work later with the original Cartds
         !The basic cartd's
         e%Particularcartd(1:e%ndime,e%pnode+1) = 0.0_rp !the additional node is defined later
         e%Particularcartd(1:e%ndime,1:e%pnode) = e%cartd(1:e%ndime,1:e%pnode)
         !Enriched cartd construction
         if(e%igaus<=ngauss_minus)then
            do inode = 1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,pStatus)
               if (pStatus > 0) then !The construction is with the other side            
                  do idime=1,e%ndime
                     e%Particularcartd(idime,e%pnode+1) = e%Particularcartd(idime,e%pnode+1) + &
                        e%cartd(idime,inode)*kshape(inode)
                  end do
               end if
            enddo  
         elseif(e%igaus>ngauss_minus)then
            do inode = 1,e%pnode
               ipoin=e%lnods(inode)
               call a%CutMesh%GetPointType(ipoin,pStatus)
               if(pStatus < 0) then       
                  do idime=1,e%ndime
                     e%Particularcartd(idime,e%pnode+1) = e%Particularcartd(idime,e%pnode+1) + &
                        e%cartd(idime,inode)*kshape(inode) 
                  end do
               end if
               
            enddo           
         end if
         !Finally we redefine the reallocated cartd
         e%cartd(1:e%ndime,e%pnode+1)=e%Particularcartd(1:e%ndime,e%pnode+1)     
      end if      
   end subroutine   

   subroutine ChangePnodeTF
      implicit none
      integer(ip)  :: elemStatus       
      call a%CutMesh%GetElementType(ielem,elemStatus)      
      if(elemStatus ==0)then
         if(ipnode0==0)then
            e%pnode=e%pnode +1
            ipnode0 = -1           
         end if
      end if
   
   end subroutine
   
   subroutine RecoverPnodeTF
      implicit none
      if(ipnode0==-1)then  !always if we change pnode
         e%pnode=e%pnode-1 !the original value    
         ipnode0=0
      end if
   end subroutine      
      
   subroutine CondensationMatrixTF
      implicit none
      integer(ip) :: idofn,jdofn,auxdimi,auxdimj,inode,jnode,idime,jdime
      !Auxiliar Arrays
      real(rp) :: A22(a%ndofn,a%ndofn),A22_Inv(a%ndofn,a%ndofn),deter
      real(rp) :: A12(a%ndofn*(e%pnode-1),a%ndofn)  
      real(rp) :: A21(a%ndofn,a%ndofn*(e%pnode-1))
      real(rp) :: RHS2(a%ndofn)
      real(rp) :: B11(a%ndofn*(e%pnode-1),a%ndofn)
      real(rp) :: C11(a%ndofn*(e%pnode-1),a%ndofn*(e%pnode-1))
      real(rp) :: R11(a%ndofn*(e%pnode-1))
      !Auxiliars elmats
      real(rp) :: Celmat(a%ndofn,(e%pnode-1),a%ndofn,(e%pnode-1))
      real(rp) :: Celrhs(a%ndofn,(e%pnode-1))
      integer(ip)  :: elemStatus 
      integer(ip)  :: ntens
      
      !--------------------------------------------------!
      !--------------------Notation----------------------!
      !                 K     U  =  RHS                  !
      !                                                  !
      !              A11 A12  U1    RHS1                 !
      !              A21 A22  U2    RHS2                 !
      !-------------------Condensation-------------------!
      !                                                  !
      ![A11-A12*(A22^-1)*A21] U1= RHS1-A12*(A22^-1)*RHS2 !
      !                                                  !
      !-------------------Nomenclature-------------------!
      ! B11=A12*(A22^-1)                                 !
      ! C11=B11*A21                                      !
      ! R11=B11*RHS2                                     !      
      !--------------------------------------------------!    
      ntens=(e%ndime-1)*(e%ndime-1)+2
      
      call a%CutMesh%GetElementType(ielem,elemStatus)      
      if(elemStatus ==0)then           
         !Type of enrichment
         if(a%kfl_EnrichVelocity==0)then
            do inode=1,e%pnode
               elmat(1:a%ndofn,inode,1:(e%ndime+ntens),e%pnode) = 0.0_rp
               elmat(1:(e%ndime+ntens),e%pnode,1:a%ndofn,inode) = 0.0_rp
            end do
            do idime=1,(e%ndime+ntens)
               elmat(idime,e%pnode,idime,e%pnode) = 1.0_rp               
            end do            
            elrhs(1:(e%ndime+ntens),e%pnode) = 0.0_rp           
         elseif(a%kfl_EnrichPressure==0)then
            do inode=1,e%pnode
               elmat(1:a%ndofn,inode,a%ndofn,e%pnode) = 0.0_rp
               elmat(a%ndofn,e%pnode,1:a%ndofn,inode) = 0.0_rp
            end do
            elmat(a%ndofn,e%pnode,a%ndofn,e%pnode) = 1.0_rp                       
            elrhs(a%ndofn,e%pnode) = 0.0_rp           
         endif
         
         !Now we can condense the enriched values         
         !Definition of the matrices in an array to use MATMUL
         
         !Initialization
         A22=0.0_rp
         A22_Inv=0.0_rp
         A21=0.0_rp
         A12=0.0_rp
         RHS2=0.0_rp 
         
         auxdimi=0
         auxdimj=0
         
         RHS2 = elrhs(:,e%pnode)
         A22 = elmat(:,e%pnode,:,e%pnode)  
         do inode = 1,e%pnode-1
            auxdimi = (inode-1)*a%ndofn
            A21(:,auxdimi+1:auxdimi+a%ndofn) = elmat(:,e%pnode,:,inode)
            A12(auxdimi+1:auxdimi+a%ndofn,:) = elmat(:,inode,:,e%pnode)
         enddo         
         
         !Now we can work with the matrix tools
         if(A22(a%ndofn,a%ndofn)==0.0_rp .and. e%ndime==3)then
         
         else   
            call invmtx(A22,A22_Inv,deter,a%ndofn)
         end if
         B11 = matmul(A12,A22_Inv)
         C11 = matmul(B11,A21)
         R11 = matmul(B11,RHS2)
         
         !Recovery of the standar elemental matrix form
         Celmat=0.0_rp
         Celrhs=0.0_rp      
         auxdimi=0
         auxdimj=0
         
         do inode=1,(e%pnode-1)
            auxdimi=(inode-1)*(a%ndofn)
            do jnode=1,(e%pnode-1)
            auxdimj=(jnode-1)*(a%ndofn)
               Celmat(:,inode,:,jnode) = C11(auxdimi+1:auxdimi+a%ndofn,auxdimj+1:auxdimj+a%ndofn)            
            end do
            Celrhs(:,inode) = R11(auxdimi+1:auxdimi+a%ndofn)            
         end do
                                 
         !Now we add the Condensed matrixes to elmat and to elrhs 
         elmat(1:a%ndofn,1:(e%pnode-1),1:a%ndofn,1:(e%pnode-1)) = elmat(1:a%ndofn,1:(e%pnode-1),1:a%ndofn,1:(e%pnode-1)) &
            - Celmat(1:a%ndofn,1:(e%pnode-1),1:a%ndofn,1:(e%pnode-1))
            
         elrhs(1:a%ndofn,1:(e%pnode-1)) = elrhs(1:a%ndofn,1:(e%pnode-1)) - Celrhs(1:a%ndofn,1:(e%pnode-1))         
      end if
   end subroutine
   
end module