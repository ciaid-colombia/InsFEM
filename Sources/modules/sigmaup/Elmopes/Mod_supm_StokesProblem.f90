module Mod_supm_StokesProblem
 use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersStokesProblem
   integer(ip), allocatable :: kfl_IsSet
   
contains
  
   subroutine SetPointersStokesProblem(itask)
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
            if(a%kfl_fsurf==1)then
               if(a%kfl_fsurfLapla==1)then
                  call PrependProcedure( ProcPointer%PreAssembly_sup,ElmatsToLapla)
               end if
               if(a%kfl_fsurfLapla==2)then
                  call ConcatenateProcedures(ProcHook_Initializations,AllocStokes) 
                  call ConcatenateProcedures(ProcHook_Finalizations,DeallocStokes)
                  call ConcatenateProcedures(ProcHook_PreGauss,InitElmuq)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,GalerkinElmuq)
                  call ConcatenateProcedures( ProcPointer%PreAssembly_sup,ElmatToStokes)
                  ProcPointer%StokesOutside => StokesOutMatrix
               end if
            end if  
         end if
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select   
      
   end subroutine   

 subroutine AllocStokes
      implicit none
      integer(ip) :: ntens

      ntens=(e%ndime-1)*(e%ndime-1)+2      
      call a%Memor%alloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuvFS,'elmuvFS','supm_elmope')      
      call a%Memor%alloc(1,e%mnode,e%ndime,e%mnode,elmuqFS,'elmuqFS','supm_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,1,e%mnode,elmpvFS,'elmpvFS','supm_elmope')
      call a%Memor%alloc(ntens,e%mnode,ntens,e%mnode,elmstFS,'elmstFS','supm_elmope')
      call a%Memor%alloc(ntens,e%mnode,e%ndime,e%mnode,elmutFS,'elmutFS','supm_elmope')      
      call a%Memor%alloc(e%mnode,auxtestf,'auxtestf','supm_elmope')      
   end subroutine
   
   subroutine DeallocStokes
      implicit none
      integer(ip) :: ntens

      ntens=(e%ndime-1)*(e%ndime-1)+2        
      call a%Memor%dealloc(e%ndime,e%mnode,e%ndime,e%mnode,elmuvFS,'elmuvFS','supm_elmope')  
      call a%Memor%dealloc(1,e%mnode,e%ndime,e%mnode,elmuqFS,'elmuqFS','supm_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,1,e%mnode,elmpvFS,'elmpvFS','supm_elmope')
      call a%Memor%dealloc(ntens,e%mnode,ntens,e%mnode,elmstFS,'elmstFS','supm_elmope')
      call a%Memor%dealloc(ntens,e%mnode,e%ndime,e%mnode,elmutFS,'elmutFS','supm_elmope')          
      call a%Memor%dealloc(e%mnode,auxtestf,'auxtestf','supm_elmope')       
   end subroutine 
   
   subroutine InitElmuq
      implicit none
      elmuvFS=0.0_rp
      elmuqFS=0.0_rp
      elmpvFS=0.0_rp
      elmstFS=0.0_rp
      elmutFS=0.0_rp   
   end subroutine
  
   subroutine GalerkinElmuq
      implicit none
      !Compute Contribution to the elemental matrix: Block U,Q -(p,divv)
      call nsm_elmbuq(e,0.0_rp,dvol,acden,LHSdtinv,AGradV,elmuqFS)
      !Compute Contribution to the elemental matrix: Block P,Q tau1(gradp,gradq)
      auxtestf=0.0_rp
      call nsm_elmbpv(e,dvol,auxtestf,elmpvFS)
      call  ProcPointer%StokesOutside
   end subroutine
   
   subroutine StokesOutMatrix
      implicit none  
      !Compute contribution to elemental matrix : Block S,T 1/2acvis(sigma,tau)
      call supm_elmbst(e,0.0_rp,0.0_rp,dvol,acvis,auxtens,a%kfl_splitOSSMomentum,elmstFS)
      !Compute contribution to elemental matrix : Block U,T (Grad_Sym,T)
      call supm_elmbut(e,0.0_rp,0.0_rp,0.0_rp,AgradV,dvol,acden,acvis,auxtens,a%kfl_splitOSSMomentum,elmutFS)  
      !Compute contribution to elemental matrix : Block U,V 2*acvis(Grad_Sym,Grad_Sym)  
      call supm_elmbuv2(e,2.0_rp*acvis,dvol,elmuvFS)       
   end subroutine
      
   subroutine ElmatToStokes
      implicit none
      integer(ip) :: ntens,elemStatus
        
      call a%CutMesh%GetElementType(ielem,elemStatus)
      ntens=(e%ndime-1)*(e%ndime-1)+2
      if(elemStatus==-1)then
         elmat(:,:,:,:) = 0.0_rp        
         elrhs(:,:)=0.0_rp
         !Momentum Equation
         ! Assembly elmuv to elmat
         elmat(ntens+1:ntens+e%ndime,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) = elmat(ntens+1:ntens+e%ndime,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) &
            + elmuvFS(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode) 
         ! Assembly elmpv to elmat
         elmat(ntens+1:ntens+e%ndime,1:e%pnode,ntens+e%ndime+1,1:e%pnode) = elmat(ntens+1:ntens+e%ndime,1:e%pnode,ntens+e%ndime+1,1:e%pnode) &
            + elmpvFS(1:e%ndime,1:e%pnode,1,1:e%pnode)   
         
         !Continuity Equation
         ! Assembly elmuq to elmat
         elmat(ntens+e%ndime+1,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) = elmat(ntens+1+e%ndime,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) &
            + elmuqFS(1,1:e%pnode,1:e%ndime,1:e%pnode)            
         ! Assembly elmpq to elmat
         elmat(ntens+e%ndime+1,1:e%pnode,ntens+e%ndime+1,1:e%pnode) = elmat(ntens+e%ndime+1,1:e%pnode,ntens+e%ndime+1,1:e%pnode) &
            + elmpq(1,1:e%pnode,1,1:e%pnode)
         
         !Constitutive equation
         ! Assembly elmst to elmat
         elmat(1:ntens,1:e%pnode,1:ntens,1:e%pnode) = elmat(1:ntens,1:e%pnode,1:ntens,1:e%pnode) + elmstFS(1:ntens,1:e%pnode,1:ntens,1:e%pnode)
         ! Assembly elmut to elmat
         elmat(1:ntens,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) = elmat(1:ntens,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) + &
               elmutFS(1:ntens,1:e%pnode,1:e%ndime,1:e%pnode)       
                  
      end if
   end subroutine
   
   subroutine ElmatsToLapla
      implicit none
      integer(ip)  :: elemStatus,inode,poinStatus,ipoin,idime,ntens 
      
      call a%CutMesh%GetElementType(ielem,elemStatus)
      ntens=(e%ndime-1)*(e%ndime-1)+2
      if(elemStatus==-1)then
         elmat(:,:,:,:) = 0.0_rp        
         elrhs(:,:)=0.0_rp
         wrmat1=0.0_rp                  
         elmuv=0.0_rp
         
         ! Viscosity terms : we only consider mu*(grad v, grad u)         
         call elmvis(e,dvolt0,acvis,wrmat1)         
         forall (idime = 1:e%ndime)
            elmuv(idime,1:e%pnode,idime,1:e%pnode) = elmuv(idime,1:e%pnode,idime,1:e%pnode) + wrmat1(1:e%pnode,1:e%pnode)
         end forall   
         ! Assembly elmuv to elmat
         elmat(ntens+1:ntens+e%ndime,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) = elmat(ntens+1:ntens+e%ndime,1:e%pnode,ntens+1:ntens+e%ndime,1:e%pnode) &
            + elmuv(1:e%ndime,1:e%pnode,1:e%ndime,1:e%pnode)
         
            do inode=1,e%pnode
               do idime=1,ntens
                  elmat(idime,inode,idime,inode) = elmat(idime,inode,idime,inode) + 1.0_rp
               end do
            end do
         do inode=1,e%pnode         
            elmat(ntens+e%ndime+1,inode,ntens+e%ndime+1,inode) =  1.0_rp
         end do
      end if      
   end subroutine 
   
end module