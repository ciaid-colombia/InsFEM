 module Mod_supm_ComputeTauSmoothing
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersComputeTauSmoothing

   
   !Tau Smoothing
   real(rp), allocatable :: elTauSmoothing(:,:)
   real(rp), allocatable :: elSmoothedTau2(:,:)
   real(rp), allocatable :: wTausmo(:,:)

   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SetPointersComputeTauSmoothing(itask)
      implicit none
      integer(ip) :: itask
      character(6) :: task
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            if(a%kfl_Tausm == 1) then
               call ConcatenateProcedures(ProcHook_Initializations,AllocTauSmoothing)
               call ConcatenateProcedures(ProcHook_Finalizations,DeallocTauSmoothing)
               call ConcatenateProcedures(ProcHook_PreGauss,ElmatsToZeroTauSmoothing)


               if (a%MatProp(imat)%lawvi >= 0) then
                  call ConcatenateProcedures(ProcHook_InGaussElmats,GpTauToElTauSmoothing)
               elseif(a%MatProp(imat)%lawvi < 0) then
                  if(a%LogFormulation==0) then
                     call ConcatenateProcedures(ProcHook_InGaussElmats,GpTauToElTauSmoothingVE_STD)
                  else
                     call ConcatenateProcedures(ProcHook_InGaussElmats,GpTauToElTauSmoothingVE_LCR)
                  end if
               end if
               
               call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyTausmo)
               call ConcatenateProcedures(ProcHook_PostLoop,SmoothTausmo)
            end if
         endif 
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   subroutine AllocTauSmoothing
      implicit none
      
      integer(ip) :: npoin
      
      call a%Memor%alloc(3,e%mnode,elTauSmoothing,'elTauSmoothing','supm_EndElmope')
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%alloc(3,npoin,wTausmo,'wTausmo','supm_EndElmope')
   end subroutine
   
   subroutine DeallocTauSmoothing
      implicit none
      
      integer(ip) :: npoin
      
      call a%Memor%dealloc(3,e%mnode,elTauSmoothing,'elTauSmoothing','supm_EndElmope')
   end subroutine
   
   subroutine ElmatsToZeroTauSmoothing
      implicit none
      
      elTauSmoothing = 0.0_rp
   end subroutine
   
   subroutine GpTauToElTauSmoothing
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none
      real(rp) :: timomaux,tidivaux,tisigaux,w2,w3
      integer(ip) :: inode
      real(rp) :: auxG
      
      
      !We compute the Taus at the Gauss point
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timomaux)
      call sup_ComputeTidiv(e,acden,acvis,gpvno,a%staco,chale,tidivaux)
      call sup_ComputeTidiv(e,acden,acvis,gpvno_VE,a%staco(4),chale,tisigaux)

      !And we assembly them
      do inode = 1,e%pnode
         elTauSmoothing(1,inode) = elTauSmoothing(1,inode) + e%shape(inode,e%igaus)*timomaux*dvol
         elTauSmoothing(2,inode) = elTauSmoothing(2,inode) + e%shape(inode,e%igaus)*tidivaux*dvol
         elTauSmoothing(3,inode) = elTauSmoothing(3,inode) + e%shape(inode,e%igaus)*tisigaux*dvol
      enddo
      
   end subroutine
   
   
   subroutine GpTauToElTauSmoothingVE_STD
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none
      
      real(rp) :: timomaux,tidivaux,tisigaux,w2,w3
      integer(ip) :: inode
      real(rp) :: auxG
      
      !We compute the Taus at the Gauss point
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timomaux)
      call sup_ComputeTidivVE(e,timomaux,a%staco,chale,tidivaux)
      auxG   = a%MatProp(imat)%LawviParam(4)/((1.0_rp-beta) + 0.00001_rp)
      call sup_ComputeTisig(e,lambda,acvis,gpvno_VE,grvel,a%staco,chale,auxG,a%PTT_model,gpsig,tisigaux)   

      !And we assembly them
      do inode = 1,e%pnode
         elTauSmoothing(1,inode) = elTauSmoothing(1,inode) + e%shape(inode,e%igaus)*timomaux*dvol
         elTauSmoothing(2,inode) = elTauSmoothing(2,inode) + e%shape(inode,e%igaus)*tidivaux*dvol
         elTauSmoothing(3,inode) = elTauSmoothing(3,inode) + e%shape(inode,e%igaus)*tisigaux*dvol
      enddo
      
   end subroutine
   
   
   subroutine GpTauToElTauSmoothingVE_LCR
      use Mod_sup_ComputeTidivVE
      use Mod_sup_ComputeTisig
      implicit none
      real(rp) :: timomaux,tidivaux,tisigaux,w2,w3
      integer(ip) :: inode
      real(rp) :: auxG

      !We compute the Taus at the Gauss point
      call ComputeTau(e,acden,acvis,gpvno,a%staco,chale,timomaux)

      call sup_ComputeTidivVE(e,timomaux,a%staco,chale,tidivaux)
      auxG   = (a%MatProp(imat)%LawviParam(4)*lambda)/(2.0_rp*(lambda0**2_ip) + 0.00001_rp) !GIESEKUSLCR
      call sup_LCR_ComputeTisig(e,lambda,lambda0,auxG,beta,acvis,gpvno_VE,grvel,ExpGpPsi_Matrix,a%staco,chale,a%PTT_model,gpsig,w2,w3,tisigaux) 

      !And we assembly them
      do inode = 1,e%pnode
         elTauSmoothing(1,inode) = elTauSmoothing(1,inode) + e%shape(inode,e%igaus)*timomaux*dvol
         elTauSmoothing(2,inode) = elTauSmoothing(2,inode) + e%shape(inode,e%igaus)*tidivaux*dvol
         elTauSmoothing(3,inode) = elTauSmoothing(3,inode) + e%shape(inode,e%igaus)*tisigaux*dvol
      enddo
      
   end subroutine
   
   subroutine AssemblyTausmo
      implicit none
      call a%Mesh%AssemblyToArray(e,size(wTausmo,1),elTauSmoothing,wTausmo)
   end subroutine
   
   subroutine SmoothTausmo
      implicit none
      integer(ip) :: npoin, nelem, ismooth, nsmooth, inode
      real(rp)    :: timom, tidiv, tisig
      real(rp)    :: timomaux,tidivaux,tisigaux
      real(rp)    :: gptau(3)
      
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supm_EnditeElmope')  
      call a%Mesh%ElementAlloc(e,a%Memor,'DefaultRule','supm_EnditeElmope') 
      
      call a%Mesh%GetNelem(nelem)

      call a%Mesh%Smooth(3,wTausmo) 
      call move_alloc(wTausmo,a%Tausmo)
      !deallocate(wTausmo)
      call a%Mesh%GetNpoin(npoin)
      call a%Memor%deallocObj(0,'wTausmo','supm_EndElmope',3*npoin*rp)

      !Re-smoothing n times
      if (a%kfl_ntausmooth>1) then
         call a%Memor%alloc(3,e%mnode,elSmoothedTau2,'elSmoothedTau2','sup_EndElmope')
         call a%Memor%alloc(3,e%mnode,elTauSmoothing,'elTauSmoothing','supm_EndElmope')
         do ismooth=2,a%kfl_ntausmooth            
            call a%Mesh%GetNpoin(npoin)
            call a%Memor%alloc(3,npoin,wTausmo,'wTausmo','supm_EndElmope')
            do ielem=1,nelem
               call a%Mesh%ElementLoad(ielem,e)
               call e%elmdcg
               call e%gather(3,elSmoothedTau2,a%Tausmo)
               
               elTauSmoothing = 0.0_rp
               

               do igaus= 1,e%pgaus
                  e%igaus=igaus
                  dvol = e%weigp(e%igaus)*e%detjm
                 
                  call e%interpg(3,elSmoothedTau2,gptau)
                  timomaux = gptau(1)
                  tidivaux = gptau(2)
                  tisigaux = gptau(3)
                  
                  !And we assembly them
                  do inode = 1,e%pnode
                     elTauSmoothing(1,inode) = elTauSmoothing(1,inode) + e%shape(inode,e%igaus)*timomaux*dvol
                     elTauSmoothing(2,inode) = elTauSmoothing(2,inode) + e%shape(inode,e%igaus)*tidivaux*dvol
                     elTauSmoothing(3,inode) = elTauSmoothing(3,inode) + e%shape(inode,e%igaus)*tisigaux*dvol
                  enddo
               
               end do
               
               call a%Mesh%AssemblyToArray(e,size(wTausmo,1),elTauSmoothing,wTausmo)
               
            end do

            a%Tausmo=0.0_rp
            call a%Mesh%Smooth(3,wTausmo) 
            call move_alloc(wTausmo,a%Tausmo)
            !deallocate(wTausmo)
            call a%Mesh%GetNpoin(npoin)
            call a%Memor%deallocObj(0,'wTausmo','supm_EndElmope',3*npoin*rp)
         end do
         
         call a%Memor%dealloc(3,e%mnode,elSmoothedTau2,'elSmoothedTau2','sup_EndElmope')
         call a%Memor%dealloc(3,e%mnode,elTauSmoothing,'elTauSmoothing','supm_EndElmope')
         
      end if
      
      call a%Mesh%ElementDeAlloc(e,a%Memor,'DefaultRule','supm_EnditeElmope')  
      call a%Mesh%ElementAlloc(e,a%Memor,a%EndLoopQuadrature,'supm_EnditeElmope') 
      
   end subroutine
   
end module

