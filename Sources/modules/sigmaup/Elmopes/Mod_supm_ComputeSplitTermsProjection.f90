module Mod_supm_ComputeSplitTermsProjection
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_ComputeGpSplitTerms
   use Mod_supm_InterpolateGradients
   implicit none
   private
   public SetPointersComputeSplitTermsProjection

   real(rp), allocatable :: elresdivu(:,:) , elresgrau(:,:) , elresexp(:,:)  
   real(rp), allocatable :: elresconvsigma(:,:) , elresdeform(:,:) 
   integer(ip), allocatable :: kfl_IsSet
   integer(ip) :: kfl_nonlinear
 
contains
   
   subroutine SetPointersComputeSplitTermsProjection(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1

            if (a%kfl_repro >= 2) then !Split OSS, ASGS and OSS term by term

               call ConcatenateProcedures(ProcHook_Initializations,InitResSplit)
               call ConcatenateProcedures(ProcHook_PreLoop,PreLoopOssSplit)
               call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroOssSplit)
            
               !We need to compute the residual at the gauss point
               call SetPointersComputeGpSplitTerms(1)
               !Now we assembly everything for computing the projection
               call ConcatenateProcedures(ProcHook_InGaussElmats,GaussPointAssemblyOssSplit)         
               call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyOssSplit)     
               
               if(a%MatProp(imat)%lawvi<0)then
                  call ConcatenateProcedures(ProcHook_PostLoop,ResidualBoundaryOssSplit)
               endif
               
               call ConcatenateProcedures(ProcHook_PostLoop,SplitTermProjection)
               call ConcatenateProcedures(ProcHook_Finalizations,FinResSplit)

               if (a%kfl_repro==4) then
                  call ConcatenateProcedures(ProcHook_Initializations,InitResSplitConstitutive)
                  call ConcatenateProcedures(ProcHook_PreLoop,PreLoopOssSplitConstitutive)
                  call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroOssSplitConstitutive)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,GaussPointAssemblyOssSplitConstitutive) 
                  call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyOssSplitConstitutive)   
                  call ConcatenateProcedures(ProcHook_PostLoop,SplitTermProjectionConstitutive)
                  call ConcatenateProcedures(ProcHook_Finalizations,FinResSplitConstitutive)
                  
                   if (a%MatProp(imat)%lawvi<0) then
                     call ConcatenateProcedures(ProcHook_Initializations,InitResSplitConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_PreLoop,PreLoopOssSplitConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroOssSplitConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_InGaussElmats,GaussPointAssemblyOssSplitConstitutive_VE) 
                     call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyOssSplitConstitutive_VE)   
                     call ConcatenateProcedures(ProcHook_PostLoop,SplitTermProjectionConstitutive_VE)
                     call ConcatenateProcedures(ProcHook_Finalizations,FinResSplitConstitutive_VE)
                   end if
               end if
               
               if (kfl_nonlinear==1)then !Add the non-linear element
                  call ConcatenateProcedures(ProcHook_Initializations,InitOssSplitNL)
                  call ConcatenateProcedures(ProcHook_PreLoop,PreLoopOssSplitNL)
                  call ConcatenateProcedures(ProcHook_ElmatsToZero,ElmatsToZeroOssSplitNL)
                  !We need to compute the residual at the gauss point
                  call SetPointersComputeGpSplitTerms(1)
                  !Now we assembly everything for computing the projection
                  call ConcatenateProcedures(ProcHook_InGaussElmats,GaussPointAssemblyOssSplitNL)         
                  call ConcatenateProcedures(ProcHook_AssemblyEndite,AssemblyOssSplitNL)     

                  if(a%MatProp(imat)%lawvi<0)then
                     call ConcatenateProcedures(ProcHook_PostLoop,ResidualBoundaryOssSplitNL)
                  endif
                  call ConcatenateProcedures(ProcHook_PostLoop,SmoothOssSplitNL)
                  call ConcatenateProcedures(ProcHook_Finalizations,FinOssSplitNL)

               end if
               
            endif
            
         endif  

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      end select
      
   end subroutine  
   
!------------------------------------------------
!SUBROUTINES


   subroutine InitResSplit
      implicit none
      
      call a%Memor%alloc(e%ndime,e%mnode,elresdivs, 'elresdivs','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elresconv, 'elresconv','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,e%mnode,elresgrap, 'elresgrap','supm_EnditeElmope')
      
   end subroutine
   
   
   subroutine InitResSplitConstitutive
      implicit none

      call a%Memor%alloc(1,e%mnode,elresdivu, 'elresdivu','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,e%mnode,elresgrau, 'elresgrau','supm_EnditeElmope')
      

   end subroutine   
   
   subroutine InitResSplitConstitutive_VE
      implicit none

      call a%Memor%alloc(auxtens,e%mnode,elresconvsigma, 'elresconvsigma','supm_elmope')
      call a%Memor%alloc(auxtens,e%mnode,elresdeform, 'elresdeform','supm_elmope')    

   end subroutine   
   
   subroutine PreLoopOssSplit
      implicit none
      a%reproSDiv   = 0.0_rp !Div(S)
      a%reproUGradU = 0.0_rp !u*grad(u) (Newton terms)
      a%reproGradP  = 0.0_rp !grad(p)
   end subroutine 
   
   subroutine PreLoopOssSplitConstitutive
      implicit none
      a%reproDivU = 0.0_rp !Div(U)
      a%reproGradU  = 0.0_rp !grad_sym(u)
   end subroutine   
   
   
   subroutine PreLoopOssSplitConstitutive_VE
      implicit none
      a%reproUGradS = 0.0_rp !u*grad(S)
      a%reproSGradU = 0.0_rp !s*grad(u)
   end subroutine  
   
   
   subroutine ElmatsToZeroOssSplit
      implicit none
      elresdivs = 0.0_rp
      elresconv = 0.0_rp
      elresgrap = 0.0_rp   
   end subroutine   
   
   
   subroutine ElmatsToZeroOssSplitConstitutive
      implicit none
      elresdivu      = 0.0_rp !div(u)
      elresgrau      = 0.0_rp !grad_sym(u)
   end subroutine   
   
      subroutine ElmatsToZeroOssSplitConstitutive_VE
      implicit none

      elresconvsigma  = 0.0_rp !u*grad(S)
      elresdeform     = 0.0_rp !s*grad(u)
   end subroutine   
   

   subroutine FinResSplit
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(e%ndime,e%mnode,elresdivs, 'elresdivs','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elresconv, 'elresconv','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,e%mnode,elresgrap, 'elresgrap','supm_EnditeElmope') 
           
   end subroutine
   
   
    subroutine FinResSplitConstitutive
      implicit none        
      
      call a%Memor%dealloc(1,e%mnode,elresdivu, 'elresdivu','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elresgrau, 'elresgrau','supm_elmope')
   end subroutine
   
   
   subroutine FinResSplitConstitutive_VE
      implicit none        
      
      call a%Memor%dealloc(auxtens,e%mnode,elresconvsigma, 'elrepconvsigma','supm_elmope')
      call a%Memor%dealloc(auxtens,e%mnode,elresdeform, 'elrepdeform','supm_elmope')  

   end subroutine
   
   subroutine GaussPointAssemblyOssSplit
      implicit none
      integer(ip) :: elemStatus
      
      !this is crucial when we have two fluids and enriched pressure      
      if(a%kfl_colev==1)then
         call a%CutMesh%GetElementType(ielem,elemStatus)
         if(elemStatus ==0)then 
            gpresdivs=0.0_rp
            gpresconv=0.0_rp
            gpresgrap=0.0_rp

         end if         
      end if      

      call supm_elmrep(e,dvol,e%ndime,gpresdivs,elresdivs)
      call supm_elmrep(e,dvol,e%ndime,gpresconv,elresconv)
      call supm_elmrep(e,dvol,e%ndime,gpresgrap,elresgrap)
   
   end subroutine   
   
   
   subroutine GaussPointAssemblyOssSplitConstitutive
      implicit none
      integer(ip) :: elemStatus
      
      !this is crucial when we have two fluids and enriched pressure      
      if(a%kfl_colev==1)then
         call a%CutMesh%GetElementType(ielem,elemStatus)
         if(elemStatus ==0)then 
            
            gpresdivu=0.0_rp
            gpresgrau=0.0_rp
            
            if (a%LogFormulation==1) then
               gpresconvsigma=0.0_rp
               gpresdeform=0.0_rp
            end if
         end if         
      end if      

      call supm_elmrep(e,dvol,auxtens,gpresgrau,elresgrau)
      call supm_elmrep(e,dvol,1,gpresdivu,elresdivu)
     

   end subroutine   
   
   
   
   subroutine GaussPointAssemblyOssSplitConstitutive_VE
      implicit none
      integer(ip) :: elemStatus
   
      call supm_elmrep(e,dvol,auxtens,gpresconvsigma,elresconvsigma)
      call supm_elmrep(e,dvol,auxtens,gpresdeform,elresdeform)

   end subroutine   
   
   subroutine AssemblyOssSplit
      implicit none    
      a%reproSDiv(:,e%lnods(1:e%pnode)) = a%reproSDiv(:,e%lnods(1:e%pnode)) + elresdivs(:,1:e%pnode)
      a%reproUGradU(:,e%lnods(1:e%pnode)) = a%reproUGradU(:,e%lnods(1:e%pnode)) + elresconv(:,1:e%pnode) 
      a%reproGradP(:,e%lnods(1:e%pnode)) = a%reproGradP(:,e%lnods(1:e%pnode)) + elresgrap(:,1:e%pnode)

   end subroutine
   
   
   subroutine AssemblyOssSplitConstitutive
      implicit none    

      a%reproGradU(:,e%lnods(1:e%pnode)) = a%reproGradU(:,e%lnods(1:e%pnode)) + elresgrau(:,1:e%pnode)
      a%reproDivU(1,e%lnods(1:e%pnode)) = a%reproDivU(1,e%lnods(1:e%pnode)) + elresdivu(1,1:e%pnode)
   end subroutine
   
   
   
   subroutine AssemblyOssSplitConstitutive_VE
      implicit none    

         a%reproUGradS(:,e%lnods(1:e%pnode)) = a%reproUGradS(:,e%lnods(1:e%pnode)) + elresconvsigma(:,1:e%pnode) 
         a%reproSGradU(:,e%lnods(1:e%pnode)) = a%reproSGradU(:,e%lnods(1:e%pnode)) + elresdeform(:,1:e%pnode)

   end subroutine
   
   subroutine ResidualBoundaryOssSplit
      call a%Mesh%Getnpoin(npoin)      
         do ipoin=1,npoin      
            call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
            if(ibopo>0 .and. a%kfl_reproBoundzero==1 .and. a%LogFormulation==0 .and. a%kfl_repro/=4) then
               a%reproSDiv(:,ipoin) = 0.0_rp
               a%reproUGradU(:,ipoin) = 0.0_rp         
               a%reproGradP(:,ipoin) = 0.0_rp           
            end if
         end do     
   end subroutine
      
   
   subroutine SplitTermProjection
      implicit none
      call a%Project(e%ndime,a%reproSDiv)
      call a%Project(e%ndime,a%reproUGradU)    
      call a%Project(e%ndime,a%reproGradP) 

   end subroutine
   

   subroutine SplitTermProjectionConstitutive
      implicit none
      call a%Project(auxtens,a%reproGradU)
      call a%Project(1,a%reproDivU)
      

      if (a%LogFormulation==1) then
         call a%Project(auxtens,a%reproExpS)
      end if
   end subroutine
   
   
   
   
   subroutine SplitTermProjectionConstitutive_VE
      implicit none
       
      if (a%MatProp(imat)%lawvi<0) then
            call a%Project(auxtens,a%reproUGradS)    
            call a%Project(auxtens,a%reproSGradU) 
      end if

   end subroutine
      
      
   !--------------------------------------   
   !Non linear elements
   !-------------------------------------
   
    subroutine PreLoopOssSplitNL
      implicit none
      a%reproLapla   = 0.0_rp ! lapla(U)
   end subroutine   
   
   subroutine InitOssSplitNL
      implicit none
      !Matrices alloc
      call a%Memor%alloc(e%ndime,e%mnode,elreslapl, 'elreslapl','supm_EnditeElmope')
   end subroutine
   
   subroutine ElmatsToZeroOssSplitNL
      implicit none
      elreslapl = 0.0_rp   
   end subroutine   

   subroutine FinOssSplitNL
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(e%ndime,e%mnode,elreplapl, 'elreslapl','supm_EnditeElmope') 
   end subroutine
      
   
   subroutine GaussPointAssemblyOssSplitNL
      implicit none 
      integer(ip) :: elemStatus
      if(a%kfl_colev==1)then
         call a%CutMesh%GetElementType(ielem,elemStatus)
         if(elemStatus ==0)then 
            gpreslapl=0.0_rp
         end if         
      end if       
      call supm_elmrep(e,dvol,e%ndime,gpreslapl,elreslapl)
   end subroutine   
   
   subroutine AssemblyOssSplitNL
      implicit none    
      a%reproLapla(:,e%lnods(1:e%pnode)) = a%reproLapla(:,e%lnods(1:e%pnode)) + elreslapl(:,1:e%pnode)
   end subroutine
   
   subroutine ResidualBoundaryOssSplitNL
      call a%Mesh%Getnpoin(npoin)      
      do ipoin=1,npoin      
         call a%Mesh%GetBoundaryPoint(ipoin,ibopo,exnor)
         if(ibopo>0)then 
            a%reproLapla(:,ipoin) = 0.0_rp
         end if
      end do    
   end subroutine     
   
   subroutine SmoothOssSplitNL
      implicit none
      call a%Project(e%ndime,a%reproLapla)
   end subroutine      
   
   
end module
