  module Mod_supm_ComputeGpSplitTerms
   use typre
   use Mod_supm_BaseElmope
   use Mod_supm_InterpolateGradients  
   implicit none
   private
   public SetPointersComputeGpSplitTerms, gpresdivs, gpresconv, gpresgrap
   public gpresdivu, gpresgrau, gpresexp, gpresconvsigma, gpresdeform, gpreslapl
   
   real(rp), allocatable :: gpresdivs(:), gpresconv(:), gpresgrap(:)
   real(rp), allocatable :: gpresdivu(:), gpresgrau(:), gpreslapl(:)
   real(rp), allocatable :: gpresexp(:), gpresconvsigma(:), gpresdeform(:)
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
   
   subroutine SetPointersComputeGpSplitTerms(itask) 
      implicit none
      integer(ip) :: itask
      integer(ip) :: kfl_nonlinear
      
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
            
            if (a%kfl_repro >=2) then
               !We need the gradients in the gauss point
               call SetPointersInterpolateGradients(1)
               
               call ConcatenateProcedures(ProcHook_Initializations,InitGpResTerm)
               call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResTerms)
               call ConcatenateProcedures(ProcHook_Finalizations,FinGpResTerm)
               
               if (a%kfl_repro==4) then
                  call ConcatenateProcedures(ProcHook_Initializations,InitGpResTermSplit)
                  call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResTermsSplit)
                  call ConcatenateProcedures(ProcHook_Finalizations,FinGpResTermSplit)
                 
                  if(a%MatProp(imat)%lawvi<0) then
                     call ConcatenateProcedures(ProcHook_Initializations,InitGpResTermSplit_VE)
                     if (a%LogFormulation==0) then
                        call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResTermsSplitSTD_VE)
                     else if (a%LogFormulation==1) then
                        call ConcatenateProcedures(ProcHook_InGaussElmats,ComputeGpResTermsSplitLCR_VE)
                     end if

                     call ConcatenateProcedures(ProcHook_Finalizations,FinGpResTermSplit_VE)
                  end if
               end if   
               
               call a%Mesh%IsNonLinearElement(ielem,kfl_nonlinear)
               if (kfl_nonlinear == 1) then
                   call ConcatenateProcedures(ProcHook_Initializations,InitOssSplitNL)
                   call ConcatenateProcedures(ProcHook_InGaussElmats,laplaToVector)
                   call ConcatenateProcedures(ProcHook_Finalizations,FinOssSplitNL)
               end if
               
            end if

         end if

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines
   !-------------------------------------------------------------------
 
    subroutine InitGpResTerm
      implicit none
      call a%Memor%alloc(e%ndime,gpresdivs, 'gpresdivs','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,gpresconv, 'gpresconv','supm_EniteElmope')
      call a%Memor%alloc(e%ndime,gpresgrap, 'gpresgrap','supm_EnditeElmope')
      call a%Memor%alloc(e%ndime,elext3,'elext3','supm_EnditeElmope')
   end subroutine
   
   
    subroutine InitGpResTermSplit
      implicit none
      call a%Memor%alloc(1,gpresdivu, 'gpresdivu','supm_EnditeElmope')
      call a%Memor%alloc(auxtens,gpresgrau, 'gpresgrau','supm_EnditeElmope')
      
    end subroutine 
    
    
    subroutine InitGpResTermSplit_VE
      implicit none
      call a%Memor%alloc(auxtens,gpresconvsigma, 'gpresconvsigma','supm_elmope')
      call a%Memor%alloc(auxtens,gpresdeform, 'gpresdeform','supm_elmope')
    end subroutine  
      
   
   
   subroutine FinGpResTerm
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(e%ndime,gpresdivs, 'gpresdivs','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,gpresconv, 'gpresconv','supm_EniteElmope')
      call a%Memor%dealloc(e%ndime,gpresgrap, 'gpresgrap','supm_EnditeElmope')
      call a%Memor%dealloc(e%ndime,elext3,'elext3','supm_EnditeElmope')
   end subroutine
   
   
   subroutine FinGpResTermSplit
      implicit none

      call a%Memor%dealloc(1,gpresdivu, 'gpresdivu','supm_elmope')
      call a%Memor%dealloc(auxtens,gpresgrau, 'gpresgrau','supm_elmope')
      
   end subroutine  
   
   
   subroutine FinGpResTermSplit_VE
      implicit none
         call a%Memor%dealloc(auxtens,gpresconvsigma, 'gpresconvsigma','supm_elmope')
         call a%Memor%dealloc(auxtens,gpresdeform, 'gpresdeform','supm_elmope')

   end subroutine  
      

   subroutine ComputeGpResTerms
      implicit none
      integer(ip) :: i
      elext3=0.0_rp
      !Compute vector of external forces. This is important in free-surface and two fluids problem
      call nsi_ComputeExternalForces(e,acden,a%grnor,a%gravi,elext3)  
      
      !SplitOSS momentum eq. terms
      call sup_termTovector(e%ndime,auxtens,gpvel,grvel,grsig,GrExpPsi,ExpGpPsi_Matrix,grpre,elext3,gpresconv,gpresdivs,gpresgrap,a%LogFormulation) 
      
   end subroutine
   
   
   subroutine ComputeGpResTermsSplit
      implicit none 
     !Three Field SplitOSS Dynamic
      call sup_termTovector2(e%ndime,auxtens,gpvel,grvel,gpresdivu,gpresgrau)
   end subroutine   
   
   subroutine ComputeGpResTermsSplitSTD_VE
      implicit none
      !Viscoelastic SplitOSS Dynamic - constitutive eq. terms
       call sup_termTovector_viscoelastic(e%ndime,auxtens,gpvel,grvel,gpsig,grsig,gpresconvsigma,gpresdeform) 
   end subroutine   
   
   
   subroutine ComputeGpResTermsSplitLCR_VE
      implicit none
      !Viscoelastic SplitOSS Dynamic - constitutive eq. terms
      call sup_termTovector_log(e%ndime,auxtens,gpvel,grvel,GrExpPsi,ExpGpPsi_Matrix,gpresconvsigma,gpresdeform)
   end subroutine   
   

   subroutine InitOssSplitNL
      implicit none
      !Matrices alloc
      call a%Memor%alloc(e%ndime,gpreslapl, 'gpreslapl','supm_EnditeElmope')
   end subroutine
   
   subroutine FinOssSplitNL
      implicit none        
      !Matrices dealloc
      call a%Memor%dealloc(e%ndime,gpreslapl, 'gpreslapl','supm_EnditeElmope')
   end subroutine
      

   subroutine laplaToVector
      implicit none
      call supm_laplaTovector(e,acvis,beta,elvel,gpreslapl)
   end subroutine
   

   subroutine supm_elmrfe_split(e,acden,veadv, &
         grapr,grave,resim)
      implicit none
      class(FiniteElement) :: e
      real(rp) :: acden,veadv(e%ndime),grapr(e%ndime),grave(e%ndime,e%ndime),resim(2*e%ndime+1)
      integer(ip) :: idime,jdime
      
      resim = 0.0_rp
      
      !Contribution from the convective term
      do jdime = 1,e%ndime
         resim(1:e%ndime) = resim(1:e%ndime) + acden*veadv(jdime)*grave(:,jdime) 
      end do

      do idime = 1,e%ndime                              ! Contribution from the divergence term
         resim(e%ndime+1) = resim(e%ndime+1) + grave(idime,idime)
      end do
      
      !Pressure Gradient, stored separately at the end in this case
      resim(e%ndime+2:2*e%ndime+1) = resim(e%ndime+2:2*e%ndime+1) + grapr     ! Contribution from pressure term

   end subroutine supm_elmrfe_split
   
end module


