module Mod_nsc_ComputeConvectionCoefficients
   use typre
   use Mod_nsc_BaseElmope
   use Mod_nsc_ArbitraryLagrangianEulerian
   implicit none
   private
   public SetPointersConvectionCoefficients
   
   integer(ip), allocatable :: kfl_IsSet

   !ALE
   logical :: isALE
   
contains

   subroutine SetPointersConvectionCoefficients(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateConvectionCoefficients)
            call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,ConvectionCoefficientsToZero)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateConvectionCoefficients)
            if (a%kfl_advec == 1) then
               call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,ComputeConvectionCoefficients)
               if (a%lawde == 1) then
                  call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,ComputeIdealConvectionCoefficients)
               end if
            end if

            call a%Mesh%GetALE(isALE)
            if (isALE) then
               call SetPointersArbitraryLagrangianEulerian(1)
               call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,ComputeALECoefficients)
            end if

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateConvectionCoefficients

      call a%Memor%alloc(e%mnode,Add,'Add','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Adm,'Adm','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Amd,'Amd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Amm,'Amm','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Ame,'Ame','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,Aed,'Aed','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Aem,'Aem','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,Aee,'Aee','nsc_elmope_im')

   end subroutine

   !Convection matrix coefficients
   subroutine ConvectionCoefficientsToZero
      implicit none    

      Add = 0.0_rp !Add(p)
      Adm = 0.0_rp !Adm(d,p)
      Amd = 0.0_rp !Amd(i,p)
      Amm = 0.0_rp !Amm(i,d,p)
      Ame = 0.0_rp !Ame(i,p)
      Aed = 0.0_rp !Aed(p)
      Aem = 0.0_rp !Aem(d,p)
      Aee = 0.0_rp !Aee(p)
                        
   end subroutine   
!Include derivatives here  A_j*d_j 
   subroutine ComputeConvectionCoefficients
      implicit none    
      integer(ip)                :: idime

      !Mass equation
      Adm = Adm + e%cartd !d_dj*d_j 
      !Momentum equation
      do idime=1,e%ndime
         Amd(idime,:) = Amd(idime,:) - gpadv(idime)*AGradV(:)!-vel_i*vel_j*d_j 
         Amm(idime,:,:) = Amm(idime,:,:) + gpadv(idime)*e%cartd(:,:)!(vel_i)*d_dj*d_j 
         Amm(idime,idime,:) = Amm(idime,idime,:) + AGradV(:)!(vel_j)*d_id*d_j 
      end do
      !Energy equation
      Aed = Aed - (gpade+gppre)*invgpd*AGradV !-(e+p)/rho vel_j*d_j   
      Aem = Aem + (gpade+gppre)*invgpd*e%cartd!(e+p)/rho d_dj*d_j

   end subroutine   

   !-------------------------------------------------------------------
   !Ideal state law convection matrix coefficients
!Include derivatives here  A_j*d_j 
   subroutine ComputeIdealConvectionCoefficients
      implicit none
      integer(ip)                :: idime
      
      !Momentum equation
      Amd = Amd + aux_d*e%cartd!(gamma-1)(vel·vel/2) d_ij*d_j    
      Ame = Ame + aux*e%cartd!(gamma-1) d_ij*d_j   
      do idime=1,e%ndime
         Amm(:,idime,:) = Amm(:,idime,:) - aux*gpadv(idime)*e%cartd(:,:)!-(gamma-1)(vel_d)*d_ij*d_j 

      !Energy equation
         Aem(idime,:) = Aem(idime,:) - aux*gpadv(idime)*AGradV(:)!-(gamma-1)*vel_d*vel_j*d_j 
      end do
      Aed = Aed + aux_d*AGradV !(gamma-1)(vel·vel/2) vel_j*d_j   
      Aee = Aee + acgamma*AGradV ! (gamma)·vel_j*d_j  

   end subroutine   

   !-------------------------------------------------------------------
   !ALE convection matrix coefficients
   subroutine ComputeALECoefficients
      implicit none    
      integer(ip)                :: idime

      !Mass equation
      Add = Add - ALEGradV !-vmesh_j*d_j 
      !Momentum equation
      do idime=1,e%ndime
         Amm(idime,idime,:) = Amm(idime,idime,:) - ALEGradV(:) !(-vmesh_j)*d_id*d_j 
      end do
      !Energy equation
      Aee = Aee - ALEGradV !-vmesh_j*d_j    

   end subroutine   

   subroutine DeallocateConvectionCoefficients

      call a%Memor%dealloc(e%mnode,Add,'Add','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Adm,'Adm','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Amd,'Amd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Amm,'Amm','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Ame,'Ame','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,Aed,'Aed','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Aem,'Aem','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,Aee,'Aee','nsc_elmope_im')

   end subroutine

end module
