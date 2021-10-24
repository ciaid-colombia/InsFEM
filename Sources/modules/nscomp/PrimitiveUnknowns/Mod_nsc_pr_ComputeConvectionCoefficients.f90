module Mod_nsc_pr_ComputeConvectionCoefficients
   use typre
   use Mod_nsc_pr_BaseElmope
   use Mod_nsc_pr_ArbitraryLagrangianEulerian
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
         
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateConvectionCoefficients)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeConvectionCoefficients,ConvectionCoefficientsToZero)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateConvectionCoefficients)
            if (a%kfl_advec == 1) then
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeConvectionCoefficients,ComputeConvectionCoefficients)
            end if

            call a%Mesh%GetALE(isALE)
            if (isALE) then
               call SetPointersArbitraryLagrangianEulerian(1)
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeConvectionCoefficients,ComputeALECoefficients)
            end if

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateConvectionCoefficients

      call a%Memor%alloc(e%mnode,Add,'Add','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Adm,'Adm','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Ade,'Ade','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Amd,'Amd','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Amm,'Amm','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Ame,'Ame','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Aed,'Aed','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Aem,'Aem','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Aee,'Aee','nsc_pr_elmope')

   end subroutine

   !Convection matrix coefficients
   subroutine ConvectionCoefficientsToZero
      implicit none    

      Add = 0.0_rp !Add(p)
      Adm = 0.0_rp !Adm(d,p)
      Ade = 0.0_rp !Adm(p)
      Amd = 0.0_rp !Amd(i,p)
      Amm = 0.0_rp !Amm(i,d,p)
      Ame = 0.0_rp !Ame(i,p)
      Aed = 0.0_rp !Aed(p)
      Aem = 0.0_rp !Aem(d,p)
      Aee = 0.0_rp !Aee(p)
                        
   end subroutine   

   subroutine ComputeConvectionCoefficients
      implicit none    
      integer(ip)                :: idime

      !Mass equation
      Add = Add + gpden*acbeta*AGradV !rho*beta*vel_j*d_j 
      Adm = Adm + gpden*e%cartd !rho*d_dj*d_j 
      Ade = Ade - gpden*acalpha*AGradV !-rho*alpha*vel_j*d_j 
      !Momentum equation
      Amd = Amd + e%cartd !d_ij*d_j 
      do idime=1,e%ndime
         Amd(idime,:) = Amd(idime,:) + gpden*acbeta*gpadv(idime)*AGradV(:)!rho*beta*vel_i*vel_j*d_j 
         Amm(idime,:,:) = Amm(idime,:,:) + gpden*gpadv(idime)*e%cartd(:,:)!rho*(vel_i)*d_dj*d_j 
         Amm(idime,idime,:) = Amm(idime,idime,:) + gpden*AGradV(:)!rho*(vel_j)*d_id*d_j 
         Ame(idime,:) = Ame(idime,:) - gpden*acalpha*gpadv(idime)*AGradV(:)!-rho*alpha*vel_i*vel_j*d_j
      !Energy equation
         Aem(idime,:) = Aem(idime,:) + gpden*gpadv(idime)*AGradV(:)!rho*vel_d*vel_j*d_j 
      end do
      !Energy equation
      Aed = Aed + aux*acbeta*AGradV !rho*(enthalpy+kineticEnerg)*beta*vel_j*d_j
      Aed = Aed - acalpha*gpadt*AGradV !-alpha*temperature*vel_j*d_j
      Aed = Aed + AGradV !vel_j*d_j   
      Aem = Aem + aux*e%cartd !rho*(enthalpy+kineticEnerg)*d_dj*d_j
      Aee = Aee - aux*acalpha*AGradV !- rho*(enthalpy+kineticEnerg)*alpha*vel_j*d_j  
      Aee = Aee + gpden*accph*AGradV !rho*cp*vel_j*d_j   

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

      call a%Memor%dealloc(e%mnode,Add,'Add','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Adm,'Adm','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Ade,'Ade','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Amd,'Amd','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Amm,'Amm','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Ame,'Ame','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Aed,'Aed','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Aem,'Aem','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Aee,'Aee','nsc_pr_elmope')

   end subroutine

end module
