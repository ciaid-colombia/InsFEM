module Mod_nsc_pr_ComputeTransientCoefficients
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersTransientCoefficients
   
   integer(ip), allocatable :: kfl_IsSet

contains

   subroutine SetPointersTransientCoefficients(itask)
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
         
            call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateTransientCoefficients)
            call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeTransientCoefficients,TransientCoefficientsToZero)
            call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateTransientCoefficients)
            if (a%kfl_timei == 1) then
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeTransientCoefficients,ComputeTransientCoefficients)
            end if

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateTransientCoefficients

      call a%Memor%alloc(e%mnode,Atdd,'Atdd','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Atde,'Atde','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Atmd,'Atmd','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Atmm,'Atmm','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Atme,'Atme','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Ated,'Ated','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,Atem,'Atem','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,Atee,'Atee','nsc_pr_elmope')

   end subroutine

   !Transient matrix coefficients
   subroutine TransientCoefficientsToZero
      implicit none    
                         
      Atdd = 0.0_rp !Atdd(p)
      Atde = 0.0_rp !Atdm(d,p)
      Atmd = 0.0_rp !Atmd(i,p)
      Atmm = 0.0_rp !Atmm(i,d,p)
      Atme = 0.0_rp !Atme(i,p)
      Ated = 0.0_rp !Ated(p)
      Atem = 0.0_rp !Atem(d,p)
      Atee = 0.0_rp !Atee(p)
                        
   end subroutine   

   subroutine ComputeTransientCoefficients
      implicit none    
      integer(ip)                :: idime

      !Mass equation
      Atdd(1:e%pnode) = Atdd(1:e%pnode) + gpden*acbeta*e%shape(1:e%pnode,e%igaus) !rho*beta
      Atde(1:e%pnode) = Atde(1:e%pnode) - gpden*acalpha*e%shape(1:e%pnode,e%igaus) !-rho*alpha
      !Momentum equation
      do idime=1,e%ndime
         Atmd(idime,1:e%pnode) = Atmd(idime,1:e%pnode) + gpden*acbeta*gpadv(idime)*e%shape(1:e%pnode,e%igaus)!rho*beta*vel_i
         Atmm(idime,idime,1:e%pnode) = Atmm(idime,idime,1:e%pnode) + gpden*e%shape(1:e%pnode,e%igaus)!rho*d_id
         Atme(idime,1:e%pnode) = Atme(idime,1:e%pnode) - gpden*acalpha*gpadv(idime)*e%shape(1:e%pnode,e%igaus)!-rho*alpha*vel_i
      !Energy equation
         Atem(idime,1:e%pnode) = Atem(idime,1:e%pnode) + gpden*gpadv(idime)*e%shape(1:e%pnode,e%igaus)!rho*vel_d 
      end do
      Ated(1:e%pnode) = Ated(1:e%pnode) + aux*acbeta*e%shape(1:e%pnode,e%igaus) !rho*(enthalpy+kineticEnerg)*beta
      Ated(1:e%pnode) = Ated(1:e%pnode) - acalpha*gpadt*e%shape(1:e%pnode,e%igaus) !-alpha*temperature
      Atee(1:e%pnode) = Atee(1:e%pnode) - aux*acalpha*e%shape(1:e%pnode,e%igaus)!- rho*(enthalpy+kineticEnerg)*alpha
      Atee(1:e%pnode) = Atee(1:e%pnode) + gpden*accph*e%shape(1:e%pnode,e%igaus) !rho*cp

   end subroutine   

   !-------------------------------------------------------------------

   subroutine DeallocateTransientCoefficients

      call a%Memor%dealloc(e%mnode,Atdd,'Atdd','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Atde,'Atde','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Atmd,'Atmd','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Atmm,'Atmm','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Atme,'Atme','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Ated,'Ated','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,Atem,'Atem','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,Atee,'Atee','nsc_pr_elmope')

   end subroutine

end module
