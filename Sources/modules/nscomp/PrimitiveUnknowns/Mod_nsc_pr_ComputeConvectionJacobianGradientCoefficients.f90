module Mod_nsc_pr_ComputeConvectionJacobianGradientCoefficients
   use typre
   use Mod_nsc_pr_BaseElmope
   implicit none
   private
   public SetPointersConvectionJacobianGradientCoefficients
   
   !External Procedures
   procedure() :: NULLSUB
   
   integer(ip), allocatable :: kfl_IsSet
   
contains

   subroutine SetPointersConvectionJacobianGradientCoefficients(itask)
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
         
            !Convection Jacobian Gradient Matrix Coefficients
            if (a%kfl_jacgr == 1) then
               call ConcatenateProcedures(ProcHook_nsc_pr_Initializations,AllocateJacobianGradientCoefficients)
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeConvectionCoefficients,JacobianGradientCoefficientsToZero)
               call ConcatenateProcedures(ProcPointer_nsc_pr_ComputeConvectionCoefficients,ComputeJacobianGradientCoefficients)
               call ConcatenateProcedures(ProcHook_nsc_pr_Finalizations,DeallocateJacobianGradientCoefficients)
            endif


         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateJacobianGradientCoefficients

      call a%Memor%alloc(e%mnode,JAdd,'JAdd','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,JAdm,'JAdm','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,JAde,'JAde','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,JAmd,'JAmd','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,JAmm,'JAmm','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,JAme,'JAme','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,JAed,'JAed','nsc_pr_elmope')
      call a%Memor%alloc(e%ndime,e%mnode,JAem,'JAem','nsc_pr_elmope')
      call a%Memor%alloc(e%mnode,JAee,'JAee','nsc_pr_elmope')

   end subroutine

   !JacobianGradient matrix coefficients
   subroutine JacobianGradientCoefficientsToZero
      implicit none    

      JAdd = 0.0_rp !JAdd(p)
      JAdm = 0.0_rp !JAdm(d,p)
      JAde = 0.0_rp !JAde(p)
      JAmd = 0.0_rp !JAmd(i,p)
      JAmm = 0.0_rp !JAmm(i,d,p)
      JAme = 0.0_rp !JAme(i,p)
      JAed = 0.0_rp !JAed(p)
      JAem = 0.0_rp !JAem(d,p)
      JAee = 0.0_rp !JAee(p)
 
   end subroutine   

   subroutine ComputeJacobianGradientCoefficients
      implicit none    
      integer(ip)                :: idime,jdime
      real(rp)                   :: grden(e%ndime),grinvden(e%ndime),vgden

      grden = (grpre/gpadp)-(grtem/gpadt)
      grinvden = (grtem/(gpden*gpadt))-(grpre/(gpden*gpadp))
      vgden = (vgpre/gpadp)-(vgtem/gpadt)

      !Mass equation
      JAdd(:) = JAdd(:) + gpden*acbeta*vgden*e%shape(:,e%igaus) !rho*beta*vel_j*(d_y rho * d_j Y)*V 
      JAdd(:) = JAdd(:) + gpden*acbeta*divvel*e%shape(:,e%igaus) !rho*beta*(d_j vel_j)*V 
      JAde(:) = JAde(:) - gpden*acalpha*vgden*e%shape(:,e%igaus) !-rho*alpha*vel_j*(d_y rho * d_j Y)*V 
      JAde(:) = JAde(:) - gpden*acalpha*divvel*e%shape(:,e%igaus) !-rho*alpha*(d_j vel_j)*V 
      do idime=1,e%ndime
         JAdm(idime,:) = JAdm(idime,:) + gpden*grden(idime)*e%shape(:,e%igaus) !rho*d_dj*(d_Y rho * d_j Y)*V 
      !Momentum equation
         JAmd(idime,:) = JAmd(idime,:) + gpden*acbeta*gpadv(idime)*vgden*e%shape(:,e%igaus)!rho*beta*vel_i*vel_j*(d_y rho * d_j Y)*V
         JAmd(idime,:) = JAmd(idime,:) + gpden*acbeta*gpadv(idime)*divvel*e%shape(:,e%igaus)!rho*beta*vel_i*(d_j vel_j)*V 
         JAmd(idime,:) = JAmd(idime,:) + gpden*acbeta*vgvel(idime)*e%shape(:,e%igaus)!rho*beta*vel_j*(d_j vel_i)*V 
         JAmd(idime,:) = JAmd(idime,:) + grden(idime)*e%shape(:,e%igaus) !d_ij*(d_Y rho * d_j Y)*V 
         JAmd(idime,:) = JAmd(idime,:) + grinvden(idime)*e%shape(:,e%igaus) !d_ij*(d_Y 1/rho * d_j Y)*V 

         JAmm(idime,idime,:) = JAmm(idime,idime,:) + gpden*vgden*e%shape(:,e%igaus)!rho*(vel_j)*d_id*(d_y rho * d_j Y)*V
         JAmm(idime,idime,:) = JAmm(idime,idime,:) + gpden*divvel*e%shape(:,e%igaus)!rho*d_j(vel_j)*d_id*V 
         do jdime=1,e%ndime
            JAmm(idime,jdime,:) = JAmm(idime,jdime,:) + gpden*gpadv(idime)*grden(jdime)*e%shape(:,e%igaus)!rho*(vel_i)*d_dj*(d_Y rho * d_j Y)*V
            JAmm(idime,jdime,:) = JAmm(idime,jdime,:) + gpden*grvel(idime,jdime)*e%shape(:,e%igaus)!rho*d_dj*d_j(vel_i)*V 
         end do
         JAme(idime,:) = Ame(idime,:) - gpden*acalpha*gpadv(idime)*vgden*e%shape(:,e%igaus)!-rho*alpha*vel_i*vel_j*(d_Y rho * d_j Y)*V 
         JAme(idime,:) = Ame(idime,:) - gpden*acalpha*gpadv(idime)*divvel*e%shape(:,e%igaus)!-rho*alpha*vel_i*(d_j vel_j)*V 
         JAme(idime,:) = Ame(idime,:) - gpden*acalpha*vgvel(idime)*e%shape(:,e%igaus)!-rho*alpha*vel_j*(d_j vel_i)*V 
          !Energy equation
         JAem(idime,:) = JAem(idime,:) + gpden*gpadv(idime)*vgden*e%shape(:,e%igaus)!rho*vel_d*vel_j*(d_Y rho * d_j Y)*V   
         JAem(idime,:) = JAem(idime,:) + gpden*gpadv(idime)*divvel*e%shape(:,e%igaus)!rho*vel_d*(d_j vel_j)*V   
         JAem(idime,:) = JAem(idime,:) + gpden*vgvel(idime)*e%shape(:,e%igaus)!rho*vel_j*(d_j vel_d)*V   
         JAem(idime,:) = JAem(idime,:) + aux*grden(idime)*e%shape(:,e%igaus) !rho*h+k*d_dj*(d_Y rho * d_j Y)*V   
      end do
      JAed(:) = JAed(:) + aux*acbeta*vgden*e%shape(:,e%igaus) !rho*h+k*beta*vel_j*(d_Y rho * d_j Y)*V   
      JAed(:) = JAed(:) + aux*acbeta*divvel*e%shape(:,e%igaus) !rho*h+k*beta*(d_j vel_j)*V   
      JAed(:) = JAed(:) - acalpha*gpadt*vgden*e%shape(:,e%igaus) !-alpha*temperature*vel_j*(d_Y rho * d_j Y)*V
      JAed(:) = JAed(:) - acalpha*gpadt*divvel*e%shape(:,e%igaus) !-alpha*temperature*(d_j vel_j)*V
      JAed(:) = JAed(:) - acalpha*vgtem*e%shape(:,e%igaus) !-alpha*vel_j*(d_j temperature)*V
      JAed(:) = JAed(:) + acalpha*gpadt*vgpre*e%shape(:,e%igaus)/gpadp !alpha*(temperature/pressure)*vel_j*(d_j press)*V
      JAed(:) = JAed(:) - acalpha*vgtem*e%shape(:,e%igaus) !alpha*vel_j*(d_j tempe)*V
      JAed(:) = JAed(:) + divvel*e%shape(:,e%igaus) !(d_j vel_j)*V   
      JAed(:) = JAed(:) - vgpre*e%shape(:,e%igaus)/gpadp !-(1/pressure)*vel_j*(d_j press)*V
      JAed(:) = JAed(:) + vgtem*e%shape(:,e%igaus)/gpadt !(1/tempe)*vel_j*(d_j tempe)*V

      JAee(:) = JAee(:) - acalpha*aux*vgden*e%shape(:,e%igaus) !-rho*alpha*h+k*vel_j*(d_Y rho * d_j Y)*V   
      JAee(:) = JAee(:) - acalpha*aux*divvel*e%shape(:,e%igaus) !-rho*alpha*h+k*(d_j vel_j)*V   
      JAee(:) = JAee(:) + gpden*accph*vgden*e%shape(:,e%igaus) !rho*accph*vel_j*(d_Y rho * d_j Y)*V   
      JAee(:) = JAee(:) + gpden*accph*divvel*e%shape(:,e%igaus) !rho*accph*(d_j vel_j)*V   

   end subroutine   

   subroutine DeallocateJacobianGradientCoefficients

      call a%Memor%dealloc(e%mnode,JAdd,'JAdd','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,JAdm,'JAdm','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,JAde,'JAde','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,JAmd,'JAmd','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,JAmm,'JAmm','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,JAme,'JAme','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,JAed,'JAed','nsc_pr_elmope')
      call a%Memor%dealloc(e%ndime,e%mnode,JAem,'JAem','nsc_pr_elmope')
      call a%Memor%dealloc(e%mnode,JAee,'JAee','nsc_pr_elmope')

   end subroutine

end module
