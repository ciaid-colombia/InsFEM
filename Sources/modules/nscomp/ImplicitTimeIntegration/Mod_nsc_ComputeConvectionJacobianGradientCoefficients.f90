module Mod_nsc_ComputeConvectionJacobianGradientCoefficients
   use typre
   use Mod_nsc_BaseElmope
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
               call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateJacobianGradientCoefficients)
               call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,JacobianGradientCoefficientsToZero)
               call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,ComputeJacobianGradientCoefficients)
               call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateJacobianGradientCoefficients)
               if (a%lawde == 1) then
                  call ConcatenateProcedures(ProcPointer_nsc_ComputeConvectionCoefficients,ComputeIdealJacobianGradientCoefficients)
               end if
            endif


         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateJacobianGradientCoefficients

      call a%Memor%alloc(e%ndime,e%mnode,JAmd,'JAmd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,JAmm,'JAmm','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,JAed,'JAed','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,JAem,'JAem','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,JAee,'JAee','nsc_elmope_im')

   end subroutine

   !JacobianGradient matrix coefficients
   subroutine JacobianGradientCoefficientsToZero
      implicit none    

      JAmd = 0.0_rp !JAmd(i,p)
      JAmm = 0.0_rp !JAmm(i,d,p)
      JAed = 0.0_rp !JAed(p)
      JAem = 0.0_rp !JAem(d,p)
      JAee = 0.0_rp !JAee(p)
 
   end subroutine   

   subroutine ComputeJacobianGradientCoefficients
      implicit none    
      integer(ip)                :: idime,jdime

      !Momentum equation
      do idime=1,e%ndime
         JAmd(idime,:) = JAmd(idime,:) + 2*invgpd*gpadv(idime)*vgden*e%shape(:,e%igaus)!(2/rho)*vel_i*vel_j*d_j(rho) 
         JAmd(idime,:) = JAmd(idime,:) - invgpd*gpadv(idime)*divmom*e%shape(:,e%igaus)!-(1/rho)*vel_i*d_j(mom_j) 
         JAmd(idime,:) = JAmd(idime,:) - invgpd*vgmom(idime)*e%shape(:,e%igaus)!-(1/rho)*vel_j*d_j(mom_i) 

         JAmm(idime,idime,:) = JAmm(idime,idime,:) - invgpd*vgden*e%shape(:,e%igaus)!-(1/rho)*(vel_j)*d_j(rho)*d_id 
         JAmm(idime,idime,:) = JAmm(idime,idime,:) + invgpd*divmom*e%shape(:,e%igaus)!-(1/rho)*d_j(mom_j)*d_id 
         do jdime=1,e%ndime
            JAmm(idime,jdime,:) = JAmm(idime,jdime,:) - invgpd*gpadv(idime)*grden(jdime)*e%shape(:,e%igaus)!-(1/rho)*(vel_i)*d_d(rho) 
            JAmm(idime,jdime,:) = JAmm(idime,jdime,:) + invgpd*grmom(idime,jdime)*e%shape(:,e%igaus)!-(1/rho)*d_d(mom_i) 
         end do

          !Energy equation
          JAem(idime,:) = JAem(idime,:) - (gpade+gppre)*sqinvgpd*grden(idime)*e%shape(:,e%igaus) !((e+p)/rho^2)*d_d(rho)   
      end do
      JAed(:) = JAed(:) + 2*(gpade+gppre)*sqinvgpd*vgden*e%shape(:,e%igaus) !2*(e+p)/rho^2 vel_j*d_j(rho)   
      JAed(:) = JAed(:) - (gpade+gppre)*sqinvgpd*divmom *e%shape(:,e%igaus) !-(e+p)/rho^2 *d_j(mom_j)   

   end subroutine   

   !-------------------------------------------------------------------
   !Ideal state law matrix coefficients
   subroutine ComputeIdealJacobianGradientCoefficients
      implicit none
      integer(ip)                :: idime,jdime
      
      !Momentum equation
      do idime=1,e%ndime
         JAmd(idime,:) = JAmd(idime,:) - 2*aux_d*invgpd*grden(idime)*e%shape(:,e%igaus)!-((gamma-1)/rho)*(vel·vel)*d_i(rho)    
         JAmd(idime,:) = JAmd(idime,:) + aux*invgpd*dot_product(gpadv(:),grmom(:,idime))*e%shape(:,e%igaus)!((gamma-1)/rho)*vel_j*d_i(mom_j)    

         do jdime=1,e%ndime
            JAmm(idime,jdime,:) = JAmm(idime,jdime,:) + aux*invgpd*gpadv(jdime)*grden(idime)*e%shape(:,e%igaus)!(gamma-1/rho)*(vel_d)*d_i(rho) 
            JAmm(idime,jdime,:) = JAmm(idime,jdime,:) - aux*invgpd*grmom(jdime,idime)*e%shape(:,e%igaus)!-(gamma-1/rho)*d_i(mom_d) 
         end do

      !Energy equation
          JAem(idime,:) = JAem(idime,:) + aux_d*invgpd*grden(idime)*e%shape(:,e%igaus) !((gamma-1)/rho)*(vel·vel/2)*d_d(rho)   
          JAem(idime,:) = JAem(idime,:) + 2*aux*invgpd*vgden*gpadv(idime)*e%shape(:,e%igaus) !((gamma-1)/rho)*2*vel_d*vel_j*d_j(rho)   
          JAem(idime,:) = JAem(idime,:) - aux*invgpd*dot_product(gpadv(:),grmom(:,idime))*e%shape(:,e%igaus) !-((gamma-1)/rho)*vel_j*d_d(mom_j)   
          JAem(idime,:) = JAem(idime,:) - aux*invgpd*vgmom(idime)*e%shape(:,e%igaus) !-((gamma-1)/rho)*vel_j*d_j(mom_d)   
          JAem(idime,:) = JAem(idime,:) - aux*invgpd*divmom*gpadv(idime)*e%shape(:,e%igaus) !-((gamma-1)/rho)*vel_d*d_j(mom_j)  
          JAem(idime,:) = JAem(idime,:) + aux*invgpd*grene(idime)*e%shape(:,e%igaus) !-((gamma-1)/rho)*d_d(ene)  
      end do
      JAed(:) = JAed(:) + 4*aux_d*invgpd*vgden*e%shape(:,e%igaus) !2*((gamma-1)/rho)*(vel·vel)*vel_j*d_j(rho)   
      JAed(:) = JAed(:) + 2*aux*invgpd*dot_product(gpadv,vgmom)*e%shape(:,e%igaus) !2*((gamma-1)/rho)*vel_i*vel_j*d_j(mom_i)   
      JAed(:) = JAed(:) + 2*aux_d*invgpd*divmom*e%shape(:,e%igaus) !((gamma-1)/rho)*(vel·vel)*d_j(mom_j)   
      JAed(:) = JAed(:) - acgamma*invgpd*vgene*e%shape(:,e%igaus) !-(gamma/rho)*vel_j*d_j(ene)   

      JAee(:) = JAee(:) + acgamma*invgpd*(divmom-vgden)*e%shape(:,e%igaus) !(gamma/rho)·(d_j(mom_j) - vel_j*d_j(rho))  

   end subroutine   

   subroutine DeallocateJacobianGradientCoefficients

      call a%Memor%dealloc(e%ndime,e%mnode,JAmd,'JAmd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,JAmm,'JAmm','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,JAed,'JAed','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,JAem,'JAem','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,JAee,'JAee','nsc_elmope_im')

   end subroutine

end module
