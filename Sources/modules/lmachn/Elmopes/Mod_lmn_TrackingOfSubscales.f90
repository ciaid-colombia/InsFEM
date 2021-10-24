module Mod_lmn_ComputeSubscales
   use typre
   use Mod_lmn_BaseElmope
   use Mod_lmn_ComputeTaus
   use Mod_lmn_SubgridSpaceResidual
   implicit none
   private
   public SetPointersComputeSubscales,GPSGS, SetPointersGetSubscales

   !SubgridScales
   real(rp) :: GpSGS(5)

   integer(ip), allocatable :: kfl_IsSet, kfl_IsSetGetSubscales

contains

   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersComputeSubscales(itask)
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

            !We need to compute the residual and project it 
            !to the subgrid scale space
            call SetPointersComputeSubgridSpaceResidual(1) 
            call SetPointersComputeTaus(1)

            if (a%kfl_tacsg == 0) then

               call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleQSS)

            elseif ( a%kfl_tacsg == 1) then
               if (a%kfl_nolsg == 1) then
                  select case (a%kfl_nolsgScheme)
                  case (1)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleNonLinear_Pic_mom)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleNonLinear_Pic_ene)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScale_con)
                  case (2)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleNonLinear_Pic_mom)
                     call ConcatenateProcedures(ProcHook%InGaussElmats_seg,ComputeSubgridScaleNonLinear_Pic_ene)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScale_con)
                  case (3)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleNonLinear_NRSeg_mom)
                     call ConcatenateProcedures(ProcHook%InGaussElmats_seg,ComputeSubgridScaleNonLinear_NRSeg_ene)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScale_con)
                  case default
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleNonLinear_NR)
                     call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScale_con)
                  end select
               else
                  call ConcatenateProcedures(ProcHook%InGaussElmats,ComputeSubgridScaleDSS) 
               endif
            endif
         endif

      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)

      end select
   end subroutine   

   !----------------------------------------------------------------------------
   !Setting Pointers
   subroutine SetPointersGetSubscales(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   

      case(0)
         allocate(kfl_IsSetGetSubscales)
         call a%Memor%allocObj(0,'kfl_IsSetGetSubscales','InitProcedurePointer',1)
         kfl_IsSetGetSubscales = -1

      case(1)

         if (kfl_IsSetGetSubscales == -1) then
            kfl_IsSetGetSubscales = 1


            if (a%kfl_trasg == 0) then
               call SetPointersComputeSubscales(1)

            else 
               !Pressure subscale is not tracked, we need to recompute it
               !I need the subgrid scale residual
               call SetPointersComputeSubgridSpaceResidual(1)

               !We need to compute the Tau values
               call SetPointersComputeTaus(1)

               call ConcatenateProcedures(ProcHook%InGaussElmats,GetSubgridScaleDSS)
            endif
         endif  
      case(100)
         deallocate(kfl_IsSetGetSubscales)
         call a%Memor%deallocObj(0,'kfl_IsSetGetSubscales','InitProcedurePointer',1)
      end select  

   end subroutine   

   !------------------------------------------------------------------
   !Tracking of Subscales
   !Dynamic subscales

   !Computes the transient stabilization parameter
   subroutine ComputeSubgridScaleDSS
      implicit none

      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = (gpSubscaleSpaceResidual(1:e%ndime) + gpden(2)*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)*timom
      a%tesgs(ielem)%a(1,e%igaus) = (gpSubscaleSpaceResidual(e%ndime+2) + gpden(2)*a%tesgs(ielem)%a(2,e%igaus)*ReferenceDtinv)*tiene
      a%prsgs(ielem)%a(e%igaus) = ticon*gpSubscaleSpaceResidual(e%ndime+1)
      GpSGS(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(e%ndime+2) = a%tesgs(ielem)%a(1,e%igaus)
      GpSGS(e%ndime+1) = a%prsgs(ielem)%a(e%igaus)

   end subroutine

   !Static subscales
   subroutine ComputeSubgridScaleQSS
      implicit none

      GpSGS(1:e%ndime) = timom*gpSubscaleSpaceResidual(1:e%ndime)
      GpSGS(e%ndime+2) = tiene*gpSubscaleSpaceResidual(e%ndime+2)
      GpSGS(e%ndime+1) = ticon*gpSubscaleSpaceResidual(e%ndime+1)
      if (a%kfl_trasg /= 0) then
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = GpSGS(1:e%ndime)
         a%tesgs(ielem)%a(1,e%igaus) = GpSGS(e%ndime+2)
         a%prsgs(ielem)%a(e%igaus) = GpSGS(e%ndime+1)
      endif

   end subroutine

   subroutine ComputeSubgridScaleNonLinear_Pic_mom
      implicit none

      real(rp)    :: Matrix(e%ndime,e%ndime),InvMatrix(e%ndime,e%ndime),RHS(e%ndime),deter
      integer(ip) :: idime,jdime
      real(rp)    :: ResidualTerm(e%ndime)

      ResidualTerm(1:e%ndime) = (gpSubscaleSpaceResidual(1:e%ndime) + gpden(2)*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)

      !System Matrix
      Matrix = 0.0_rp

      !Implementation 1
      !Identity term
      do idime = 1,e%ndime
         Matrix(idime,idime) = 1.0_rp
      enddo

      !Tau_t times convective term
      do idime = 1,e%ndime
         do jdime = 1,e%ndime
            Matrix(idime,jdime) = Matrix(idime,jdime) + acden*timom*grvel(idime,jdime)
         enddo
      enddo

      !RHS
      RHS(1:e%ndime) = timom*ResidualTerm(1:e%ndime) - GpSGS(1:e%ndime)

      !Now we invert and solve
      call invmtx(Matrix,InvMatrix,deter,e%ndime)
      a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = matmul(InvMatrix,RHS) + GpSGS(1:e%ndime)

   end subroutine 

   subroutine ComputeSubgridScaleNonLinear_Pic_ene
      implicit none

      real(rp)    :: ResidualTerm

      ResidualTerm = (gpSubscaleSpaceResidual(e%ndime+2) + gpden(2)*a%tesgs(ielem)%a(2,e%igaus)*ReferenceDtinv) - acden*dot_product(a%vesgs(ielem)%a(:,1,e%igaus),grtem(1,:))
      a%tesgs(ielem)%a(1,e%igaus) = ResidualTerm*tiene

   end subroutine
      
   subroutine ComputeSubgridScale_con
      implicit none
      
      a%prsgs(ielem)%a(e%igaus) = ticon*gpSubscaleSpaceResidual(e%ndime+1)
   end subroutine

   subroutine ComputeSubgridScaleNonLinear_NRSeg_mom
      implicit none

      real(rp) :: Matrix(e%ndime,e%ndime),InvMatrix(e%ndime,e%ndime), RHS(e%ndime),deter
      integer(ip) :: idime,jdime
      real(rp) :: ResidualTerm(e%ndime),TauDeriv(e%ndime),Tau_1Deriv(e%ndime)

      if (gpvno > 0.0_rp) then
         ResidualTerm(1:e%ndime) = (gpSubscaleSpaceResidual(1:e%ndime) + gpden(2)*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
         Tau_1Deriv = acden*a%staco(2)/chale(1)*gpadv/gpvno
         TauDeriv = timom*Tau_1Deriv

         !System Matrix
         Matrix = 0.0_rp

         !Implementation 1
         !Identity term
         do idime = 1,e%ndime
            Matrix(idime,idime) = 1.0_rp
         enddo

         !Tau_t times convective term
         do idime = 1,e%ndime
            do jdime = 1,e%ndime
               Matrix(idime,jdime) = Matrix(idime,jdime) + acden*timom*grvel(idime,jdime)
            enddo
         enddo

         !Tau_t times residual term
         do idime = 1,e%ndime
            do jdime = 1,e%ndime
               Matrix(idime,jdime) = Matrix(idime,jdime) + TauDeriv(jdime)*GpSGS(idime)
            enddo
         enddo

         !RHS
         RHS(1:e%ndime) = timom*ResidualTerm(1:e%ndime) + dot_product(TauDeriv,GpSGS(1:e%ndime))*GpSGS(1:e%ndime) 

         !Now we invert and solve
         call invmtx(Matrix,InvMatrix,deter,e%ndime)
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = matmul(InvMatrix,RHS)
      end if   

   end subroutine

   subroutine ComputeSubgridScaleNonLinear_NRSeg_ene
      implicit none

      real(rp) :: ResidualTerm,aux,tau_aux

      aux = (a%staco(2)/chale(1)*gpvno*GpSGS(e%ndime+2) + dot_product(gpadv,grtem(1,:)) + ReferenceDtinv*(gptem(1)-gptem(2) + GpSGS(e%ndime+2))) * acden / (gptem(1)+GpSGS(e%ndime+2))
      ResidualTerm = (gpSubscaleSpaceResidual(e%ndime+2) + gpden(2)*a%tesgs(ielem)%a(2,e%igaus)*ReferenceDtinv) - acden*dot_product(a%vesgs(ielem)%a(:,1,e%igaus),grtem(1,:)) - aux*GpSGS(e%ndime+2)
      tau_aux = 1.0_rp / (acden*ReferenceDtinv + 1/tiene - aux)
      a%tesgs(ielem)%a(1,e%igaus) = ResidualTerm*tau_aux 

   end subroutine

   subroutine ComputeSubgridScaleNonLinear_NR
      implicit none

      real(rp) :: Matrix(e%ndime+1,e%ndime+1),InvMatrix(e%ndime+1,e%ndime+1),RHS(e%ndime+1),deter,vecsol(e%ndime+1)
      integer(ip) :: idime,jdime
      real(rp) :: ResidualTerm(e%ndime+1),Tau_1Deriv(e%ndime),TauDeriv_mm(e%ndime),TauDeriv_me(e%ndime),TauDeriv_em(e%ndime),TauDeriv_ee

      if (gpvno > 0.0_rp) then
         ResidualTerm(1:e%ndime) = (gpSubscaleSpaceResidual(1:e%ndime) + gpden(2)*a%vesgs(ielem)%a(1:e%ndime,2,e%igaus)*ReferenceDtinv)
         ResidualTerm(e%ndime+1) = (gpSubscaleSpaceResidual(e%ndime+2) + gpden(2)*a%tesgs(ielem)%a(2,e%igaus)*ReferenceDtinv)
         Tau_1Deriv  = acden*a%staco(2)/chale(1)*gpadv/gpvno
         TauDeriv_mm = timom*Tau_1Deriv
         TauDeriv_me = tiene*Tau_1Deriv
         TauDeriv_ee = (a%staco(2)/chale(1)*gpvno*GpSGS(e%ndime+2) + dot_product(gpadv,grtem(1,:)) + ReferenceDtinv*(gptem(1)-gptem(2) + GpSGS(e%ndime+2))) * acden* tiene / (gptem(1)+GpSGS(e%ndime+2))
         TauDeriv_em(1:e%ndime) = (a%staco(2)/chale(1)*gpvno*GpSGS(1:e%ndime) + matmul(gpadv,grvel) - a%gravi(1:e%ndime) + ReferenceDtinv*(gpvel(1:e%ndime,1)-gpvel(1:e%ndime,2) + GpSGS(1:e%ndime))) * acden*timom / (gptem(1)+GpSGS(e%ndime+2))

         !System Matrix
         Matrix = 0.0_rp

         !Implementation 1
         !Identity term
         do idime = 1,e%ndime+1
            Matrix(idime,idime) = 1.0_rp
         enddo

         !Tau_t times convective term
         do idime = 1,e%ndime
            do jdime = 1,e%ndime
               Matrix(idime,jdime) = Matrix(idime,jdime) + acden*timom*grvel(idime,jdime)
            enddo
            Matrix(e%ndime+1,idime) = Matrix(e%ndime+1,idime) + acden*tiene*grtem(1,idime)
         enddo

         !Tau_t times residual term
         do idime = 1,e%ndime
            do jdime = 1,e%ndime
               Matrix(idime,jdime) = Matrix(idime,jdime) + TauDeriv_mm(jdime)*GpSGS(idime)
            enddo
            Matrix(e%ndime+1,idime) = Matrix(e%ndime+1,idime) + TauDeriv_me(idime)*GpSGS(e%ndime+2)
            Matrix(idime,e%ndime+1) = Matrix(idime,e%ndime+1) - TauDeriv_em(idime)
         enddo
         Matrix(e%ndime+1,e%ndime+1) = Matrix(e%ndime+1,e%ndime+1) - TauDeriv_ee

         !RHS
         RHS(1:e%ndime) = timom*ResidualTerm(1:e%ndime) + dot_product(TauDeriv_mm,GpSGS(1:e%ndime))*GpSGS(1:e%ndime) - GpSGS(e%ndime+2)*TauDeriv_em(1:e%ndime)
         RHS(e%ndime+1) = tiene*ResidualTerm(e%ndime+1) - GpSGS(e%ndime+2)*TauDeriv_ee + dot_product(TauDeriv_me,GpSGS(1:e%ndime))*GpSGS(e%ndime+2)

         !Now we invert and solve
         call invmtx(Matrix,InvMatrix,deter,e%ndime+1)
         vecsol = matmul(InvMatrix,RHS)
         a%vesgs(ielem)%a(1:e%ndime,1,e%igaus) = vecsol(1:e%ndime)
         a%tesgs(ielem)%a(1,e%igaus) = vecsol(e%ndime+1)
      end if   

   end subroutine

   subroutine GetSubgridScaleDSS
      implicit none

      GpSGS(1:e%ndime) = a%vesgs(ielem)%a(1:e%ndime,1,e%igaus)
      GpSGS(e%ndime+2) = a%tesgs(ielem)%a(1,e%igaus)
      GpSGS(e%ndime+1) = ticon*gpSubscaleSpaceResidual(e%ndime+1)
   end subroutine
end module 
