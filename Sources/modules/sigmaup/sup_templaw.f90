subroutine sup_templaw_WLF(e, gptemp, ReferenceTemp,c1_WLF, c2_WLF, LawViParam, acvis, lambda)
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement) :: e
   real(rp),    intent(in)  :: gptemp,ReferenceTemp,c1_WLF,c2_WLF, LawViParam(*)
   real(rp),    intent(out) :: acvis, lambda
   real(rp) :: lambda0, visc0,tempe0, c1, c2, func,diftemp,func2
   
   tempe0=ReferenceTemp
   c1=c1_WLF
   c2=c2_WLF

   visc0=LawViParam(1)
   lambda0=LawViParam(3)
   
   diftemp=gptemp-tempe0

   
   if (diftemp>-4) then
      func= exp(-c1*(diftemp)/(c2+diftemp))
   else
      func= exp(-c1*(-4)/(c2-4))
   end if
   
   
   acvis=visc0*func
   lambda=lambda0*func
     
end subroutine



subroutine sup_templaw_Arrhenius(e, gptemp, tempe0, alpha, LawViParam, acvis, lambda)
   use typre
   use Mod_Element
   implicit none
   class(FiniteElement) :: e
   real(rp),    intent(in)  :: gptemp,tempe0,alpha, LawViParam(*)
   real(rp),    intent(out) :: acvis, lambda
   real(rp)                 :: lambda_ini, visc0, func, diftemp, diftemp0
   

   visc0=LawViParam(1)
   lambda_ini=LawViParam(3)
   
   !Warning!! Temperature solution in ºC, not ºK.
   diftemp =gptemp+273.15 !=> to pass to Kelvin grades
   diftemp0=tempe0+273.15 !=> to pass to Kelvin grades

   
   if (diftemp==0 .or. diftemp0==0) then
      func= exp(alpha*(1.0_rp/(diftemp+0.00000001) -1.0_rp/(diftemp0+0.00000001)))
   else 
      func= exp(alpha*(1.0_rp/diftemp -1.0_rp/diftemp0))
   end if
   
   
   acvis=visc0*func
   lambda=lambda_ini*func
     
end subroutine