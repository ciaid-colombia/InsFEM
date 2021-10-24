module Mod_nsm_Viscosity
   use typre
   use Mod_Element 
contains
   
   subroutine nsi_vislaw(e,grvel,lawvi,LawViParam,acvis)
      implicit none
      class(FiniteElement)     :: e
      integer(ip), intent(in)  :: lawvi
      real(rp),    intent(in)  :: LawViParam(*)
      real(rp),    intent(in)  :: grvel(e%ndime,e%ndime)
      real(rp),    intent(out) :: acvis
      integer(ip) :: idime,jdime,inode,vislaw 
      real(rp)    :: vis_a,vispolaw,polawn,viszer,visinf,vislamb
      real(rp)    :: I2,strain(e%ndime,e%ndime),facvis1,facvis2,heavi
   
      !********************************************************************
      !vislawParam(1)=newtonian viscosity (mu) or consistency index (m) 
      !visco(2)=infinity viscosity Carreau-Yasuda model
      !visco(3)=power law index in non newtonian models
      !visco(4)=lambda parameter in Carreau Yasuda model
      !visco(5)=a parameter in carreau yasuda model if a=2 carreau model
      !********************************************************************
   
      vispolaw=LawViParam(1)
      polawn=LawViParam(2)
      viszer=LawViParam(1)
      visinf=LawViParam(3)
      vislamb=LawViParam(4)
      vis_a=LawViParam(5)
      vislaw=lawvi
      
      !rate of strain tensor to use in non Newtonian models
      I2 = 0.0_rp
      strain = grvel + transpose(grvel)
      call tracedoubledot(e%ndime,strain,I2)
      !call trace(e%ndime,strain,tr_strain1)
      !call trace(e%ndime,matmul(strain,strain),tr_strain2)
      !I2 = tr_strain1**2 + tr_strain2
      
      if (vislaw == 1) then
         acvis = vispolaw*((0.5_rp*I2)**(polawn*0.5_rp-0.5_rp)) !Power Law viscosity
         if (polawn .gt. 1.0_rp) then
            facvis1 = LawViParam(1)/10000_rp
            facvis2 = LawViParam(1)*10000_rp
         end if
         if (polawn .le. 1.0_rp) then
            facvis1 = LawViParam(1)/10000_rp
            facvis2 = LawViParam(1)*10000_rp
         end if
         if (acvis <= facvis1) acvis = facvis1
         if (acvis >= facvis2) acvis = facvis2
      elseif (vislaw == 2) then
         acvis = visinf+(viszer-visinf)*(1.0_rp+(0.5_rp*vislamb*I2)**vis_a)**((polawn-1.0_rp)/vis_a)
      end if
   
   end subroutine nsi_vislaw

   subroutine nsm_smago(e,grvel,denac,turbu,vista)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)  :: e
      real(rp), intent(in)  :: turbu
      real(rp), intent(in)  :: grvel(e%ndime,e%ndime)
      real(rp), intent(in)  :: denac
      real(rp), intent(out) :: vista
      real(rp)    :: seci4,divel
      integer(ip) :: igaus,inode,idime,jdime
   
      vista = 0.0_rp
      !Smagorinsky
      seci4=0.0_rp
      do idime=1,e%ndime
         do jdime=1,e%ndime
            seci4=seci4+grvel(idime,jdime)*(grvel(idime,jdime)+grvel(jdime,idime))
         end do
      end do
      divel=0.0_rp
      do idime=1,e%ndime
         divel=divel+grvel(idime,idime)
      end do
      vista=2.0_rp*denac*turbu*e%hleng(e%ndime)*e%hleng(e%ndime)*sqrt(seci4-2.0_rp*divel*divel/3.0_rp)

   end subroutine nsm_smago

   !> This is the WALE LES model based on  Subgrid-Scale Stress Modelling Based on the
   !>   Square of the Velocity Gradient Tensor. F. NICOUD and F. DUCROS Flow, Turbulence and Combustion 62: 183–200, 1999.
   subroutine nsm_wale(e,grvel,denac,turbu,vista)
      use typre
      use Mod_Element
      implicit none
      class(FiniteElement)  :: e
      real(rp), intent(in)  :: turbu
      real(rp), intent(in)  :: grvel(e%ndime,e%ndime)
      real(rp), intent(in)  :: denac
      real(rp), intent(out) :: vista
      real(rp) :: CaligraphicS(e%ndime,e%ndime)
      real(rp) :: S(e%ndime,e%ndime), g2(e%ndime,e%ndime)
      real(rp) :: MeanStressG2, OP1,OP2, TraceG2
      real(rp) :: CSDoubleContractionCS, SDoubleContractionS
      integer(ip) :: idime,jdime
      
      !Default value
      vista = 0.0_rp
      !Symmetric gradient
      S = 0.5_rp*(grvel + transpose(grvel))
      !Gradient times gradient
      g2 = matmul(grvel,grvel)
      
      TraceG2 = 0.0_rp
      do idime = 1,e%ndime
         TraceG2 = TraceG2 + g2(idime,idime)
      enddo
      MeanStressG2 = 1.0_rp/3.0_rp*TraceG2
      
      !CaligraphicS
      CaligraphicS = 0.5_rp*(g2 + transpose(g2)) 
      do idime = 1,e%ndime
         CaligraphicS(idime,idime) = CaligraphicS(idime,idime) - MeanStressG2
      enddo
      
      !Symmetric gradient against itself
      SDoubleContractionS = 0.0_rp
      do idime = 1,e%ndime
         do jdime = 1,e%ndime
            SDoubleContractionS = SDoubleContractionS + S(idime,jdime)*S(idime,jdime)
         enddo
      enddo
      
      !Caligraphic S against itself
      CSDoubleContractionCS = 0.0_rp
      do idime = 1,e%ndime
         do jdime = 1,e%ndime
            CSDoubleContractionCS = CSDoubleContractionCS + CaligraphicS(idime,jdime)*CaligraphicS(idime,jdime)
         enddo
      enddo
      
      !OP1 and OP2 as defined in the Nicoud paper
      OP1 = CSDoubleContractionCS**(3.0_rp/2.0_rp)
      OP2 = SDoubleContractionS**(5.0_rp/2.0_rp) + CSDoubleContractionCS**(5.0_rp/4.0_rp)
      
      if (OP2 /= 0.0_rp) vista = denac*((turbu*e%hleng(e%ndime))**2.0_rp)*OP1/OP2
      
   end subroutine nsm_wale

end module
