module Mod_nsc_ComputeReactionCoefficients
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersReactionCoefficients
   
   integer(ip), allocatable :: kfl_IsSet
   integer(ip) ::  kfl_damp 
   real(rp)    :: dampc


contains

   subroutine SetPointersReactionCoefficients(itask)
      implicit none
      integer(ip) :: itask
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
         kfl_damp = 0
      
      case(1)

         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
         
            call ConcatenateProcedures(ProcHook_nsc_Initializations,AllocateReactionCoefficients)
            call ConcatenateProcedures(ProcHook_nsc_Finalizations,DeallocateReactionCoefficients)
            call ConcatenateProcedures(ProcPointer_nsc_ComputeTransportCoefficients,ReactionCoefficientsToZero)
            !Sources Exists
            if (a%kfl_sourc == 1) then
                call ConcatenateProcedures(ProcPointer_nsc_ComputeTransportCoefficients,ComputeSourcesCoefficients)
            endif
            !Rectangular Damping Exists
            if (a%ndamp > 0) then
               kfl_damp = 1
               call ConcatenateProcedures(ProcPointer_nsc_ComputeTransportCoefficients,ComputeRectangularDamping)
            endif

            !Radial Damping Exists
            if (a%rdamp > 0) then
               kfl_damp = 1
               call ConcatenateProcedures(ProcPointer_nsc_ComputeTransportCoefficients,ComputeRadialDamping)
            endif

            if (kfl_damp == 1) then
                call ConcatenateProcedures(ProcPointer_nsc_ComputeTransportCoefficients,ComputeDampingCoefficients)
            endif

         endif
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine
   
   
   !-----------------------------------------------------------------------
   subroutine AllocateReactionCoefficients

      call a%Memor%alloc(e%mnode,Sdd,'Sdd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Smd,'Smd','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%ndime,e%mnode,Smm,'Smm','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,Sed,'Sed','nsc_elmope_im')
      call a%Memor%alloc(e%ndime,e%mnode,Sem,'Sem','nsc_elmope_im')
      call a%Memor%alloc(e%mnode,See,'See','nsc_elmope_im')

   end subroutine

   subroutine ReactionCoefficientsToZero

      Sdd = 0.0_rp !Smd(p)
      Smd = 0.0_rp !Smd(i,p)
      Smm = 0.0_rp !Smm(i,d,p)
      Sem = 0.0_rp !Sem(d,p)
      Sed = 0.0_rp !Sed(p)
      See = 0.0_rp !See(p)

   end subroutine   

   subroutine ComputeSourcesCoefficients
     implicit none
     integer(ip) :: idime

      !Momentum equation
      do idime=1,e%ndime
         Smd(idime,:) = Smd(idime,:) + a%gravi(idime)*a%grnor*e%shape(:,e%igaus)!f
      !Energy equation
         Sem(idime,:) = Sem(idime,:) + a%gravi(idime)*a%grnor*e%shape(:,e%igaus)!f
      end do
      Sed(:) = Sed(:) + a%srce*e%shape(:,e%igaus) !r   

   end subroutine   

   subroutine ComputeRectangularDamping
     use Mod_Mesh    
     implicit none      
     real(rp)    :: gpcod(e%ndime)

     !Interpolate
     call e%interpg(e%ndime,e%elcod,gpcod)

     !Compute coefficient of damping
     call damping%nsc_ComputeRectangularDamping(e%ndime,gpcod,a)
     call damping%nsc_GetDampingCoefficient(dampc)
     
   end subroutine   

   subroutine ComputeRadialDamping
     use Mod_Mesh    
     implicit none      
     real(rp)    :: gpcod(e%ndime)

     !Interpolate
     call e%interpg(e%ndime,e%elcod,gpcod)

     !Compute coefficient of damping
     call damping%nsc_ComputeRadialDamping(e%ndime,gpcod,a)
     call damping%nsc_GetDampingCoefficient(dampc)
     
   end subroutine   

   subroutine ComputeDampingCoefficients
     implicit none      
     integer(ip) :: idime
     
     !Mass equation
     Sdd(:) = Sdd(:) - dampc*e%shape(:,e%igaus)*LHSdtinv
     !Momentum equation
     do idime=1,e%ndime
        Smm(idime,idime,:) = Smm(idime,idime,:) - dampc*e%shape(:,e%igaus)*LHSdtinv
     end do
     !Energy equation
     See(:) = See(:) - dampc*e%shape(:,e%igaus)*LHSdtinv

   end subroutine   

   subroutine DeallocateReactionCoefficients

      call a%Memor%dealloc(e%mnode,Sdd,'Sdd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Smd,'Smd','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%ndime,e%mnode,Smm,'Smm','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,Sed,'Sed','nsc_elmope_im')
      call a%Memor%dealloc(e%ndime,e%mnode,Sem,'Sem','nsc_elmope_im')
      call a%Memor%dealloc(e%mnode,See,'See','nsc_elmope_im')

   end subroutine

end module
