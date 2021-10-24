 module Mod_nsc_ExternalForces
   use typre
   use Mod_nsc_BaseElmope
   implicit none
   private
   public SetPointersExternalForces

   integer(ip), allocatable :: kfl_IsSet 
   integer(ip) ::  kfl_damp 
   

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
 
   subroutine SetPointersExternalForces(itask)
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
            
            !ExactSol
            if (a%kfl_exacs/=0) then  
               call ConcatenateProcedures(ProcPointer_nsc_ExternalForces,ExactSolutionForces)
            endif 

            !Rectangular Damping Exists
            if (a%ndamp > 0) then
               kfl_damp = 1
               call ConcatenateProcedures(ProcPointer_nsc_ExternalForces,RectangularDamping)
            endif

            !Radial Damping Exists
            if (a%rdamp > 0) then
               kfl_damp = 1
               call ConcatenateProcedures(ProcPointer_nsc_ExternalForces,RadialDamping)
            endif

            if (kfl_damp == 1) then
               call ConcatenateProcedures(ProcPointer_nsc_ExternalForces,DampingForces)
            endif

         endif  
      
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   
   
   
   
   !---------------------------------------------------------------------------
   !Computation Subroutines

   !Exact solution forces 
   subroutine  ExactSolutionForces    
      use Mod_Mesh    
      implicit none      
      real(rp)    :: gpcod(e%ndime)
      !Interpolate
      call e%interpg(e%ndime,e%elcod,gpcod)

      !Compute vector of external forces  

      call exacso%nsc_ComputeSolution(e%ndime,gpcod,a)
      call exacso%nsc_GetConservativeForce(e%ndime,elexd,elexm,elexe,a)        
    
   end subroutine

   !Rectangular Damping
   subroutine  RectangularDamping
      use Mod_Mesh    
      implicit none      
      real(rp)    :: gpcod(e%ndime)
      !Interpolate
      call e%interpg(e%ndime,e%elcod,gpcod)

      !Compute damping  

      call damping%nsc_ComputeRectangularDamping(e%ndime,gpcod,a)
    
   end subroutine

   !Radial Damping forces 
   subroutine  RadialDamping
      use Mod_Mesh    
      implicit none      
      real(rp)    :: gpcod(e%ndime)
      !Interpolate
      call e%interpg(e%ndime,e%elcod,gpcod)

      !Compute damping  

      call damping%nsc_ComputeRadialDamping(e%ndime,gpcod,a)
    
   end subroutine

   !Damping forces 
   subroutine  DampingForces    
      implicit none      
      !Compute vector of damping forces  

      call damping%nsc_GetDampingForce(e%ndime,a%dtinv,elexd,elexm,elexe)        
    
   end subroutine

end module

