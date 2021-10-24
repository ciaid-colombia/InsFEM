 module Mod_supm_PostGauss
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersPostGauss
   
   integer(ip), allocatable :: kfl_IsSet
   
contains
   
   subroutine SetPointersPostGauss(itask)
      implicit none
      integer(ip) :: itask
      !procedure() :: NULL()
      select case (itask)   
      
      case(0)
         allocate(kfl_IsSet)
         call a%Memor%allocObj(0,'kfl_IsSet','InitProcedurePointer',1)
         kfl_IsSet = -1
      
      case(1)
         if (kfl_IsSet == -1) then
            kfl_IsSet = 1
             ProcPointer%PostGaussElmats_sup => PostGaussElmats
         endif  
         
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select   
   end subroutine   
   
!------------------------------------------------
!SUBROUTINES
   subroutine PostGaussElmats
      implicit none
      
      if(a%MatProp(imat)%lawvi>=0)then     
         ! tau2*(div v, div u)
         call nsm_elmdiv(e,dvolt2,elmuv)
                     
         !BLOCK P,Q : tau1*(graq q, grad p)
         call nsm_elmbpq(e,dvolt1,elmpq)
      end if
      
      !If you want the complete term div(e·grad u) uncomment next line 
      if(a%MatProp(imat)%lawvi<0)then
         auxvis=beta*acvis
          if ( a%fvins > zensi ) then
             call nsm_elmvis_div(e,dvolt0,auxvis,elmuv)
          endif
      end if
      
   end subroutine
   
end module
