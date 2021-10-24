module Mod_supm_PressTempSubscale
   use typre
   use Mod_supm_BaseElmope
   implicit none
   private
   public SetPointersPressTempSubscale
   
   integer(ip), allocatable :: kfl_IsSet

contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SetPointersPressTempSubscale(itask)
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
         
            !Tracking of the subscales
            !Dynamic subscales- Pressure
            if (a%kfl_PressTempSubscale== 1) then
               call ConcatenateProcedures(ProcPointer%PostGaussElmats_sup,PostPressTempSubscale)
            endif
            
         endif  
      case(100)
         deallocate(kfl_IsSet)
         call a%Memor%deallocObj(0,'kfl_IsSet','InitProcedurePointer',1)
      
      end select
   end subroutine   

   !-------------------------------------------------------------------
   !Computations subroutines
   subroutine PostPressTempSubscale
      implicit none
      ! tau2*(div v, div u)
      call nsm_elmdiv(e,dvolt2*(dvolt1/a%dtinv),elmuv)
      
      call supm_elmrhu_dynamicPressSubscale(e,dvolt2*(dvolt1/a%dtinv),grvel,elrhu)
      
   end subroutine 
   
      
   subroutine supm_elmrhu_dynamicPressSubscale(e,dvol,grvel,elrhs)
      implicit none

      class(FiniteElement) :: e
      real(rp),    intent(in)    :: dvol
      real(rp),    intent(in)    :: grvel(e%ndime,e%ndime)
      real(rp),    intent(inout) :: elrhs(e%ndime,e%pnode)

      integer(ip)                :: inode,idime

      do inode=1,e%pnode
         do idime=1,e%ndime
            elrhs(:,inode) =  dvol*e%cartd(:,inode)*grvel(idime,idime) + elrhs(:,inode)
         end do   
      end do

   end subroutine supm_elmrhu_dynamicPressSubscale
   
end module
