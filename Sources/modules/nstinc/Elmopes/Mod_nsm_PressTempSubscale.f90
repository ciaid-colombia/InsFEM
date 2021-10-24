module Mod_nsm_PressTempSubscale
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersPressTempSubscale

   type, extends(PointerSetter) :: SPPressTempSubscale
contains
      procedure :: SpecificSet => SpecificSetPressTempSubscale
   end type
   type(SPPressTempSubscale) :: SetPointersPressTempSubscale
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetPressTempSubscale(d)
      implicit none
      class(SPPressTempSubscale) :: d
         
      !Tracking of the subscales
      !Dynamic subscales- Pressure
      if (a%kfl_PressTempSubscale== 1) then
         call ConcatenateProcedures(ProcPointer_PostGaussElmats,PostPressTempSubscale)
      endif
            
   end subroutine   

   !-------------------------------------------------------------------
   !Computations subroutines
   subroutine PostPressTempSubscale
      implicit none
      ! tau2*(div v, div u)
      call nsm_elmdiv(e,dvolt2*(dvolt1/a%dtinv),elmuv)
      
      call nsm_elmrhu_dynamicPressSubscale(e,dvolt2*(dvolt1/a%dtinv),grvel,elrhu)
      
   end subroutine 
   
      
   subroutine nsm_elmrhu_dynamicPressSubscale(e,dvol,grvel,elrhs)
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

   end subroutine nsm_elmrhu_dynamicPressSubscale
   
end module
