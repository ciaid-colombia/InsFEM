module Mod_nsm_NewtonRaphson
   use typre
   use Mod_PointerSetter
   use Mod_nsm_BaseElmope
   implicit none
   private
   public SetPointersNewtonRaphson
   
   type, extends(PointerSetter) :: SPNewtonRaphson
contains
      procedure :: SpecificSet => SpecificSetNewtonRaphson
   end type
   type(SPNewtonRaphson) :: SetPointersNewtonRaphson
   
contains

   subroutine SpecificSetNewtonRaphson(d)
      implicit none
      class(SPNewtonRaphson) :: d
         
      call ConcatenateProcedures(ProcHook_InGaussElmats,AddNRterms)
   end subroutine
   
   !-----------------------------------------------------------------------
   subroutine AddNRterms
      implicit none
      integer(ip) :: inode,jnode,idime,jdime
      
      if (a%istep>2) then
         do jnode=1,e%pnode
            do inode=1,e%pnode
               do idime=1,e%ndime
                  do jdime=1,e%ndime
                     elmat(idime,inode,idime,jnode) = elmat(idime,inode,idime,jnode) + a%MatProp(imat)%densi*e%shape(inode,e%igaus)*e%shape(jnode,e%igaus)*grvel(jdime,idime)*dvol
                  enddo
               enddo
            enddo
         enddo

         do inode=1,e%pnode
            do idime=1,e%ndime
               do jdime=1,e%ndime
                  elrhs(idime,inode) = elrhs(idime,inode) + a%MatProp(imat)%densi*e%shape(inode,e%igaus)*gpadv(idime)*grvel(jdime,idime)*dvol
               enddo
            enddo
         enddo
      endif
   end subroutine

end module
