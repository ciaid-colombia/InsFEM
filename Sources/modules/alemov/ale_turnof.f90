subroutine ale_turnof(a)
   use typre
   use def_parame
   use Mod_iofile
   use Mod_Alemov
   use Mod_sldAlemov
   implicit none
   class(AlemovProblem) :: a

   integer(ip) :: npoin,ncomp,oldnpoin

   call a%Mesh%GetNpoin(npoin)

   ncomp = a%ncomp

   call a%Memor%dealloc(a%ndofbc,npoin,a%ncomp,a%Displacement,'Displacement','ale_turnof')
   call a%Memor%dealloc(a%ndofbc,npoin,a%ncomp,a%Velocity,'Velocity','ale_turnof')
   call a%Memor%dealloc(a%ndofbc,npoin,a%kfl_fixno0,'kfl_fixno0','ale_turnof')

   if (a%kfl_rstar == 1) then
      if (a%kfl_inter == 1) then
         call a%OldMesh%GetNpoin(oldnpoin)
      else
         oldnpoin = npoin
      endif

      call a%Memor%dealloc(a%ndofbc,oldnpoin,a%oldncomp,a%restar_displ,'restar_displ','ale_restar')
      call a%Memor%dealloc(a%ndofbc,oldnpoin,a%ncomp   ,a%restar_veloc,'restar_veloc','ale_restar')
   endif

   select type (a)
   type is (sldAlemovProblem)
       call a%Memor%dealloc(a%ndofbc,npoin,2,a%bdisp,'bdisp','ale_turnof')
       call a%Memor%dealloc(a%ndofbc,npoin,3,a%bres,'bres','ale_turnof')
   end select

   if (a%nptra/=0) then
       if (a%MPIrank == a%MPIroot) then
           call iofile(two,a%lun_trap,a%fil_trap,'ALE TRACKING OF POINTS')
       endif
   endif

   call a%deferredTurnof
end subroutine
