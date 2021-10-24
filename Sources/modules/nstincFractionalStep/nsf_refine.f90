subroutine nsf_refine(a,itask)
   use typre
   use Mod_NSFractionalStep
   implicit none
   class(NSFractionalStepProblem) :: a
   character(6) :: itask

   integer(ip) :: oldnpoin,oldnboun,newnpoin,newnboun,icomp,ndime
   real(rp), allocatable :: auxbvpress(:)
   integer(ip), allocatable :: auxkfl_fixpr(:)
   
   interface
      subroutine nsi_refine(a,itask) 
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         character(6) :: itask
      end subroutine
   end interface
   
   !This subroutine transforms the tempe arrays from one mesh to the other 
   !Adaptive Mesh Refinement
   !Old Dimensions
   oldnpoin = size(a%veloc,2)
   oldnboun = size(a%kfl_bours)
   
   !First we do as nsi
   call nsi_refine(a,itask)
   
   !Then additional things
   
   !We need to modify all the arrays
   !We assume that the mesh has already been updated
   call a%Mesh%GetNpoin(newnpoin)
   call a%Mesh%GetNboun(newnboun)
   call a%Mesh%GetNdime(ndime)
   
   !kfl_fixpr
   call a%Memor%alloc(newnpoin,auxkfl_fixpr,'kfl_fixpr','nsf_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%kfl_fixpr,auxkfl_fixpr,'minim')
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%kfl_fixpr,auxkfl_fixpr)
   endif   
   call move_alloc(auxkfl_fixpr,a%kfl_fixpr)
   call a%Memor%deallocObj(0,'kfl_fixpr','nsf_Refine',ip*oldnpoin)
   
   !bvpress
   call a%Memor%alloc(newnpoin,auxbvpress,'bvpress','nsf_Refine')
   if (itask == 'Refine') then
      call a%Refiner%UpdateVariable(1_ip,a%bvpress,auxbvpress)
   elseif (itask == 'Rebala') then
      call a%Refiner%RebalanceVariable(1_ip,a%bvpress,auxbvpress)
   endif   
   call move_alloc(auxbvpress,a%bvpress)
   call a%Memor%deallocObj(0,'bvpress','nsf_Refine',rp*oldnpoin)
   

   
end subroutine 