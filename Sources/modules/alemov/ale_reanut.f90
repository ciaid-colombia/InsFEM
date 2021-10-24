subroutine ale_reanut(a,itask)
   use typre
   use Mod_Alemov
   implicit none
   
   class (AlemovProblem) :: a
   integer(ip)           :: itask
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   
   call a%Listener%getarrs(words,param,nnpar,nnwor)
   
   !Initializations
   if (itask == 0) then
      a%kfl_timei = 1
      a%kfl_isFixedMesh = 0
      a%kfl_RemeshingCriteria = 0
      
   elseif (itask == 1) then
      if(words(1)=='FIXED') then
         if(a%Listener%exists('EXTER')) a%kfl_isFixedMesh = 1
      
      elseif (words(1) == 'REMES') then
         if (a%Listener%exists('FOLDI')) then
            a%kfl_RemeshingCriteria = 0
         
         elseif (a%Listener%exists('ALWAY')) then
            a%kfl_RemeshingCriteria = 2
         
         endif

      endif
      
   elseif (itask == 100) then !Final operations

   endif

end subroutine   
