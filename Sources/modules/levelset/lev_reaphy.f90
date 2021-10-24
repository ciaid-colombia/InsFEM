subroutine lev_reaphy(a,itask)
   use typre
   use Mod_Listen
   use Mod_Memor
   use Mod_MPIObject
   use Mod_LevelSet
   implicit none
   
   integer(ip) :: itask
   class(LevelSetProblem) :: a
	
	integer(ip) :: imate,istat,npara,ifunp,nfunp
	!For a%Listener%listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
 
   call a%Listener%getarrs(words,param,nnpar,nnwor)
      
   if (itask == 0) then
   
   !Initializations
   a%kfl_advec = 1                                    ! Default is convection on
   a%kfl_ExactLevel = 0                               ! Exact LevelSet initialization
   !Problem data
   elseif (itask == 1) then 
         if(words(1)=='TEMPO') then                    ! Temporal evolution
            if(a%Listener%exists('ON   ')) a%kfl_timei = 1 
         elseif(words(1)=='EXACT') then
            if(words(2)== 'ON   ')then
               a%kfl_ExactLevel =1
            end if
         elseif(words(1)=='ADVEC') then                    !Advection
            if(a%Listener%exists('ON   '))  then
               a%kfl_advec = 1    
            elseif(a%Listener%exists('OFF  ')) then
               a%kfl_advec = 0    
            endif
            
         end if         
   
   !Properties
   elseif(itask == 2) then   
 
   
   !Other initializations
   elseif (itask == 100) then
   
   
   endif

end subroutine lev_reaphy

