subroutine  sup_reabcs(a,itask,kflag)

   use typre
   use Mod_Listen   
   use Mod_Postpr
   use Mod_Mesh
   use Mod_PhysicalProblem
   use Mod_NavierStokes
   use Mod_ThreeField
   implicit none
  
   class(ThreeFieldNSProblem) :: a
   integer(ip) :: itask
   integer(ip), optional :: kflag   
   integer(ip) :: ndime,ntens  
   
   interface
      subroutine nsi_Reabcs(a,itask,kflag)
         use typre
         use Mod_NavierStokes
         implicit none
         class(NavierStokesProblem) :: a
         integer(ip) :: itask
         integer(ip), optional :: kflag
      
      end subroutine
   end interface    
       
   if (itask == 0) then
      
      call nsi_Reabcs(a,0)
      
   elseif(itask == 1) then     

      call nsi_Reabcs(a,1)  
    
   elseif(itask == 100) then
   
      call nsi_Reabcs(a,100)
      
      !Dimensions
      call a%Mesh%GetNdime(ndime)  
      
      ntens=(ndime-1)*(ndime-1) + 2
      
      if(a%ndofbc == (ndime + ntens))then
         call a%FilePostpr%postpr(a%bvess(1:ntens,:,1),adjustl(trim(a%exmod))//'_SBC',0_ip,0.0_rp,a%Mesh)
         call a%FilePostpr%postpr(a%kfl_fixno(1:ntens,:),adjustl(trim(a%exmod))//'_SFIXNO',0_ip,0.0_rp,a%Mesh)
         call a%FilePostpr%postpr(a%bvess(ntens + 1:a%ndofbc,:,1),adjustl(trim(a%exmod))//'_UBC',0_ip,0.0_rp,a%Mesh)
         call a%FilePostpr%postpr(a%kfl_fixno(ntens+1:a%ndofbc,:),adjustl(trim(a%exmod))//'_UFIXNO',0_ip,0.0_rp,a%Mesh)   
      endif 
   
   endif  
      
end subroutine


subroutine sup_ReadOnNodes(a)
   use typre
   use Mod_Mesh
   use Mod_ThreeField
   use Mod_SupExacso  
   implicit none
   class(ThreeFieldNSProblem) :: a
   integer(ip)       :: ntens,ndime
   real(rp)          :: uex,vex
   integer(ip)       :: aux1,aux2
   logical           :: isALE
   
   !Dimensions
   call a%Mesh%GetALE(isALE)   
   call a%Mesh%GetNdime(ndime)  
      
   if (isALE) a%bvess(:,a%gipoin,2)=a%bvess(:,a%gipoin,1)

end subroutine

