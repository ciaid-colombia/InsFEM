subroutine nsm_InitStats(a)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   
   !Write info
   a%tamin=1.0e20_rp                  ! Minimum tau
   a%tamax=0.0_rp                     ! Maximum tau
   a%tamea=0.0_rp                     ! Mean tau
   a%remin=1.0e20_rp                  ! Minimum Re
   a%remax=0.0_rp                     ! Maximum Re
   a%remea=0.0_rp                     ! Mean Re
   a%nmean=0
end subroutine

subroutine nsm_InGaussStats(a,acden,acvis,gpvno,chale,timom)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   real(rp) :: acden,acvis,gpvno,chale(2),timom
   
   real(rp) :: reyno
      
   ! Compute cell reynolds number and store values.
   a%nmean = a%nmean+1
   reyno = acden*gpvno*chale(2)/acvis
   a%remin = min(a%remin,reyno)
   a%remax = max(a%remax,reyno)
   a%remea = a%remea + reyno

   a%tamin = min(a%tamin,timom)
   a%tamax = max(a%tamax,timom)
   a%tamea = a%tamea + timom
end subroutine

subroutine nsm_FinalizeStats(a)
   use typre
   use Mod_NavierStokes
   implicit none
   class(NavierStokesProblem) :: a
   
   a%tamea=a%tamea/real(a%nmean,rp)
   a%remea=a%remea/real(a%nmean,rp)
   
end subroutine