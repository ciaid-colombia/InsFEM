subroutine rom_postpr(a)
   use typre
   use Mod_Postpr
   use Mod_PodRom
   implicit none
   class(PodRomProblem) :: a
   character(200)       :: nameP
   integer(ip)          :: idofn,idofr

   do idofr = 1,a%ndofr
      do idofn = 1,a%ndofn
         write(nameP,"(A6,A150)") 'Basis ',a%nameBasis(idofn)
         call a%Problem%FilePostpr%postpr(a%Basis(idofn,:,idofr),trim(adjustl(nameP)),idofr-1,real(idofr-1,rp),a%Mesh)
      end do
   end do

!   do idofn = 1,a%ndofn
!      write(nameP,"(A5,A150)") 'Mean ',a%nameBasis(idofn)
!      call a%Problem%FilePostpr%postpr(a%SnapMean(idofn,:),trim(adjustl(nameP)),1,0.0_rp,a%Mesh)
!   end do
   

end subroutine rom_postpr
