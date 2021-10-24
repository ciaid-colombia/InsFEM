subroutine sld_bouope(a,ndime,iboun,e,sz,bc,elmat,elrhs,ndofn,npoinLocal)
!This subroutine computes the boundary contribution for the Cauchy eqn
!-----------------------------------------------------------------------
   use typre
   use Mod_Element
   use Mod_Solids
   implicit none
   class(SolidsProblem),target :: a
   class(FiniteElement)        :: e
   integer(ip),intent(in)      :: npoinLocal,iboun,ndime,ndofn,bc
   real(rp),intent(inout)      :: elmat(ndofn,e%mnode,ndofn,e%mnode)              ! Element matrices
   real(rp),intent(inout)      :: elrhs(ndofn,e%mnode)
   integer(ip) :: igaub,inodb,inode,ilnod,idime,sz,ielem
   real(rp)    :: tract(ndime,e%pnodb),tract_b(ndime)
   real(rp)    :: dsurf,dsurf0
   real(rp)    :: c_el(sz,sz),w
      
      call e%elmdel
      
      !Gather operations
      tract=0.0_rp
      tract_b=0.0_rp

      if(associated(a%trac)) then

          call e%gatherb(ndime,tract,a%trac(:,:))

      else if(a%kfl_fixbo(iboun) == 2) then

          tract_b(1:e%ndime) = a%bvnat(iboun)%a(1:e%ndime)

      end if

      dsurf0 = 0.0_rp

      !Gauss-Point Loop
      do igaub=1,e%pgaub
         e%igaub = igaub

         ielem = e%lboel(e%pnodb+1) 
         a%btraction(ielem)%a(:,igaub) = 0.0_rp
   
         !Calculate exterior Normal
         call e%bounor

         !Derivatives at the boundary
         call e%elmderb
         
         dsurf=e%weigb(e%igaub)*e%eucta
         dsurf0 = dsurf0 + dsurf
         

         do inodb=1,e%pnodb
             inode = e%lboel(inodb)
             ilnod = e%lnods(inode)

             do idime=1,ndime
                 !elrhs= elrhs + N_i*t_f*ds
                 elrhs(bc+idime,inode)=elrhs(bc+idime,inode)    &
                     & - e%shapb(inodb,e%igaub)*tract(idime,inodb)*dsurf  & !FSI (or ext distributed) nodal tractions
                     & + e%shapb(inodb,e%igaub)*tract_b(idime)*dsurf        !Prescribed boundary tractions

             end do   

             !We store the above for output
             a%btraction(ielem)%a(:,igaub) = a%btraction(ielem)%a(:,igaub)  &
                 & - e%shapb(inodb,e%igaub)*tract(:,inodb) 

         end do

         !We store the above for output
         a%btraction(ielem)%a(:,igaub) = a%btraction(ielem)%a(:,igaub) + tract_b(:) 
     end do

 end subroutine sld_bouope
