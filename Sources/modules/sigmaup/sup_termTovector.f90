   subroutine sup_termTovector(ndime,auxtens,gpvel,grvel,grsig,GrExpPsi,ExpPosMatrix,grpre,elext2,gpconv,gpdivs,gpgrap,auxLCR)
      use typre

      implicit none   
      integer(ip), intent(in)   :: ndime,auxtens,auxLCR
      real(rp), intent(in)      :: gpvel(ndime),grsig(auxtens,ndime),grvel(ndime,ndime),grpre(ndime),elext2(ndime)
      real(rp), intent(in)      :: ExpPosMatrix(ndime,ndime) ,GrExpPsi(auxtens,ndime)
      real(rp), intent(out)     :: gpconv(ndime),gpdivs(ndime),gpgrap(ndime)
      integer(ip)               :: idime, i, j, k
      real(rp)                  :: grsigmat(ndime,ndime,ndime),aux_gpdivs(ndime)
      
      if(ndime==2)then

         gpconv(1)=gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2)
         gpconv(2)=gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2)
         
         if (auxLCR==0) then
            gpdivs(1)=grsig(1,1) + grsig(3,2)                          
            gpdivs(2)=grsig(3,1) + grsig(2,2)     
         else if (auxLCR==1) then
            aux_gpdivs(1)=grsig(1,1) + grsig(3,2)                          
            aux_gpdivs(2)=grsig(3,1) + grsig(2,2) 
            gpdivs=0.0_rp
            do k=1,ndime
               do j=1,ndime
                  gpdivs(k)=ExpPosMatrix(k,j)*aux_gpdivs(j) + gpdivs(k)
               end do
            end do

!             gpdivs(1)=GrExpPsi(1,1) + GrExpPsi(3,2)                          
!             gpdivs(2)=GrExpPsi(3,1) + GrExpPsi(2,2)
!             
         end if
         
         gpgrap(1)=grpre(1) - elext2(1)                         
         gpgrap(2)=grpre(2) - elext2(2)                                 
         
         
       elseif(ndime==3)then
       
         gpconv(1)=gpvel(1)*grvel(1,1) + gpvel(2)*grvel(1,2) + gpvel(3)*grvel(1,3)
         gpconv(2)=gpvel(1)*grvel(2,1) + gpvel(2)*grvel(2,2) + gpvel(3)*grvel(2,3) 
         gpconv(3)=gpvel(1)*grvel(3,1) + gpvel(2)*grvel(3,2) + gpvel(3)*grvel(3,3)

         if (auxLCR==0) then
            gpdivs(1)=grsig(1,1) + grsig(6,2) + grsig(5,3)                          
            gpdivs(2)=grsig(6,1) + grsig(2,2) + grsig(4,3)                         
            gpdivs(3)=grsig(5,1) + grsig(4,2) + grsig(3,3)
         else if (auxLCR==1) then
            gpdivs(1)=GrExpPsi(1,1) + GrExpPsi(6,2) + GrExpPsi(5,3)                          
            gpdivs(2)=GrExpPsi(6,1) + GrExpPsi(2,2) + GrExpPsi(4,3)                         
            gpdivs(3)=GrExpPsi(5,1) + GrExpPsi(4,2) + GrExpPsi(3,3)
         end if   

         gpgrap(1)=grpre(1) - elext2(1)                         
         gpgrap(2)=grpre(2) - elext2(2)   
         gpgrap(3)=grpre(3) - elext2(3)
         
      end if
          

   end subroutine sup_termTovector


   subroutine sup_termTovector2(ndime,auxtens,gpvel,grvel,divu,gradu_vector)
      use typre
      use Mod_supm_StressGenerator
      implicit none   
      integer(ip), intent(in)   :: ndime,auxtens
      real(rp), intent(in)      :: gpvel(ndime),grvel(ndime,ndime)
      real(rp), intent(out)     :: divu(1),gradu_vector(auxtens)
      real(rp)                  :: gradu(ndime,ndime)
      integer(ip)               ::  i, j

      divu=0.0_rp
      gradu=0.0_rp
      do i=1,ndime
         divu=grvel(i,i)+divu
         
         do j=1,ndime
            gradu(i,j)=0.5_rp*(grvel(i,j)+grvel(j,i))
         end do   

         call PassSymMatrixToVector(ndime,auxtens,gradu,gradu_vector)
         
      end do
   end subroutine sup_termTovector2
   
   
   subroutine sup_termTovector_viscoelastic(ndime,auxtens,gpvel,grvel,gpsig,grsig,gpconvsigma,gpdeform)
      use Mod_supm_StressGenerator
      use typre
      implicit none   
      integer(ip), intent(in)   :: ndime,auxtens
      real(rp), intent(in)      :: gpvel(ndime),grsig(auxtens,ndime),grvel(ndime,ndime),gpsig(auxtens)
      real(rp), intent(out)     :: gpconvsigma(auxtens), gpdeform(auxtens)
      integer(ip)               :: idime, i, j, k
      real(rp)                  :: grsigMat(ndime,ndime,ndime), gpsigMat(ndime,ndime)
      real(rp)                  :: aux_gpconvsigma(ndime,ndime), aux_gpdef(ndime,ndime)

      grsigMat=0.0_rp
      gpsigMat=0.0_rp
      
      call PassGradStressToTensor3B(ndime,auxtens,grsig,grsigMat)
      call sup_SigmaMatrix(ndime,auxtens,gpsig,gpsigMat)
      
      aux_gpconvsigma=0.0_rp
      aux_gpdef=0.0_rp
      
      do i=1,ndime
         do j=1,ndime
         
            do k=1,ndime
               !Convective Term
               aux_gpconvsigma(i,j)=gpvel(k)*grsigMat(i,j,k) + aux_gpconvsigma(i,j)
               
               !Deformation Terms
               aux_gpdef(i,j)=gpsigMat(i,k)*grvel(j,k) +grvel(i,k)*gpsigMat(k,j) + aux_gpdef(i,j)
            end do
            
         end do   
      end do
      
      call PassSymMatrixToVector(ndime,auxtens,aux_gpconvsigma,gpconvsigma)
      call PassSymMatrixToVector(ndime,auxtens,aux_gpdef,gpdeform)
   end subroutine sup_termTovector_viscoelastic


   subroutine sup_termTovector_log(ndime,auxtens,gpvel,grvel,GrExpPsi,ExpPosMatrix,gpconvsigma,gpdeform)
      use Mod_supm_StressGenerator
      use typre
      implicit none   
      integer(ip), intent(in)   :: ndime,auxtens
      real(rp), intent(in)      :: gpvel(ndime),grvel(ndime,ndime)
      real(rp), intent(in)      :: ExpPosMatrix(ndime,ndime) ,GrExpPsi(auxtens,ndime)
      real(rp), intent(out)     :: gpconvsigma(auxtens), gpdeform(auxtens)
      integer(ip)               :: idime, i, j, k
      real(rp)                  :: GrExpPsimat(ndime,ndime,ndime)
      real(rp)                  :: aux_gpconvsigma(ndime,ndime), aux_gpdef(ndime,ndime)

      GrExpPsimat=0.0_rp
      call PassGradStressToTensor3B(ndime,auxtens,GrExpPsi,GrExpPsimat)
      
      aux_gpconvsigma=0.0_rp
      aux_gpdef=0.0_rp
      
      do i=1,ndime
         do j=1,ndime
            do k=1,ndime
               !u*grad(exp)
               aux_gpconvsigma(i,j)=gpvel(k)*GrExpPsimat(i,j,k) + aux_gpconvsigma(i,j)
               
               !Exp(S)*gradu
               aux_gpdef(i,j)=ExpPosMatrix(i,k)*grvel(j,k) + grvel(i,k)*ExpPosMatrix(k,j) + aux_gpdef(i,j)
               
            end do
            
            aux_gpdef(i,j)= -(grvel(i,j)+grvel(j,i)) + aux_gpdef(i,j)
            
         end do   
      end do
      
      call PassSymMatrixToVector(ndime,auxtens,aux_gpconvsigma,gpconvsigma)
      call PassSymMatrixToVector(ndime,auxtens,aux_gpdef,gpdeform)
      
end subroutine sup_termTovector_log