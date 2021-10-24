subroutine trace(ndime,tensor,res)
 use typre
 implicit none
 integer(ip), intent(in) :: ndime
 real(rp)   , intent(in) :: tensor(ndime,ndime)
 real(rp)   , intent(out):: res
 integer(ip)             :: idime

 res = 0.0_rp

 do idime = 1,ndime 
     res = res + tensor(idime,idime)
 enddo

end subroutine trace

subroutine get2ndIITensor_voigt(ndime,sz,ii)
 use typre
 implicit none
 integer(ip), intent(in) :: ndime,sz
 real(rp), intent(inout) :: ii(sz)

 if (ndime.eq.2) then

     ii =  reshape([ 1.0_rp, 1.0_rp, 0.0_rp], &
         [sz])

 elseif (ndime.eq.3) then

     ii =  reshape([ 1.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
         [sz])

 endif

end subroutine get2ndIITensor_voigt

subroutine get2ndIITensor(ndime,ii)
 use typre
 implicit none
 integer(ip), intent(in)    :: ndime
 real(rp)   , intent(inout) :: ii(ndime,ndime)

 ii = 0.0_rp
 ii(1,1)=1.0_rp
 ii(2,2)=1.0_rp

 if (ndime.eq.3) then
     ii(3,3)=1.0_rp
 end if

end subroutine get2ndIITensor

subroutine getDiagTensor(ndime,sz,ii)
 use typre
 implicit none
 integer(ip), intent(in)    :: ndime,sz
 real(rp)   , intent(inout) :: ii(sz,sz)

 if (ndime.eq.2) then

      ii =  reshape([ 1.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 1.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 1.0_rp], &
             [3,3])

 elseif (ndime.eq.3) then

      ii =  reshape([ 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], &
             [6,6])

     endif
end subroutine getDiagTensor

subroutine getDiagTensor2Txy(ndime,sz,ii)
 use typre
 implicit none
 integer(ip), intent(in)    :: ndime,sz
 real(rp)   , intent(inout) :: ii(sz,sz)

 if (ndime.eq.2) then

      ii =  reshape([ 1.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 1.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 2.0_rp], &
             [3,3])

 elseif (ndime.eq.3) then

      ii =  reshape([ 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 2.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 2.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 2.0_rp], &
             [6,6])

     endif

end subroutine getDiagTensor2Txy

subroutine get4thIITensor(ndime,sz,ii)
 use typre
 implicit none
 integer(ip), intent(in)    :: ndime,sz
 real(rp)   , intent(inout) :: ii(sz,sz)

 if (ndime.eq.2) then

      ii =  reshape([ 1.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 1.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 1.0_rp], &
             [3,3])

 elseif (ndime.eq.3) then

      ii =  reshape([ 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], &
             [6,6])

     endif

end subroutine get4thIITensor

subroutine get4thIIVolumetricTensor(ndime,sz,ii)
 use typre
 implicit none
 integer(ip), intent(in)    :: ndime,sz
 real(rp)   , intent(inout) :: ii(sz,sz)

 if (ndime.eq.2) then

      ii =  reshape([ 1.0_rp, 1.0_rp, 0.0_rp, &
                      1.0_rp, 1.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp], &
             [3,3])

 elseif (ndime.eq.3) then

      ii =  reshape([ 1.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      1.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      1.0_rp, 1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                      0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
             [6,6])

     endif
     
     ii = (1.0_rp/3.0_rp)*ii

end subroutine get4thIIVolumetricTensor

subroutine pushForward(ndime,sz,F,tensor,res)
 use typre
 implicit none
 integer(ip), intent(in)  :: ndime,sz
 real(rp)   , intent(in)  :: F(ndime,ndime),tensor(ndime,ndime)
 real(rp)   , intent(out) :: res(ndime,ndime)
 real(rp)             :: auxmat(ndime,ndime)

 auxmat = 0.0_rp
 res = 0.0_rp

 auxmat = matmul(tensor,transpose(F))
 res    = matmul(F,auxmat)

end subroutine pushForward

subroutine getVoigtStrain(sz,ndime,t,voigt)
 use typre
 implicit none
 integer(ip), intent(in)  :: ndime,sz
 real(rp)   , intent(in)  :: t(ndime,ndime)
 real(rp)   , intent(out) :: voigt(sz,1)

 if (ndime.eq.2) then

     voigt =  reshape([ t(1,1), t(2,2), 2.0_rp*t(1,2)], &
         [sz,1])

 elseif (ndime.eq.3) then

     voigt =  reshape([ t(1,1), t(2,2), t(3,3), 2.0_rp*t(2,3), 2.0_rp*t(1,3), 2.0_rp*t(1,2)], &
         [sz,1])

 endif

end subroutine getVoigtStrain 

subroutine getVoigtStress(sz,ndime,t,voigt)
 use typre
 implicit none
 integer(ip), intent(in)  :: ndime,sz
 real(rp)   , intent(in)  :: t(ndime,ndime)
 real(rp)   , intent(out) :: voigt(sz)

 if (ndime.eq.2) then

     voigt =  reshape([ t(1,1), t(2,2), t(1,2)], &
         [sz])

 elseif (ndime.eq.3) then

     voigt =  reshape([ t(1,1), t(2,2), t(3,3), t(2,3), t(1,3), t(1,2)], &
         [sz])

 endif

end subroutine getVoigtStress

subroutine getStrainTensor(sz,ndime,voigt,t)
 use typre
 implicit none
 integer(ip), intent(in)  :: ndime,sz
 real(rp)   , intent(in)  :: voigt(sz)
 real(rp)   , intent(out) :: t(ndime,ndime)

 t=0.0_rp

 if (ndime.eq.2) then

     t =  reshape([      voigt(1),0.5_rp*voigt(3),&
                  0.5_rp*voigt(3),voigt(2)],&
         [ndime,ndime])

 elseif (ndime.eq.3) then

     t =  reshape([ voigt(1)     ,0.5_rp*voigt(6),0.5_rp*voigt(5),&
                  0.5_rp*voigt(6),       voigt(2),0.5_rp*voigt(4),&
                  0.5_rp*voigt(5),0.5_rp*voigt(4),       voigt(3)],&
         [ndime,ndime])

 endif

end subroutine getStrainTensor

subroutine getStressTensor(sz,ndime,voigt,t)
 use typre
 implicit none
 integer(ip), intent(in)  :: ndime,sz
 real(rp)   , intent(in)  :: voigt(sz)
 real(rp)   , intent(out) :: t(ndime,ndime)

 t=0.0_rp

 if (ndime.eq.2) then

     t =  reshape([voigt(1),voigt(3),&
                   voigt(3),voigt(2)],&
         [ndime,ndime])

 elseif (ndime.eq.3) then

     t =  reshape([ voigt(1),voigt(6),voigt(5),&
                    voigt(6),voigt(2),voigt(4),&
                    voigt(5),voigt(4),voigt(3)],&
         [ndime,ndime])

 endif

end subroutine  

subroutine matvec(s1,s2,lmat,rvec,res)
 use typre
 implicit none
 integer(ip), intent(in)  :: s1,s2
 real(rp)   , intent(in)  :: lmat(s1,s2)
 real(rp)   , intent(in)  :: rvec(s2)
 real(rp)   , intent(out) :: res(s1)
 integer(ip) :: j,k
 
 res=0.0
 do j = 1,s1
   do k = 1,s2
      res(j) = res(j) + lmat(j,k)*rvec(k)
   enddo
enddo

 
end subroutine

subroutine matvecmult(s1,s2,lmat,rvec,res)
 use typre
 implicit none
 integer(ip), intent(in)  :: s1,s2
 real(rp)   , intent(in)  :: lmat(s1*s2,s1*s2)
 real(rp)   , intent(in)  :: rvec(s1*s2)
 real(rp)   , intent(out) :: res(s1*s2)
 integer(ip) :: j,k
 !real(rp),    intent(in)    :: lmat(s1,s2,s1,s2)
 !real(rp),    intent(in)    :: rvec(s1,s2)
 !real(rp),    intent(out)   :: res(s1,s2)
 !integer(ip) :: j,k,l,m
 
 res=0.0

 !Both forms produce the same result, second one
 !makes use of the way fortran saves matrices into
 !memory, has two less for loops, therefore efficient
 
 !do j=1,s2
 !    do k=1,s2
 !        do l=1,s1
 !           do m=1,s1
 !               res(l,m) = res(l,m) + lmat(m,j,l,k)*rvec(l,k)
 !           end do
 !        end do
 !    end do
 !end do

 do j = 1,s1*s2
   do k = 1,s1*s2
      res(j) = res(j) + lmat(j,k)*rvec(k)
   enddo
enddo

 
end subroutine

subroutine matmatmulti(s1,s2,lmat,rmat,res)
 use typre
 implicit none
 integer(ip), intent(in)  :: s1,s2
 real(rp)   , intent(in)  :: lmat(s1*s2,s1*s2)
 real(rp)   , intent(in)  :: rmat(s1*s2,s1*s2)
 real(rp)   , intent(out) :: res(s1*s2,s1*s2)
 integer(ip) :: i,j,k
 
 res=0.0

 do i = 1,s1*s2
     do j = 1,s1*s2
         do k = 1,s1*s2
             res(i,j) = res(i,j) + lmat(i,k)*rmat(k,j)
         enddo
     enddo
enddo

 
end subroutine

subroutine tracedoubledot(ndime,tensor,res)
 !This routine calculates the following
 !tr(A**2)
 use typre
 implicit none
 integer(ip), intent(in)  :: ndime
 real(rp)   , intent(in)  :: tensor(ndime,ndime)
 real(rp)   , intent(out) :: res
 integer(ip)              :: i,j

 res = 0.0_rp

      if (ndime.eq.2) then
          res = res + tensor(1,1)*tensor(1,1)&
              & + tensor(2,2)*tensor(2,2) + 2.0_rp*tensor(1,2)*tensor(1,2)
              
      else
          res = res + tensor(1,1)*tensor(1,1)&
              & + tensor(2,2)*tensor(2,2)+ tensor(3,3)*tensor(3,3) &
              & + 2.0_rp*tensor(1,2)*tensor(1,2)+2.0_rp*tensor(1,3)*tensor(1,3)& 
              & + 2.0_rp*tensor(2,3)*tensor(2,3)
      endif

end subroutine

subroutine getTensorInvariants(ndime,tensor,invar)
 !This routine calculates the invariants
 !of a tensor
 use typre
 implicit none
 integer(ip), intent(in)    :: ndime
 real(rp)   , intent(in)    :: tensor(ndime,ndime)
 real(rp)   , intent(inout) :: invar(ndime)
 real(rp)  :: auxvar

 !Invariant 1
 auxvar = 0.0_rp
 call trace(ndime,tensor,auxvar)
 invar(1)= auxvar

 !Invariant 2
 auxvar = 0.0_rp
 call tracedoubledot(ndime,tensor,auxvar)
 invar(2) = 0.5*(invar(1)*invar(1)-auxvar)

 if (ndime.eq.3) then
  !Invariant 3
  auxvar = 0.0_rp
  call deter(tensor,auxvar,ndime)
  invar(3)= auxvar
 endif

end subroutine
