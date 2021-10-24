subroutine sld_calcVelandAccel(a,dtinv)
    !This routine calculates vel and accel for solid probs
    use typre
    use Mod_TimeIntegrator
    use Mod_php_SetTimeIntegrator
    use Mod_Solids
    implicit none
    class(SolidsProblem)    :: a
    type(TimeIntegratorDt1) :: Integrator1
    type(TimeIntegratorDt2) :: Integrator2
    real(rp)                :: c1(10),c2(10)
    integer(ip)             :: ndime,idime,nsteps,nsteps2,i
    real(rp)                :: dtinv,dtinv2,dtime,beta,omega
    real(rp)                :: LHSdtinv2,LHSdtinv

    dtinv2=dtinv*dtinv
    dtime = 1.0_rp/dtinv

    call a%Mesh%GetNdime(ndime)

    if (a%kfl_tsche_2nd_current == 'NEWMA') then

        beta  = a%beta
        omega = a%omega

        do idime = 1,ndime

            a%accel(idime,:,1) = (dtinv2/beta)*(a%disp(idime,:,1) - a%disp(idime,:,3)) -&
                    (dtinv/beta)*a%veloc(idime,:,2) - ((1.0_rp/(2.0_rp*beta))-1.0_rp)*a%accel(idime,:,2)

            a%veloc(idime,:,1) = a%veloc(idime,:,2) + dtime*((1.0_rp-omega)*a%accel(idime,:,2) +&
                    omega*a%accel(idime,:,1))
        enddo

    else 

        call php_SetTimeIntegrator(a,Integrator1,LHSdtinv,nsteps)
        call php_SetTimeIntegrator(a,Integrator2,LHSdtinv2,nsteps2)
        call Integrator1%GetTimeSchemeCoefficients(c1)
        call Integrator2%GetTimeSchemeCoefficients(c2)

        do idime = 1,ndime
            a%accel(idime,:,1) = dtinv2*c2(1)*a%disp(idime,:,1) 
            a%veloc(idime,:,1) =  dtinv*c1(1)*a%disp(idime,:,1) 
        end do

        do i = 2,nsteps2
            do idime = 1,ndime

                a%accel(idime,:,1) = a%accel(idime,:,1) - dtinv2*c2(i)*a%disp(idime,:,i+1)
            enddo
        enddo

        do i = 2,nsteps
            do idime = 1,ndime

                a%veloc(idime,:,1) = a%veloc(idime,:,1) - dtinv*c1(i)*a%disp(idime,:,i+1)
            enddo
        enddo

    endif

end subroutine 
