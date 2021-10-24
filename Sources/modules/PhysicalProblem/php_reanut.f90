subroutine php_reanut(a)
!Generic Physical Problem read numerical treatment
   use typre
   use Mod_PhysicalProblem
   implicit none

   class(PhysicalProblem) :: a

   !For Listener
   real(rp), pointer     :: param(:) => NULL()
   character(5), pointer :: words(:) => NULL()
   integer(ip), pointer  :: nnpar => NULL(),nnwor => NULL()
   character(150) :: outstr

   call a%Listener%SetLunits(a%lun_pdata,a%lun_outpu)
   call a%Listener%getarrs(words,param,nnpar,nnwor)

   outstr = adjustl(trim(a%exmod))//'_REANUT'

   a%kfl_linea = 1                                    !Linearization (RHS=0, Picard=1, Newton=2)
   a%npica     = 1                                    !Number of Picard iteration (Newton's lin.)
   a%sstol     = 1.0d-5                               !Steady-satate tolerance
   a%kfl_normc = 2                                    !L2 norm for convergence
   a%maxit     = 1                                    !One internal iteration
   a%cotol     = 1.0_rp                               !Internal tolerance
   a%cptol     = 1.0_rp                               !Internal tolerance
   a%safet     = 1.0d10                               !Safety factor
   a%kfl_tsche_1st_datafile = 'UNDEF'  !first time derivative du/dt    scheme
   a%kfl_tsche_2nd_datafile = 'BDF1 ' !second time derivative d2u/dt2 scheme
   a%neule = 0
   a%neule_2nd = 0
   a%kfl_elmat_datafile = 1
   a%kfl_elmat_current  = 1
   a%kfl_twost = 0                                    !# Steps for time derivative (BDF1)  OBSOLETE
   a%kfl_eigensolve = .FALSE.
   a%kfl_error      = .FALSE.
   a%extramatrix    = 'NONE '
   a%kfl_initRelax  = 2                                    !# Steps for time derivative (BDF1)  OBSOLETE

   a%kfl_ProjectionType = 0                           !Default projection scheme is by using the lumped mass matrix
   a%EndLoopQuadrature = 'ForceClosedRule'        !Default quadrature for endsteElmope and EnditeElmope is a closed quadrature

   a%RefinerErrorEstimator = 'ZZ'
   a%RefinerErrorCriteria  = 'Tolerance'
   a%RefinerErrorLimits = (/ 1e-3, 1e-4 /)


   call a%SpecificReanut(0)  !Initializations

   !Reach the section
   call a%Listener%listen(outstr)
   do while(words(1)/='NUMER')
      call a%Listener%listen(outstr)
   end do

   !Begin to read data
   do while(words(1)/='ENDNU')
      call a%Listener%listen(outstr)

      if(words(1)=='STEAD') then
        a%sstol = param(1)
      elseif(words(1)=='NORMO') then
         if(a%Listener%exists('L1   ')) then
            a%kfl_normc = 1
         elseif(a%Listener%exists('L-inf')) then
            a%kfl_normc = 0
         endif
      elseif(words(1)=='MAXIM') then
         a%maxit = int(param(1))
      elseif(words(1)=='CONVE') then
         a%cotol = param(1)
      elseif(words(1)=='CPCON') then
         a%cptol = param(1)
      elseif(words(1)=='LINEA') then
         if(a%Listener%exists('RHS  ')) then
            a%kfl_linea = 0
         else if(a%Listener%exists('PICAR')) then
            a%kfl_linea = 1
         elseif(a%Listener%exists('NEWTO')) then
            a%kfl_linea = 2
            a%npica     = a%Listener%getint('STEPS',0,'#Number of Picard iterations')
         endif
      else if(words(1)=='SAFET') then
        a%safet = param(1)
      else if(words(1)=='T1SCH') then
        a%kfl_tsche_1st_datafile = words(2)
      else if(words(1)=='T1EUL') then
        a%neule = int(param(1))
      else if(words(1)=='T2SCH') then
          a%kfl_tsche_2nd_datafile = words(2)
      else if(words(1)=='T2EUL') then
          a%neule_2nd = int(param(1))
      elseif (words(1)=='ELMAT') then
          a%kfl_elmat_datafile = int(param(1))
      else if(words(1)=='TIMEA') then ! OBSOLETE 2014-06-05
          a%kfl_tiacc = int(param(1))
          a%neule = a%Listener%getint('EULER',0,'#Number of Euler time steps')
          if(a%kfl_timei==0) a%kfl_tiacc = 1
      else if(words(1)=='TIMEI') then ! OBSOLETE 2014-06-05
          if(a%Listener%exists('BDF  ')) a%kfl_twost=2            !temporary meaning only to indicate gear OBSOLETE 2014-06-05
      else if(words(1)=='PROJE') then
          if(a%Listener%exists('LUMPE')) then
              a%kfl_ProjectionType = 0 !Lumped mass matrix for projections
              a%EndLoopQuadrature = 'ForceClosedRule'
          elseif(a%Listener%exists('L2   ')) then
              a%kfl_ProjectionType = 1 !Consistent L2 mass matrix for projections
              a%EndLoopQuadrature = 'DefaultRule'
          endif
      else if(words(1)=='ENDQU') then
          if(a%Listener%exists('DEFAU')) then
              a%EndLoopQuadrature = 'DefaultRule' !Default quadrature rule for endite and endste loops
          elseif(a%Listener%exists('CLOSE')) then
              a%EndLoopQuadrature = 'ForceClosedRule' !Force a closed quadrature rule for endite and endste loops
          endif
      else if(words(1)=='ERROR') then
          a%RefinerErrorEstimator = trim(words(2))
          if (a%Listener%exists('TOLER')) then
              a%RefinerErrorCriteria = 'Tolerance'
          elseif (a%Listener%exists('ELEME')) then
              a%RefinerErrorCriteria = 'Elements'
          endif
          a%RefinerErrorLimits(1) = a%Listener%getrea('MAXER',1e-3_rp,'Error Estimator Maximum')
          a%RefinerErrorLimits(2) = a%Listener%getrea('MINER',1e-4_rp,'Error Estimator Minimum')
          if (a%RefinerErrorEstimator == 'SUBSC') a%kfl_error = .true. 
      else if(words(1)=='DOEPS') then
          if(a%Listener%exists('ON   ')) then
              a%kfl_eigensolve = .TRUE.
              a%extramatrix = 'EIGEN'
          endif
      elseif (words(1)=='DIRIC') then
         if (a%Listener%exists('DELET')) then
            a%kfl_DeleteDirichletColumns = .true.
         elseif (a%Listener%exists('ROWSO') .or. a%Listener%exists('NO   ') .or. a%Listener%exists('OFF  ')) then
            a%kfl_DeleteDirichletColumns = .false.
         endif
      endif

      !Inside the Numerical Treatment block
      call a%SpecificReanut(1)

   end do

   !Obtain the correct value for kfl_twost
   if(a%kfl_twost==2.and.a%kfl_tiacc==2) then
      a%kfl_twost=1
   else if(a%kfl_twost==2.and.a%kfl_tiacc==3) then
      a%kfl_twost=2
   else
      a%kfl_twost=0
   end if

   !conversion from kfl_twost and kfl_tiacc to kfl_tsche_1st_datafile  FUTURE: change a%kfl_twost and a%kfl_tiacc by kfl_twost and kfl_tiacc which are local to this subroutine
   if(a%kfl_timei==1) then !kfl_timei: existence of time derivatives du/dt or d2u/dt2
      if (a%kfl_tsche_1st_datafile=='UNDEF') then
         !write(*,*) 'kfl_twost: ',a%kfl_twost
         !write(*,*) 'kfl_tiacc: ',a%kfl_tiacc
         if     ( a%kfl_twost == 1 .and. a%kfl_tiacc == 2 ) then
            a%kfl_tsche_1st_datafile = 'BDF2 '
         elseif ( a%kfl_twost == 2 .and. a%kfl_tiacc == 3 ) then
            a%kfl_tsche_1st_datafile = 'BDF3 '
         elseif ( a%kfl_twost == 0 .and. a%kfl_tiacc == 2 ) then
            a%kfl_tsche_1st_datafile = 'CNOBS'
         else
            a%kfl_tsche_1st_datafile = 'BDF1 '
         end if
      end if
   else
      a%kfl_tsche_1st_datafile = 'NONE ' !maybe can be 'BDF1 ', I have not thought about that
      a%kfl_tsche_2nd_datafile = 'NONE ' !maybe can be 'BDF1 ', I have not thought about that
   end if

   if (a%neule < 0) then !automatic euler 1st time derivative
      if     (a%kfl_tsche_1st_datafile == 'NONE ') then
         a%neule = 0
      elseif (a%kfl_tsche_1st_datafile == 'BDF1 ') then
         a%neule = 0
      elseif (a%kfl_tsche_1st_datafile == 'BDF2 ') then
         a%neule = 1
      elseif (a%kfl_tsche_1st_datafile == 'BDF3 ') then
         a%neule = 2
      elseif (a%kfl_tsche_1st_datafile == 'CNOBS') then
         a%neule = 0
      elseif (a%kfl_tsche_1st_datafile == 'CN   ') then
         a%neule = 0
      elseif (a%kfl_tsche_1st_datafile == 'TR   ') then
         a%neule = 0
      else
         call runend('php_reanut: not ready for kfl_tsche_1st_datafile')
      endif
      write(*,*) 'WARNING: neule was modified to ',a%neule
   end if

   if (a%neule_2nd < 0) then !automatic euler 2nd time derivative
      if     (a%kfl_tsche_2nd_datafile == 'NONE ') then
         a%neule_2nd = 0
      elseif (a%kfl_tsche_2nd_datafile == 'BDF1 ') then
         a%neule_2nd = 0
      elseif (a%kfl_tsche_2nd_datafile == 'BDF2 ') then
         a%neule_2nd = 1
      elseif (a%kfl_tsche_2nd_datafile == 'NEWMA') then
         a%neule_2nd = 0
      else
         call runend('php_reanut: not ready for kfl_tsche_2nd_datafile')
      endif
      write(*,*) 'WARNING: neule_2nd was modified to ',a%neule_2nd
   end if

   if (a%kfl_elmat_datafile == 0) then
      write(*,*) 'WARNING: kfl_elmat_datafile = 0 does not work for some time-dependent boundary conditions'
   end if
   a%kfl_tsche_1st_current = a%kfl_tsche_1st_datafile
   a%kfl_tsche_2nd_current = a%kfl_tsche_2nd_datafile
   !write(*,*) 'php_reanut: kfl_tsche_1st_datafile     = ',a%kfl_tsche_1st_datafile
   !write(*,*) 'php_reanut: kfl_tsche_2nd_datafile = ',a%kfl_tsche_2nd_datafile
   !write(*,*) 'php_reanut: neule         = ',a%neule
   !write(*,*) 'php_reanut: neule_2nd     = ',a%neule_2nd
   !call runend('php_reanut: provisional stop for testing')

   call a%SpecificReanut(100)     !Final operations

end subroutine
