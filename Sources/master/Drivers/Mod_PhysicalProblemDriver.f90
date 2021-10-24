module Mod_PhysicalProblemDriver
    use typre
    use Mod_caseVariables
    use Mod_BroadCastBuffer
    use Mod_DriverInterface
    use Mod_PhysicalProblem
    use def_parame
    implicit none

contains

      subroutine physical_Turnon(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          !Physical Problem data
          call Physical%SetMPI(c%masterVars%MPIcomm,c%masterVars%MPIsize,c%masterVars%MPIroot,c%masterVars%MPIrank)
          call Physical%SetReadType(c%masterVars%ReadTypeString)
          call Physical%SetInputFolder(c%masterVars%BaseDataFolder)
          call Physical%SetRestartFolder(c%masterVars%RestartFolder)
          call Physical%SetOldRestartFile(c%masterVars%RestartFolder)
          call Physical%SetInputFile(c%masterVars%namda)
          call Physical%SetOutputFolder(c%masterVars%ResultsFolder)
          call Physical%SetOutputFiles(c%masterVars%lun_memor,c%masterVars%lun_outpu)
          call Physical%SetMesh(c%domainVars%Mesh)
          call Physical%SetPostprFile(c%masterVars%FilePostpr)
          call Physical%SetIOPointers(c%masterVars%Writerpr,c%masterVars%Readerpr)
          call Physical%SetFlush(c%masterVars%kfl_flush)
          call Physical%SetIOFormat(c%masterVars%kfl_iofor)
          call Physical%SetMPICommunicationsType(c%masterVars%kfl_MPIComType)
          call Physical%SetParallelLibrary(c%masterVars%ParallelLibrary)
          call Physical%SetCoupledConvFlag(c%masterVars%kfl_doCoupledConv)
          call Physical%SetAdaptiveFlag(c%adaptiveVars%kfl_AdaptiveRefinement)
          call Physical%SetMulticommData(c%masterVars%kfl_multicomm,c%masterVars%multicommColor)
          call Physical%Turnon

          !Adaptive Mesh Refinement
          if (c%adaptiveVars%kfl_AdaptiveRefinement == 1 .or. c%adaptiveVars%NumberOfInitialUniformRefinementSteps > 0) then
              call Physical%SetAdaptiveRefiner(c%adaptiveVars%Refiner)
          endif

          !Restart with interpolation
          if (a%kfl_rstar==1) then
            call Physical%SetRestart
            if (a%kfl_inter==1) then
               call Physical%SetInterpolation
               call Physical%SetOldInputFile(c%masterVars%oldnamda)
               call Physical%SetOldRestartFile(c%masterVars%OldRestartFolder)
               call c%domainVars%OldMesh%SetALE(0_ip)
               call Physical%SetOldMesh(c%domainVars%OldMesh)
               call Physical%SetInterp(c%masterVars%Int_Restart)
               call Physical%Restart(one)
            else
               call Physical%SetOldInputFile(c%masterVars%namda)
               call Physical%Restart(one)
            endif
          endif

      end subroutine

      subroutine physical_GetTimeStep(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          real(rp) :: aux_dtinv = 0.0_rp

          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%GetSte(aux_dtinv)

          aux_dtinv = max(c%masterVars%dtinv,aux_dtinv)
          c%masterVars%dtinv = aux_dtinv

      end subroutine

      subroutine physical_SetTimeStep(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical


          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%Setste(c%masterVars%dtinv,c%masterVars%ctime,c%masterVars%timef,c%masterVars%nsmax)

      end subroutine

      subroutine physical_GetRefinementCriteria(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical


          call Physical%GetRefinementCriteria(c%adaptiveVars%RefinerMarkel)
      end subroutine

      subroutine physical_PreRefine(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          if (c%adaptiveVars%kfl_AdaptiveRefinement == 1) then
              call Physical%PreRefine
          end if

      end subroutine

      subroutine physical_Refine(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          if (c%adaptiveVars%kfl_AdaptiveRefinement == 1) then
              call Physical%Refine('Refine')
          end if

      end subroutine

      subroutine physical_Rebalance(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          if (c%adaptiveVars%kfl_AdaptiveRefinement == 1) then
              call Physical%Refine('Rebalance')
          end if
      end subroutine

      subroutine physical_Begste(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical


          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%Begste
      end subroutine

      subroutine physical_MeshProjections(a,c,Physical,itask)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical
         integer(ip) :: itask

         if (a%ndela > c%masterVars%istep - 1) return

         call Physical%ProjectArraysOntoUndeformedMesh(c%domainVars%FMALEInterpolator,itask)

      end subroutine
      
      subroutine physical_MeshAdvections(a,c,Physical,itask)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical
         integer(ip) :: itask

         if (a%ndela > c%masterVars%istep - 1) return

         call Physical%AdvectArraysOntoUndeformedMesh(c%domainVars%FMALEAdvector,itask)

      end subroutine

      subroutine physical_DoIter(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%SetGlobalIteration(c%masterVars%iiter)
          call Physical%SetCouplingIteration(c%masterVars%cpiter)
          call Physical%Doiter

      end subroutine

      subroutine physical_Convergence(a,c,Physical,glres)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical
          real(rp) :: glres

          real(rp) :: residual

          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%GetResidual(glres)

      end subroutine

      subroutine physical_CoupledConvergence(a,c,Physical,kfl_flag,cpres)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical
          logical :: kfl_flag
          real(rp):: cpres

          kfl_flag = .true.
          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%GetCoupledConvFlag(kfl_flag)
          call Physical%GetCouplingResidual(cpres)

      end subroutine

      subroutine physical_Endste(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%Endste(c%masterVars%kfl_gotim)

      end subroutine

      subroutine physical_Output(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical

          if (a%ndela > c%masterVars%istep - 1) return

          call Physical%lev1Output

          if(a%kfl_preli==1.and.(mod(c%masterVars%istep,a%nprit)==0.or.((c%masterVars%ctime-c%masterVars%timef).gt.-1e-10).or.c%masterVars%istep>=c%masterVars%nsmax)) then
              Physical%oldnamda = Physical%namda
              call Physical%Restart(two)
          end if


          !GetData for statistics
          call Physical%GetTime(a%OutData%cputim)
          call Physical%GetMemo(a%OutData%Current,a%OutData%MaxMemo,a%OutData%TotalMemo,a%OutData%TotalMax)

      end subroutine

      subroutine physical_Turnof(a,c,Physical)
          class(DriverInterface) :: a
          type(caseVariables) :: c
          class(PhysicalProblem) :: Physical
          integer(8)  :: ToMem

          if (a%kfl_preli==1) then
            Physical%oldnamda = Physical%namda
            call Physical%Restart(two)
          endif

          call Physical%Turnof

          !GetData for statistics
          call Physical%GetTime(a%OutData%cputim)
          call Physical%GetMemo(a%OutData%Current,a%OutData%MaxMemo,a%OutData%TotalMemo,a%OutData%TotalMax)
      end subroutine

      subroutine physical_WriteTimes(a,Physical)
         implicit none
         class(DriverInterface) :: a
         class(PhysicalProblem) :: Physical

         call Physical%WriteTimes
      end subroutine

end module Mod_PhysicalProblemDriver
