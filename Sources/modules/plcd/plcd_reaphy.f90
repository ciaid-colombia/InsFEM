subroutine plcd_reaphy(a,itask)
   use typre
   use Mod_Listen
   use Mod_PLCD
   use Mod_plcd_Material
   use Mod_plcd_MaterialFactory
   use Mod_plcd_ReadMaterials
   implicit none
   class(PLCDProblem), target :: a
   integer(ip) :: itask

   !For Listener
   real(rp), pointer     :: param(:)
   character(5), pointer :: words(:)
   integer(ip), pointer  :: nnpar,nnwor

   integer(ip), save :: idnum, idnums, idnumss,ndime

   call a%Listener%getarrs(words,param,nnpar,nnwor)

   !Initializations
   if (itask == 0) then


   elseif (itask == 1) then
      if (words(1) == 'NUMBE' .and. words(3) == 'STAGE') then
         a%NumberOfStages = param(3)
         allocate(a%Stages(a%NumberOfStages))
         call a%Memor%allocObj(0,'Stages','plcd_memall',1*a%NumberOfStages)
         idnums = 0
      elseif (words(1) == 'BEGIN' .and. words(2) == 'STAGE') then
         idnums = idnums + 1

         do while(words(1)/='ENDST')
            call a%Listener%listen('plcd_reaphy')
            if (words(1) == 'NUMBE' .and. words(3) == 'SUBST') then
               a%Stages(idnums)%NumberOfSubstages = param(3)
               allocate(a%Stages(idnums)%Substages(a%Stages(idnums)%NumberOfSubstages))
               call a%Memor%allocObj(0,'SubStages','plcd_memall',1*a%Stages(idnums)%NumberOfSubstages)
               idnumss = 0
            elseif(words(1) == 'BEGIN' .and. words(2) == 'SUBST') then
               idnumss = idnumss + 1
               a%css => a%Stages(idnums)%Substages(idnumss)

               do while(words(1)/='ENDSU')
                  call a%Listener%listen('plcd_reaphy')
                  if(words(1) == 'INITI') then
                     a%css%InitialLoadFactor = param(1)
                  elseif(words(1) == 'FINAL') then
                      a%css%FinalLoadFactor = param(1)
                  elseif(words(1) == 'TIMES') then
                     a%css%TimeStep = param(1)
                  elseif(words(1) == 'TIMEI') then
                     a%css%TimeInterval = param(1)
                  elseif(words(1) == 'MAXIM') then
                     a%css%MaximumNonLinearIterations = int(param(1))
                  elseif(words(1) == 'ITERA') then
                     a%css%IterationsTolerance = param(1)
                  elseif(words(1) == 'NONLI') then
                     a%css%NonLinearAlgorithm = int(param(1))
                  endif
               enddo
            endif
         enddo
      elseif (words(1) == 'NEWMA') then
         a%kfl_TransientProblem = 1
         a%Beta = param(1)
         a%Gamma = param(2)
      elseif (words(1) == 'FORMU') then
         if (words(2) == 'UPDAT') then
            a%kfl_LargeStrains = 1
         elseif (words(2) == 'TOTAL') then
            a%kfl_LargeStrains = 2
         endif
      elseif (words(1) == 'TOPOL') then
         if (words(2) == 'SIMP ') then
            a%kfl_TopologyOptimization = 1
            call a%SIMPData%ReadData(a%Listener)
         elseif (words(2) == 'TOPOL') then
            a%kfl_TopologyOptimization = 2
            if (a%Listener%exists('STOCH')) then
               a%kfl_StochasticTopologyOptimization = 1
               call a%STOData%ReadData(a%Listener)
            endif
            call a%TDData%ReadData(a%Listener)
            
         endif
      elseif(words(1).eq.'GRAVI') then
         a%kfl_GravityForce = 1
         a%gravity(1) = a%Listener%getrea('GX   ',0.0_rp,'#x-component of g')
         a%gravity(2) = a%Listener%getrea('GY   ',0.0_rp,'#y-component of g')
         a%gravity(3) = a%Listener%getrea('GZ   ',0.0_rp,'#z-component of g')
      elseif(words(1).eq.'ANGVE') then
         if (words(2) == 'NEWMA') then
            a%kfl_RotatingFrame = 1
         elseif (words(2) == 'EXPLI') then
            a%kfl_RotatingFrame = 2
         else  
            call runend('There is no other time integrator scheme implemented for rotating frame of reference')
         endif
         a%angvelocity = param(2:4)
         a%angvelocitynorm2 = dot_product(a%angvelocity,a%angvelocity)
       endif

   elseif (itask == 2) then

      call a%Mesh%GetNdime(ndime)
      call ReadMaterials_Root(a%Listener,a%NumberOfMaterials,a%Materials,a%AuxRead_MaterialTypeList,a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank,ndime,a%kfl_LargeStrains,a%Memor)


   !Final operations
    elseif(itask==100) then
      !Time integration always on for plcd
      a%kfl_timei = 1_ip

    endif

end subroutine plcd_reaphy
