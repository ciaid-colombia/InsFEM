subroutine ApplyErrorCriteria(MPIComm,nelem,RefinerErrorCriteria,RefinerErrorLimits,error,markel)
   use typre
   use MPI
   implicit none
   integer(ip)  :: nelem,MPIcomm
   character(9) :: RefinerErrorCriteria
   real(rp)     :: RefinerErrorLimits(2),error(nelem)
   integer(ip)  :: markel(nelem)
   
   integer(ip) :: ielem
   integer(ip), parameter :: numboxes = 10
   real(rp) :: lmaxerror,lminerror,maxerror,minerror,lerrorboxes(numboxes),errorboxes(numboxes)
   real(rp) :: ratio,errorinterval
   integer(ip) :: ibox,gnelem
   integer(ip) :: ObjectiveNumberOfElements,ierr
   real(rp) :: BetaParameter
   real(rp) :: NelemRatio,chi,epsilon1,epsilon2,ErrorValue1,ErrorValue2
   real(rp) :: MeanValue,Deviation,RefinementRatio,ErrorCenter
   
   integer(ip), save :: ipass = 0
   
   if (RefinerErrorCriteria == 'Tolerance') then
      do ielem = 1,nelem
         if (error(ielem) > RefinerErrorLimits(1)) then
            markel(ielem) = 1
         elseif (error(ielem) < RefinerErrorLimits(2)) then
            markel(ielem) = -1
         else
            markel(ielem) = 0
         endif
      enddo
   elseif (RefinerErrorCriteria == 'Elements') then
      !This refinement criteria follows Algorithm 3 in:
      !Badia Baiges "ADAPTIVE FINITE ELEMENT SIMULATION OF INCOMPRESSIBLE FLOWS BY
      !HYBRID CONTINUOUS-DISCONTINUOUS GALERKIN FORMULATIONS"
      
      call runend('This is not working')

      ObjectiveNumberOfElements = RefinerErrorLimits(1)
      BetaParameter = RefinerErrorLimits(2)
      
      !First, determine the maximum and minimum error amongst all processors
      lmaxerror = maxval(error)
      lminerror = minval(error)
      
      call  MPI_AllREDUCE( lmaxerror, maxerror, 1, MPI_REAL8, MPI_MAX,MPIcomm, ierr )
      call  MPI_AllREDUCE( lminerror, minerror, 1, MPI_REAL8, MPI_MIN,MPIcomm, ierr )
      
      !Secondly split that interval into 10 parts and count how many elemnts 
      !there are in each part
      !This allows us to build a piecewise linear model of the number of elements below a certain error threshold
      lerrorboxes(1:numboxes) = 0
      errorinterval = maxerror-minerror
      
      if (errorinterval == 0.0_rp) then
         markel = 0
         return
      endif
      
      do ielem = 1,nelem
         !Some elements are repeated in some processors, results will not be exactly the same when run with different number of processors
         ratio = (error(ielem)-minerror)/errorinterval
         ibox = min(numboxes,max(1,ceiling(ratio*numboxes)))
         lerrorboxes(ibox) = lerrorboxes(ibox)+1
         
      enddo   
      
      call  MPI_AllREDUCE( lerrorboxes, errorboxes, numboxes, MPI_REAL8, MPI_SUM,MPIcomm, ierr )
      gnelem = sum(errorboxes(1:numboxes))
      
      !Cumulative proportion
      errorboxes = errorboxes/gnelem
      do ibox = 2,numboxes
         errorboxes(ibox) = errorboxes(ibox-1)+errorboxes(ibox)
      enddo
      
!       NelemRatio = BetaParameter*(gnelem**2)/(ObjectiveNumberOfElements**2)
!       chi = 1.0_rp/(1.0_rp+NelemRatio)
!       epsilon1 = chi/2.0_rp
!       epsilon2 = (1+chi)/2.0_rp
!       call LookErrorValue(numboxes,1-epsilon1,maxerror,minerror,errorboxes,ErrorValue1)
!       call LookErrorValue(numboxes,epsilon2,maxerror,minerror,errorboxes,ErrorValue2)

      ipass = ipass +1
      if (ipass > 15) ObjectiveNumberOfElements = 1000

      BetaParameter = 0.5;
      RefinementRatio = 8;
      NelemRatio = 1.0_rp*(gnelem**2)/(ObjectiveNumberOfElements**2)
      MeanValue = NelemRatio/(NelemRatio+1.0_rp)
      Deviation = BetaParameter*0.5_rp
      epsilon1 = MeanValue+RefinementRatio*Deviation
      epsilon2 = MeanValue-Deviation
      call LookErrorValue(numboxes,MeanValue,maxerror,minerror,errorboxes,ErrorCenter)
      call LookErrorValue(numboxes,epsilon1,maxerror,minerror,errorboxes,ErrorValue1)
      call LookErrorValue(numboxes,epsilon2,maxerror,minerror,errorboxes,ErrorValue2)

      
      
      
      ErrorValue1 = max(ErrorValue1,ErrorCenter*2)


      
      
      
      do ielem = 1,nelem
         if (error(ielem) > ErrorValue1) then
            markel(ielem) = 1
         elseif (error(ielem) < ErrorValue2) then
            markel(ielem) = -1
         else
            markel(ielem) = 0
         endif
      enddo
      
   else
      call runend('Selected Error criteria not implemented')
   endif
   
 
   
end subroutine



subroutine LookErrorValue(numboxes,seekedproportion,maxerror,minerror,errorboxes,ErrorValue)
   use typre
   implicit none
   integer(ip) :: numboxes
   real(rp), intent(in) :: SeekedProportion
   real(rp) :: maxerror,minerror,errorboxes(*),Errorvalue

   real(rp) :: Errorinterval,ProportionLow,ProportionTop,PropotionInterval,error0,error1,errorinterval2
   integer(ip) :: ibox
   
   errorinterval = maxerror - minerror
   ProportionLow = 0.0_rp
   do ibox = 1,numboxes
      if (SeekedProportion < errorboxes(ibox)) then
         ProportionTop = ErrorBoxes(ibox)
         exit
      endif
      ProportionLow = ErrorBoxes(ibox)
   enddo
   PropotionInterval = ProportionTop-ProportionLow
   
   error0 = minerror + errorinterval/numboxes*(ibox-1)
   error1 = minerror + errorinterval/numboxes*ibox
   errorinterval2 = error1-error0
   
   ErrorValue = error0 + (SeekedProportion-ProportionLow)/PropotionInterval*errorinterval2
   
end subroutine