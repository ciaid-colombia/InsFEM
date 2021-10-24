module Mod_plcd_ComputeMean
   use typre
   use Mod_plcd_BaseElmope
   use MPI
   implicit none
   private
   public SetPointersComputeMean
   
   type, extends(PointerSetter) :: SPComputeMean
contains
      procedure :: SpecificSet => SpecificSetComputeMean
   end type
   type(SPComputeMean) :: SetPointersComputeMean
   
   Real(rp) :: tvol,gtvol
   
   
contains
   
   !----------------------------------------------------------------------------
   !Setting Pointers
  
   subroutine SpecificSetComputeMean(d)
      implicit none
      class(SPComputeMean) :: d
   
      !Post Process Strain
      if(a%ComputeMeanflag) then
         call ConcatenateProcedures(ProcHook%Initializations,zerotvol)
         call ConcatenateProcedures(ProcHook%Initializations,zeroMeanStress)
         call ConcatenateProcedures(ProcHook%InGaussElmats,Sumtvol)
         call ConcatenateProcedures(ProcHook%InGaussElmats,SumMeanStress)
         call ConcatenateProcedures(ProcHook%Finalizations,ComputeMeanStress)
         call ConcatenateProcedures(ProcHook%Finalizations,ComputeMeanConst)
      endif

   end subroutine
   
   subroutine zerotvol
      tvol = 0.0_rp
      gtvol = 0.0_rp
   
   end subroutine
   
   subroutine zeroMeanStress
      a%MeanStress = 0.0_rp
      a%GMeanStress = 0.0_rp
   
   end subroutine
     
   subroutine Sumtvol
      implicit none
      real(rp) :: weightfactor
      
      call e%GetWeightFactor(weightfactor)
      tvol = tvol + dvol*weightfactor

   end subroutine
      
   subroutine SumMeanStress
      implicit none
      real(rp) :: weightfactor
      
      call e%GetWeightFactor(weightfactor)
      a%MeanStress = a%MeanStress + dvol*stress*weightfactor

   end subroutine
   
   subroutine ComputeMeanStress
      implicit none
      integer(ip) :: ierr
      
      call  MPI_REDUCE(tvol,gtvol,1,MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
      call  MPI_REDUCE(a%MeanStress,a%GMeanStress,size(a%MeanStress),MPI_REAL8,MPI_SUM,a%MPIroot,a%MPIcomm,ierr)
      
      if (a%MPIrank == a%MPIroot) a%GMeanStress = a%GMeanStress/gtvol
   
   
 !     a%MeanStress = a%MeanStress/tvol

   end subroutine
   
   subroutine ComputeMeanConst
   
      if (a%MPIrank == a%MPIroot) a%MeanConst(:,a%CurrentStage) = a%GMeanStress

   end subroutine   
  
    
end module 
