subroutine lmn_reampi(a)
   use MPI
   use typre
   use Mod_BroadCastBuffer
   use Mod_LowMach
   implicit none
   class(LowMachProblem) :: a
   type(BroadCastBuffer) :: BBuffer
   integer(ip) :: ierr

   call BBuffer%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call BBuffer%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call BBuffer%Initialize(100,500)
   
   !Communicate lmn_reaphy
   call BBuffer%add(a%kfl_eqnst)
   call BBuffer%add(a%kfl_advec)
   call BBuffer%add(a%kfl_confi)
   call BBuffer%add(a%kfl_visco)
   call BBuffer%add(a%kfl_sourc)
   
   call BBuffer%add(a%fvins)
   call BBuffer%add(a%densi)
   call BBuffer%add(a%itpre)
   call BBuffer%add(a%kfl_pther)
   call BBuffer%add(a%visco)
   call BBuffer%add(a%tcond)
   call BBuffer%add(a%cphea)
   call BBuffer%add(a%texpc)
   call BBuffer%add(a%sgasc)
   call BBuffer%add(a%grnor)
   call BBuffer%add(a%gravi)


   call BBuffer%add(a%react)
   call BBuffer%add(a%sourc)
   
   !Communicate lmn_reanut
   call BBuffer%add(a%kfl_repro)
   call BBuffer%add(a%kfl_trasg)
   call BBuffer%add(a%kfl_nolsg)
   call BBuffer%add(a%kfl_nolsgScheme)
   call BBuffer%add(a%mtrit)
   call BBuffer%add(a%kfl_tacsg)
   call BBuffer%add(a%kfl_stabm)
   call BBuffer%add(a%kfl_adapsgs)
        
   call BBuffer%add(a%staco)
   call BBuffer%add(a%tosgs)
   call BBuffer%add(a%relsg)
   call BBuffer%add(a%epspe)
  
   !Forces and moments
   call BBuffer%add(a%kfl_outfm)
   call BBuffer%add(a%adimf)
   call BBuffer%add(a%adimm)
   call BBuffer%add(a%adimh)
   call BBuffer%add(a%origm) 
   
   call BBuffer%BroadCast
   call BBuffer%Dealloc
   
   !For fixing CN
   if (a%kfl_tsche_1st_datafile == 'CNOBS') a%kfl_tsche_1st_datafile = 'CN   '

end subroutine
