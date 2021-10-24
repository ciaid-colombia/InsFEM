subroutine nsc_reampi(a)
   use typre
   use MPI
   use Mod_NSCompressible
   use Mod_BroadCastBuffer
   implicit none
   class(NSCompressibleProblem) :: a
   type(BroadCastBuffer) :: BBuffer

   !Initialize BroadcastBuffer
   call BBuffer%SetMPI(a%MPIcomm,a%MPIsize,a%MPIroot,a%MPIrank)
   call BBuffer%SetLOutputFile(a%lun_memo,a%lun_outpu)
   call BBuffer%Initialize(200,200)

   !Communicate nsc_reaphy

   call BBuffer%Add(a%lawde)
   call BBuffer%Add(a%kfl_outfm)

   call BBuffer%Add(a%kfl_advec)
   call BBuffer%Add(a%kfl_confi)
   call BBuffer%Add(a%kfl_visco)
   call BBuffer%Add(a%kfl_nsirstar)
   call BBuffer%Add(a%visco)
   call BBuffer%Add(a%tcond)
   call BBuffer%Add(a%cphea)
   call BBuffer%Add(a%cvhea)
   call BBuffer%Add(a%grnor)
   call BBuffer%Add(a%gravi)
   call BBuffer%Add(a%srce)
   call BBuffer%Add(a%lawdep)
   call BBuffer%Add(a%relpre)
   call BBuffer%Add(a%reltem)
   call BBuffer%Add(a%ndamp)
   call BBuffer%Add(a%dampco)
   call BBuffer%Add(a%dampxo)
   call BBuffer%Add(a%dampxf)
   call BBuffer%Add(a%dampff)
   call BBuffer%Add(a%rdamp)
   call BBuffer%Add(a%rdampco)
   call BBuffer%Add(a%rdampxo)
   call BBuffer%Add(a%rdampro)
   call BBuffer%Add(a%rdampex)
   call BBuffer%Add(a%rdampff)
   call BBuffer%Add(a%speed)
   call BBuffer%Add(a%molmass)
   call BBuffer%Add(a%densi)

   !Communicate nsc_reanut

   call BBuffer%Add(a%kfl_repro)
   call BBuffer%Add(a%kfl_shock)
   call BBuffer%Add(a%kfl_sctyp)
   call BBuffer%Add(a%kfl_trasg)
   call BBuffer%Add(a%kfl_nolsg)
   call BBuffer%Add(a%kfl_tacsg)
   call BBuffer%Add(a%kfl_stabm)
   call BBuffer%Add(a%kfl_jacgr)
   call BBuffer%Add(a%staco)
   call BBuffer%Add(a%shock)
   call BBuffer%Add(a%epspe)
   call BBuffer%Add(a%ErrorEstimatorTypeOfSubscales)

   !Communicate nsc_reaous
   !Forces and moments
   call BBuffer%Add(a%adimf)
   call BBuffer%Add(a%adimm)
   call BBuffer%Add(a%origm)

   call BBuffer%BroadCast
   call BBuffer%Dealloc
   
   if (a%kfl_advec == 0) a%kfl_jacgr = 0
   if (a%kfl_shock > 0) a%kfl_visco = 1
   if((a%grnor)>(zensi).or.(a%srce)>(zensi)) a%kfl_sourc = 1
   if(a%ndamp>0 .or. a%rdamp>0 .or.a%kfl_sourc==1) a%kfl_react = 1
        

   call a%SpecificNSCompReaMPI

end subroutine
