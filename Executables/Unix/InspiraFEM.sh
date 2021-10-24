#!/bin/bash

#Put all arguments into the arguments array
declare -a arguments
count=0; 
for i in $*
do 
   arguments[$count]=$i;
   count=$((count+1));
done
narg=$count;

nprocs=1;
nomerge=1;
nomemcheck=0;
softlink=1;
runtype="";
steppost=1;

#Defaults
#Get Script directory name
pushd `dirname $0` > /dev/null
basedir=`pwd`
popd > /dev/null

#echo "$basedir";
casename="";
executable=$basedir"/InspiraFEM.O$COMPILER";
runtype="";
after="";
vtune=0;
inspe=0;
gprof="0"
#POSTPROCESSFOLDER="."
POSTPROCESSFOLDER="post"
RESTARTFOLDER="rst"
DATAFOLDER="data"
RESULTSFOLDER="results"
SCRIPT=$basedir"/../../Scripts"
SOURCE=$basedir"/../../Sources"


count=0;
casecount=0;
casenames="";
numcases=0;
while [ $count -lt $narg ];
do
   if [ $count -eq 0 ]; then
      casename=${arguments[$count]};
   elif [ ${arguments[$count]} = "-p" ]; then
      vtune=1
      executable=$basedir"/InspiraFEM.p$COMPILER"
      runtype="amplxe-runss";
      nomerge=1;  
   elif [ ${arguments[$count]} = "-i" ]; then
      inspe=1
      executable=$basedir"/InspiraFEM.p$COMPILER"
      runtype="inspector";
      nomerge=1;  
   elif [ ${arguments[$count]} = "-pg" ]; then
      executable=$basedir"/InspiraFEM.p$COMPILER"
      runtype=""
      gprof="1"
   elif [ ${arguments[$count]} = "-g" ]; then
      executable=$basedir"/InspiraFEM.g$COMPILER"
   elif [ ${arguments[$count]} = "-o" ]; then
      executable=$basedir"/InspiraFEM.O$COMPILER"
   elif [ ${arguments[$count]} = "-t" ]; then 
      executable=$basedir"/InspiraFEM.g$COMPILER"
      runtype="totalview"; 
      nomerge=1;
   elif [ ${arguments[$count]} = "-to" ]; then 
      executable=$basedir"/InspiraFEM.O$COMPILER"
      runtype="totalview"; 
      nomerge=1;   
   elif [ ${arguments[$count]} = "-v" ]; then 
      executable=$basedir"/InspiraFEM.g$COMPILER"
      runtype="valgrind"; 
      nomerge=1;
   elif [ ${arguments[$count]} = "-vo" ]; then 
      executable=$basedir"/InspiraFEM.O$COMPILER"
      runtype="valgrind"; 
      nomerge=1;
   elif [ ${arguments[$count]} = "-np" ]; then
      runtype="mpirun -np ${arguments[$((count+1))]}"  
      nprocs=${arguments[$((count+1))]}; 
   elif [ ${arguments[$count]} = "-nomemcheck" ]; then
      nomemcheck=1;
   elif [ ${arguments[$count]} = "-steppost" ]; then
      steppost=${arguments[$((count+1))]};    
   elif [ ${arguments[$count]} = "-info" ]; then
      after="$after -info"
   elif [ ${arguments[$count]} = "-nomerge" ]; then
      nomerge=1;
   elif [ ${arguments[$count]} = "-merge" ]; then
      nomerge=0;
   elif [ ${arguments[$count]} = "-mergeonly" ]; then
      runtype="mergeonly";   
      echo "Warning: mergeonly"
   elif [ ${arguments[$count]} = "-cases" ]; then
       after="$after -cases ${arguments[$count+1]}"
       numcases=${arguments[$((count+1))]};
       #softlink=0;
           #echo "numcases =$numcases"
           #echo "casecount =$casecount"
       while [ $casecount -lt $numcases ];
       do
           count=$((count+1));
           #echo "casecount =$casecount"
           #echo "${arguments[$count+1]}"
           after="$after ${arguments[$count+1]}"
           casecount=$((casecount+1));
       done
   elif [ ${arguments[$count]} = "-datafolder" ]; then
      after="$after -datafolder ${arguments[$count+1]}"
      $DATAFOLDER="${arguments[$count+1]}"
   elif [ ${arguments[$count]} = "-olddatafolder" ]; then
      after="$after -olddatafolder ${arguments[$count+1]}"
   elif [ ${arguments[$count]} = "-oldrestarfolder" ]; then
      after="$after -oldrestarfolder ${arguments[$count+1]}"
   elif [ ${arguments[$count]} = "-resultsfolder" ]; then
      $RESULTSFOLDER="${arguments[$count+1]}" 
      after="$after -resultsfolder ${arguments[$count+1]}"
   elif [ ${arguments[$count]} = "-restartfolder" ]; then
      #after="$after -restartfolder ${arguments[$count+1]}"
      RESTARTFOLDER="${arguments[$count+1]}"
   elif [ ${arguments[$count]} = "-postprocessfolder" ]; then
      #after="$after -postprocessfolder ${arguments[$count+1]}"
      POSTPROCESSFOLDER="${arguments[$count+1]}"
   elif [ ${arguments[$count]} = "-readtype" ]; then
      after="$after -readtype ${arguments[$count+1]}"   
   elif [ ${arguments[$count]} = "-nosoftlink" ]; then
      softlink=0;
   elif [ ${arguments[$count]} = "-multicomm" ]; then
      after="$after -multicomm ${arguments[$count+1]}" 
      nomerge=1;
   elif [ ${arguments[$count]} = "-memtrack" ]; then
      after="$after -memtrack" 
   fi
   count=$((count+1));
done

after="$after -postprocessfolder $POSTPROCESSFOLDER"
after="$after -restartfolder $RESTARTFOLDER"

#set runtype
if [ "$runtype" = "" ]; then
   if [ $nprocs -ne 1 ]; then
      runtype="mpirun -np $nprocs";
      
      if [ -f "$basedir/OMPI_hostfile" ]; then
         runtype="$runtype --hostfile $basedir/OMPI_hostfile";
      else
         echo "Warning: InspiraFEM was not able to find the OMPI hostfile"
      fi
      
   fi
elif [ "$runtype" = "totalview" ]; then
   runtype="totalview -mpi \"MPICH3\" -nodes $nprocs -np $nprocs";
   #runtype="totalview mpiexec -a -np $nprocs ";
   inbetween="-a";
elif [ "$runtype" = "valgrind" ]; then
   LD_PRELOAD=$prefix/lib/valgrind/libmpiwrap-amd64-linux.so  
   runtype="mpirun -np $nprocs $VALGRIND_DIR/bin/valgrind  --leak-check=full --undef-value-errors=yes --error-limit=no";
elif [ "$vtune" = "1" ]; then
   echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope;
   rm -rf r000
   runtype="amplxe-cl -collect hotspots -result-dir r000 mpirun -np $nprocs ";
elif [ "$inspe" = "1" ]; then
   echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope;
   rm -rf r001
   runtype="inspxe-cl -c mi2 -search-dir all:rp=$basedir -result-dir r001 mpirun -np $nprocs ";

fi

if [ "$runtype" != "mergeonly" ]; then
   #Remove the results and postprocess directory
   if [ -d ./post/ ]; then 
      find post/ ! -type d -exec rm '{}' \;
   fi
   if [ -d ./results/ ]; then    
      find results/ ! -type d -exec rm '{}' \;
   fi

   
   
   #Run the case
   runstring="$runtype $executable $inbetween $casename $after";
   echo "$runstring"
   
   eval $runstring;
fi

#Open amplxe-runss amplifier results
if [ "$vtune" = "1" ]; then
   eval amplxe-gui r000/r000.amplxe
fi

if [ "$inspe" = "1" ]; then
   eval inspxe-gui r001/r001.inspxe
fi

if [ "$gprof" = "1" ]; then
   gprof -a $executable gmon.out &>gmon.txt
fi

#Error code from the run
ERRORCODE=$?;

#Run the memory check only if the run was successful
if [ $ERRORCODE -eq 0 ]  && [ $nomemcheck -eq 0 ]; then
echo "Memchecking";
   pushd `dirname $0` > /dev/null
   basedir=`pwd`
   popd > /dev/null
   #basedir="$(dirname "$basedir")";
   basedir="$(dirname "$basedir")";
   basedir="$(dirname "$basedir")";
   $basedir/Scripts/Memcheck/Memchecker $casename;
fi
