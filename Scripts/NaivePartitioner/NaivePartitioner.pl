#!/usr/bin/perl -w
use FileHandle;

# quit unless we have the correct number of command-line args
$num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: naivepartition.pl problem.gid numberofprocessors\n";
    exit;
}

my $start = time;

#Remove all backup files
system("find ./ -name '*~' | xargs rm");


# we got two command line args
# problem and number of processors
$problem=$ARGV[0];
$nproc=$ARGV[1];
 
$name = $problem;
$name =~ s/.gid//;   #remove .gid

#------------------------------------------------------------------------------------
#Directories and files
print "Start Directories and files \n";
@alldirs="$problem";
foreach $dir (@alldirs) {
  opendir(SDIR,$dir);
  @entries=readdir(SDIR);
  foreach (@entries) {
     $file = $dir."/".$_;
     if( -d $file ){
       if( $_ =~ /data/  ) { 
          opendir(DDIR,$file);
          @dataentries=readdir(DDIR);
	  $modulecount = 0;
          foreach (@dataentries) {
            if( ! ($_=~/^\./)  ) {
               $suff = $_;
               $suff =~ s/$name//;   #leave suffix
	       $modname = "$dir/data/$name$suff";
	       $rstname = "$dir/rst/$name$suff";
	            if( $rstname =~ /.rst/  ) {$rstfile = 1;}
               if( $modname =~ /.geo.dat/  ) {$geofile = $modname;}
               if($modname =~ /dom.dat/) { $domdatfile = $modname;}
               if( $modname =~ /dom.fix/  ) {$domfile = $modname;}
               if( ($modname =~ /.dat/) && !($modname =~ /.geo.dat/) && !($modname =~ /.fix/)  && !($modname =~ /.sol/)&& !($modname =~ /$name.dat/)  && !($modname =~ /.bod/)  && !($modname =~ /.dom/) ) {
		       @modwords = split m/\./, $modname;
		       $module = $modwords[2];
                       push(@savedmods,$module);}
	               $modulecount++;
            }
          }
       }
     }
  }
  closedir(SDIR);
}
print "End Directories and files \n";

print "Start DomDat \n";
#dom.dat reading
open(IN,"$domdatfile");
$reachedgeometry = 0;
foreach $stringin (<IN>) {
   $stringaux = $stringin;
   $stringaux =~ s/^\s+//; 
   @words = (split(/ {1,}/, $stringaux )) ;
   if (exists($words[0])) {
      if ($words[0]=~/NODAL_POINTS/i) {$nodalpoints = $words[1];}     
   }
   push(@saveddomdat,$stringin);
   #for non-root processes, we only want the geometry
   if ( $stringin =~ /GEOMETRY/  ) {$reachedgeometry = 1};
   if ($reachedgeometry == "1") {push(@saveddomdat1,$stringin)};

}
close(IN);
print "End DomDat \n";

print "Start ModDat \n";
#.mod.dat reading
$modco = 0;
foreach $module (@savedmods) {
   $reachedboundaryconditions = 0;
   open(IN,"$problem/data/$name.$module.dat");
   foreach $stringin (<IN>) {#read all lines in order
      push(@savedmoddat,$stringin);
      #for non root processes we want to copy only boundary conditions
      if ( $stringin =~ /BOUNDARY_CONDITIONS/ ) {$reachedboundaryconditions = 1};
      if ($reachedboundaryconditions == "1") {push(@savedmoddat1,$stringin)};
   }
   #Just in case we did not have boundary conditions
   if ($reachedboundaryconditions /= "1") {
      push(@savedmoddat1,"BOUNDARY_CONDITIONS \nEND_BOUNDARY_CONDITIONS");
   }   
   
   close(IN);
   $modulesdat[$modco]=[@savedmoddat];
   $modulesdat1[$modco]=[@savedmoddat1];
   undef(@savedmoddat);
   undef(@savedmoddat1);   
   $modco++;
}
print "End ModDat \n";

# modified problem directories and dat files
print "Start Modified Problem Directories and Data Files \n";
$problem_out = "$name.gid.$nproc";
system ("rm -r $problem_out");
system ("mkdir $problem_out");
system ("mkdir $problem_out/results");
system ("mkdir $problem_out/post");
system ("mkdir $problem_out/data");
system ("mkdir $problem_out/rst");
for ($iproc=0; $iproc<$nproc;$iproc++){
	$localdir = "$problem_out/data/data$iproc/";
   #$localdir = "$problem_out/data$iproc/";
        system ("mkdir $localdir");
        if ($iproc == 0){ 
            system ("cp $problem/data/$name.dat $localdir$name.dat");  
            system ("cp $problem/data/$name.dom.sol $localdir$name.dom.sol");
        }
	     $odomd = FileHandle->new(">$localdir$name.dom.dat");
        $printline = 1;
        if ($iproc == 0) {
           @auxsaveddomdat = @saveddomdat;
        } else {
           @auxsaveddomdat = @saveddomdat1;
        }
        foreach $stringin (@auxsaveddomdat) {
           $stringaux = $stringin;
           if($printline) {print $odomd "$stringaux";}
           $stringaux =~ s/^\s+//;
           @words = (split(/ {1,}/, $stringaux )) ;
           if (exists($words[0])) {
              if ($words[0] =~ /GEOMETRY,/i) {
		           print $odomd " INCLUDE  $name.geo.dat\n INCLUDE  $name.dom.fix\nEND_GEOMETRY\n\$-------------------------------------------------------------";
		           $printline = 0;
             }
          }
        }
        $odomd->close;
	$modco = 0;
	#Restart folders
	$rstdir = "$problem_out/rst/rst$iproc/";
   system ("mkdir $rstdir");
	foreach $module (@savedmods) {
           if ($iproc == 0){ system ("cp $problem/data/$name.$module.sol $localdir$name.$module.sol"); }
           
           if (defined ($rstfile)){ system ("cp $problem/rst/$name$iproc.$module.rst $rstdir$name$iproc.$module.rst 2>/dev/null"); }
           if (defined ($rstfile)){ system ("cp $problem/rst/$name$iproc.$module.rst2 $rstdir$name$iproc.$module.rst2 2>/dev/null"); }
           $omodd = FileHandle->new(">$localdir$name.$module.dat");
           if ($iproc == 0) {
      	    $savedmoddat = $modulesdat[$modco];
      	  } else {
      	    $savedmoddat = $modulesdat1[$modco];
      	  }
           foreach $stringin (@$savedmoddat) {
              $printline = 1;
              $stringaux = $stringin;
              $stringaux =~ s/^\s+//; 
              @words = (split(/ {1,}/, $stringaux )) ;
	      if (defined $words[0]){
                 if ($words[0]=~/INCLUDE/i) {
         	        print $omodd " INCLUDE  $name.$module.fix\n";
               	   $printline = 0;
                 }
              }
              if($printline) {print $omodd "$stringaux";}
           }
           $omodd->close;
	   $modco++;
	}
}
print "End Modified Problem Directories and Data Files \n";

$npoinlocal = int(($nodalpoints/$nproc)+0.99);

#geo.dat 
print "Start Geo \n";
open(IN,"$geofile");
my @elempart = (-1) x $nproc; 
$elempart[0]=-1;
#open files
for ($iproc=0; $iproc<$nproc;$iproc++){
   $localdir = "$problem_out/data/data$iproc/";
   #$localdir = "$problem_out/data$iproc/";
   $ogeod[$iproc] = FileHandle->new(">$localdir$name.geo.dat");
}
foreach $stringin (<IN>){ 
   $saveline = 1;
   $stringaux = $stringin;
   $stringaux =~ s/^\s+//; 
   @words = (split(/ {1,}/, $stringaux )) ;
   if ($words[0]=~/END_ELEMENTS/i) {$inelements = 0;$stringin = " -1 0 0\n$stringin"; }     
   if ($words[0]=~/END_COORDINATES/i) {$innodes = 0;}     
   if ($inelements) {
      if ($words[0]%100000==0) { 
         print "Element: $words[0] \n";
      }
      $pnode = $words[1]; 
      for($inode=1; $inode<$pnode+1; $inode++){
         $ipoin = $words[1+$inode];
         $iproc = int(($ipoin -1)/$npoinlocal);
         if ($elempart[$iproc]!=$words[0]){ 
            if (defined $ogeod[$iproc]) {
               $out = $ogeod[$iproc];
               print $out "$stringin";
           }
   	     $elempart[$iproc]=$words[0];
        }           
      }
	   $saveline =0;
   }    
   if ($innodes) {
	   $ipoin = $words[0];
	   if ($ipoin%100000==0) { 
         print "Point: $words[0] \n";
      }
      $iproc = int(($ipoin -1)/$npoinlocal); 
      if (defined $ogeod[$iproc]) {
         $out = $ogeod[$iproc];
         print $out "$stringin";
      }
	   $saveline =0;
   }     
   if (($words[0]=~/ELEMENTS/i)&&!($words[0]=~/END_ELEMENTS/i)) {$inelements = 1;}     
   if (($words[0]=~/COORDINATES/i)&&!($words[0]=~/END_COORDINATES/i)) {$innodes = 1;} 
   if ($saveline) {
      for ($iproc=0; $iproc<$nproc;$iproc++){
          if (defined $ogeod[$iproc]) {
             $out = $ogeod[$iproc];
             print $out "$stringin";
          }
      }  
   }
}
close(IN);
for ($iproc=0; $iproc<$nproc;$iproc++){
	$ogeod[$iproc]->close;
}
undef @ogeod;
print "End Geo \n";


#dom.fix 
print "Start Dom.fix \n";
#open files
for ($iproc=0; $iproc<$nproc;$iproc++){
	$localdir = "$problem_out/data/data$iproc/";
	#$localdir = "$problem_out/data$iproc/";
        $odomf[$iproc] = FileHandle->new(">$localdir$name.dom.fix");
	$bodypart[$iproc] = 0;
}
open(IN,"$domfile");
foreach $stringin (<IN>) {#read all lines in order
   $saveline = 1;
   $stringaux = $stringin;
   $stringaux =~ s/^\s+//; #remove leading spaces
   @words = (split(/ {1,}/, $stringaux )) ; # words in the line, space separator
   if ($words[0]=~/BODY_END/i) {$inbody = 0;}     
   if ($words[0]=~/END_EXTERNAL_NORMAL/i) {$innor = 0;}     
   if ($inbody) {
     if ($words[0]%100000==0) { 
         print "Body: $words[0] \n";
     }
     $pnode = $words[1]; 
     for($inode=1; $inode<$pnode+1; $inode++){
        $ipoin = $words[1+$inode];
        $iproc = int(($ipoin -1)/$npoinlocal);
        if ($bodypart[$iproc]!=$words[0]){ 
           if (defined $odomf[$iproc]) {
              $out = $odomf[$iproc];
              print $out "$stringin";
           }
           $bodypart[$iproc]=$words[0];
        }           
     }
	  $saveline =0;
   }    
   if ($innor) {
      if ($words[0]%100000==0) { 
         print "Normal: $words[0] \n";
      }
	   $ipoin = $words[0];
        $iproc = int(($ipoin -1)/$npoinlocal); 
        if (defined $odomf[$iproc]) {
           $out = $odomf[$iproc];
           print $out "$stringin";
        }
	   $saveline =0;
   }     
   if ($words[0]=~/BODY_DEFINITION/i) {$inbody = 1;}     
   if (($words[0]=~/EXTERNAL_NORMAL/i)&&!($words[0]=~/END_EXTERNAL_NORMAL/i)) {$innor = 1;} 
   if ($saveline) {
      for ($iproc=0; $iproc<$nproc;$iproc++){
          if (defined $odomf[$iproc]) {
             $out = $odomf[$iproc];
             print $out "$stringin";
          }
      }  
   }
}
close(IN);
for ($iproc=0; $iproc<$nproc;$iproc++){
	$odomf[$iproc]->close;
}
undef @odomf;
print "End Dom.fix \n";

#mod.fix
print "Start Mod.fix \n";
$modco = 0;
foreach $module (@savedmods) {
   #open files
   my @bounfixpart = (-1) x $nproc; 
   $bounfixpart[0]=-1;
   open(IN,"$problem/data/$name.$module.fix");
   for ($iproc=0; $iproc<$nproc;$iproc++){
	   $localdir = "$problem_out/data/data$iproc/";
	   #$localdir = "$problem_out/data$iproc/";
           $omodf[$iproc] = FileHandle->new(">$localdir$name.$module.fix");
   }
   $inbound = 0;
   $innodfix = 0;
   foreach $stringin (<IN>) {
      $saveline = 1;
      $stringaux = $stringin;
      $stringaux =~ s/^\s+//; #remove leading spaces
      @words = (split(/ {1,}/, $stringaux )) ; # words in the line, space separator
      if ($words[0]=~/END_ON_BOUNDARIES/i) {$inbound = 0;}     
      if ($words[0]=~/END_ON_NODES/i) {$innodfix = 0;}     
      if ($inbound) {
           if ($words[0]%100000==0) { 
               print "Mod $module, Boundary: $words[0] \n";
           }
           $pnode = $words[1]; 
           for($inode=1; $inode<$pnode+1; $inode++){
              $ipoin = $words[1+$inode];
              $iproc = int(($ipoin -1)/$npoinlocal);
              if ($bounfixpart[$iproc]/=$words[0]){ 
                 if (defined $omodf[$iproc]) {
                    $out = $omodf[$iproc];
                    print $out "$stringin";
                 }
   	         $bounfixpart[$iproc]=$words[0];
              }           
           }
   	   $saveline =0;
      }    
      if ($innodfix) {
         if ($words[0]%100000==0) { 
            print "Mod $module, Point: $words[0] \n";
         }
         $ipoin = $words[0];
         $iproc = int(($ipoin -1)/$npoinlocal); 
         if (defined $omodf[$iproc]) {
            $out = $omodf[$iproc];
            print $out "$stringin";
         }
         $saveline =0;
      }     
      if (($words[0]=~/ON_BOUNDARIES/i)&&!($words[0]=~/END_ON_BOUNDARIES/i)) {$inbound = 1;}     
      if (($words[0]=~/ON_NODES/i)&&!($words[0]=~/END_ON_NODES/i)) {$innodfix = 1;} 
      if ($saveline) {
         for ($iproc=0; $iproc<$nproc;$iproc++){
             if (defined $omodf[$iproc]) {
                $out = $omodf[$iproc];
                print $out "$stringin";
             }
         }  
      }  
   }
   #close files
   close(IN);
   for ($iproc=0; $iproc<$nproc;$iproc++){
      $omodf[$iproc]->close;
   }
   undef @omodf;
   $modco++;
}

print "End Mod.fix \n";

my $duration = time - $start;
print "Execution time: $duration s\n";
